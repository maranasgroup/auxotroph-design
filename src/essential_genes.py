"""
Essential Gene Analysis Module

This module provides functionality for identifying conditionally essential genes
under different carbon substrate conditions and finding potential rescue reactions
through gap-filling analysis.

Author: Roghaye Mohammadbeygi
Date: 2025
"""

import logging
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple, Union
from dataclasses import dataclass, field

import cobra
from cobra.io import read_sbml_model, load_json_model
from cobra.flux_analysis import find_essential_genes, find_essential_reactions
from cobra.flux_analysis.gapfilling import GapFiller
import pandas as pd
import numpy as np

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


@dataclass
class CarbonSubstrate:
    """Represents a carbon substrate with its exchange reaction."""
    name: str
    exchange_reaction: str


@dataclass
class GeneEssentialityResult:
    """Results from gene essentiality analysis for a single condition."""
    substrate: str
    essential_genes: List[str]
    growth_rate: float


@dataclass
class ConditionalEssentialityResult:
    """Results from conditional essentiality analysis."""
    gene_id: str
    gene_name: str
    substrate: str
    rescue_reactions: List[str]
    
    
class EssentialGeneAnalyzer:
    """
    Main class for analyzing essential genes under different carbon conditions.
    
    This class provides methods to:
    1. Load and configure metabolic models
    2. Identify essential genes for different carbon substrates
    3. Find conditionally essential genes
    4. Perform gap-filling to identify rescue reactions
    
    Attributes:
        model: The COBRApy metabolic model
        carbon_substrates: List of carbon substrates to analyze
        solver: The solver to use (default: 'cplex' or 'glpk')
    """
    
    def __init__(
        self,
        model_path: Union[str, Path],
        solver: str = 'glpk',
        processes: int = 4
    ):
        """
        Initialize the analyzer with a metabolic model.
        
        Args:
            model_path: Path to the SBML or JSON model file
            solver: Solver to use ('cplex', 'glpk', 'gurobi')
            processes: Number of parallel processes for analysis
        """
        self.model_path = Path(model_path)
        self.solver = solver
        self.processes = processes
        self.model = self._load_model()
        self.carbon_substrates: List[CarbonSubstrate] = []
        self.off_reactions: List[str] = []
        self._essential_genes_cache: Dict[str, Set[str]] = {}
        
        logger.info(f"Loaded model: {self.model.id}")
        logger.info(f"Reactions: {len(self.model.reactions)}, "
                   f"Metabolites: {len(self.model.metabolites)}, "
                   f"Genes: {len(self.model.genes)}")
    
    def _load_model(self) -> cobra.Model:
        """Load the metabolic model from file."""
        if self.model_path.suffix == '.xml':
            model = read_sbml_model(str(self.model_path))
        elif self.model_path.suffix == '.json':
            model = load_json_model(str(self.model_path))
        else:
            raise ValueError(f"Unsupported model format: {self.model_path.suffix}")
        
        model.solver = self.solver
        return model
    
    def load_carbon_substrates(
        self,
        csv_path: Union[str, Path],
        name_column: str = 'cmpd name',
        reaction_column: str = 'EX_ rxn ?',
        filter_column: Optional[str] = 'finding',
        filter_value: str = 'growth'
    ) -> None:
        """
        Load carbon substrates from a CSV file.
        
        Args:
            csv_path: Path to CSV file containing carbon substrate information
            name_column: Column name for substrate names
            reaction_column: Column name for exchange reaction IDs
            filter_column: Optional column to filter substrates
            filter_value: Value to filter by in filter_column
        """
        df = pd.read_csv(csv_path)
        
        if filter_column and filter_column in df.columns:
            df = df[df[filter_column] == filter_value]
        
        self.carbon_substrates = []
        for _, row in df.iterrows():
            rxn_name = self._normalize_reaction_name(row[reaction_column])
            substrate = CarbonSubstrate(
                name=row[name_column],
                exchange_reaction=rxn_name
            )
            self.carbon_substrates.append(substrate)
        
        logger.info(f"Loaded {len(self.carbon_substrates)} carbon substrates")
    
    def load_off_reactions(
        self,
        csv_path: Union[str, Path],
        column: str = 'reactions'
    ) -> None:
        """
        Load reactions to turn off from a CSV file.
        
        Args:
            csv_path: Path to CSV file containing reaction IDs to turn off
            column: Column name containing reaction IDs
        """
        df = pd.read_csv(csv_path)
        self.off_reactions = list(df[column])
        logger.info(f"Loaded {len(self.off_reactions)} reactions to turn off")
    
    @staticmethod
    def _normalize_reaction_name(rxn_name: str) -> str:
        """Normalize reaction names for consistency between models."""
        rxn_name = (rxn_name
                   .replace("(e)", "_e")
                   .replace("_D", "__D")
                   .replace("_L", "__L")
                   .replace("EX_glc_e", "EX_glc__D_e"))
        return rxn_name
    
    def _apply_off_reactions(self, model: cobra.Model) -> List[str]:
        """
        Turn off specified reactions in the model.
        
        Returns:
            List of reactions that couldn't be found in the model
        """
        not_found = []
        for rxn_id in self.off_reactions:
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            except KeyError:
                not_found.append(rxn_id)
        
        if not_found:
            logger.warning(f"Reactions not found in model: {not_found}")
        
        return not_found
    
    def _set_carbon_conditions(
        self,
        model: cobra.Model,
        substrate: CarbonSubstrate,
        carbon_uptake: float = -10.0,
        oxygen_uptake: float = -20.0
    ) -> None:
        """
        Set medium conditions for a specific carbon substrate.
        
        Args:
            model: The model to modify
            substrate: The carbon substrate to use
            carbon_uptake: Lower bound for carbon uptake (negative for uptake)
            oxygen_uptake: Lower bound for oxygen uptake (negative for uptake)
        """
        # Turn off all carbon substrates
        for sub in self.carbon_substrates:
            try:
                rxn = model.reactions.get_by_id(sub.exchange_reaction)
                rxn.lower_bound = 0
                rxn.upper_bound = 0
            except KeyError:
                pass
        
        # Turn on the specific carbon substrate
        try:
            rxn = model.reactions.get_by_id(substrate.exchange_reaction)
            rxn.lower_bound = carbon_uptake
        except KeyError:
            logger.error(f"Exchange reaction not found: {substrate.exchange_reaction}")
            raise
        
        # Set oxygen conditions
        try:
            model.reactions.get_by_id("EX_o2_e").lower_bound = oxygen_uptake
        except KeyError:
            logger.warning("Oxygen exchange reaction (EX_o2_e) not found")
    
    def find_essential_genes_for_substrate(
        self,
        substrate: CarbonSubstrate,
        threshold: float = 1e-6
    ) -> GeneEssentialityResult:
        """
        Find essential genes for a specific carbon substrate.
        
        Args:
            substrate: The carbon substrate to analyze
            threshold: Growth threshold below which a gene is considered essential
            
        Returns:
            GeneEssentialityResult with essential genes and growth rate
        """
        with self.model.copy() as model:
            self._set_carbon_conditions(model, substrate)
            
            # Optimize to get baseline growth
            solution = model.optimize()
            growth_rate = solution.objective_value if solution.status == 'optimal' else 0.0
            
            if growth_rate < threshold:
                logger.warning(f"No growth on {substrate.name} (rate: {growth_rate})")
                return GeneEssentialityResult(
                    substrate=substrate.name,
                    essential_genes=[],
                    growth_rate=growth_rate
                )
            
            # Find essential genes
            essential = find_essential_genes(model, processes=self.processes)
            essential_ids = [gene.id for gene in essential]
            
            logger.info(f"{substrate.name}: {len(essential_ids)} essential genes, "
                       f"growth rate: {growth_rate:.4f}")
            
            return GeneEssentialityResult(
                substrate=substrate.name,
                essential_genes=essential_ids,
                growth_rate=growth_rate
            )
    
    def analyze_all_substrates(
        self,
        progress_callback: Optional[callable] = None
    ) -> Dict[str, GeneEssentialityResult]:
        """
        Find essential genes for all loaded carbon substrates.
        
        Args:
            progress_callback: Optional callback function for progress updates
            
        Returns:
            Dictionary mapping substrate names to essentiality results
        """
        results = {}
        
        for i, substrate in enumerate(self.carbon_substrates):
            logger.info(f"Analyzing {substrate.name} ({i+1}/{len(self.carbon_substrates)})")
            
            try:
                result = self.find_essential_genes_for_substrate(substrate)
                results[substrate.name] = result
                self._essential_genes_cache[substrate.name] = set(result.essential_genes)
            except Exception as e:
                logger.error(f"Error analyzing {substrate.name}: {e}")
                continue
            
            if progress_callback:
                progress_callback(i + 1, len(self.carbon_substrates))
        
        return results
    
    def find_conditionally_essential_genes(
        self,
        reference_substrate: str = 'α-D-Glucose'
    ) -> Dict[str, Set[str]]:
        """
        Find genes that are essential only under specific carbon conditions.
        
        Conditionally essential genes are those that are essential for growth
        on a specific substrate but not on the reference substrate (usually glucose).
        
        Args:
            reference_substrate: The reference substrate for comparison
            
        Returns:
            Dictionary mapping substrate names to sets of conditionally essential genes
        """
        if reference_substrate not in self._essential_genes_cache:
            raise ValueError(f"Reference substrate '{reference_substrate}' not analyzed")
        
        reference_genes = self._essential_genes_cache[reference_substrate]
        conditionally_essential = {}
        
        for substrate, genes in self._essential_genes_cache.items():
            if substrate == reference_substrate:
                conditionally_essential[substrate] = set()
                continue
            
            # Genes essential for this substrate but not for reference
            cond_genes = genes - reference_genes
            conditionally_essential[substrate] = cond_genes
            
            if cond_genes:
                logger.info(f"{substrate}: {len(cond_genes)} conditionally essential genes")
        
        return conditionally_essential
    
    def identify_rescue_reactions(
        self,
        gene_id: str,
        substrate: CarbonSubstrate,
        universal_model: Optional[cobra.Model] = None
    ) -> List[str]:
        """
        Identify reactions that can rescue growth after gene knockout.
        
        Uses gap-filling to find reactions from a universal model that can
        restore growth when a gene is knocked out.
        
        Args:
            gene_id: ID of the gene to knock out
            substrate: Carbon substrate for the growth condition
            universal_model: Universal model for gap-filling (optional)
            
        Returns:
            List of reaction IDs that can rescue growth
        """
        with self.model.copy() as model:
            self._set_carbon_conditions(model, substrate)
            
            # Knock out the gene
            try:
                gene = model.genes.get_by_id(gene_id)
                gene.knock_out()
            except KeyError:
                logger.error(f"Gene not found: {gene_id}")
                return []
            
            # Check if knockout affects growth
            solution = model.optimize()
            if solution.status == 'optimal' and solution.objective_value > 1e-6:
                logger.info(f"Gene {gene_id} knockout doesn't prevent growth")
                return []
            
            # Get reactions associated with the gene
            associated_reactions = []
            for rxn in model.reactions:
                if gene_id in [g.id for g in rxn.genes]:
                    associated_reactions.append(rxn.id)
            
            return associated_reactions
    
    def run_gapfilling(
        self,
        gene_id: str,
        substrate: CarbonSubstrate,
        universal_model: cobra.Model,
        demand_reactions: bool = False,
        exchange_reactions: bool = True
    ) -> List[cobra.Reaction]:
        """
        Run gap-filling to find rescue reactions for a gene knockout.
        
        Args:
            gene_id: ID of the gene to knock out
            substrate: Carbon substrate for the growth condition
            universal_model: Universal model containing potential rescue reactions
            demand_reactions: Whether to include demand reactions
            exchange_reactions: Whether to include exchange reactions
            
        Returns:
            List of reactions that can rescue growth
        """
        with self.model.copy() as model:
            self._set_carbon_conditions(model, substrate)
            
            # Knock out the gene
            try:
                gene = model.genes.get_by_id(gene_id)
                gene.knock_out()
            except KeyError:
                logger.error(f"Gene not found: {gene_id}")
                return []
            
            try:
                gapfiller = GapFiller(
                    model,
                    universal=universal_model,
                    demand_reactions=demand_reactions,
                    exchange_reactions=exchange_reactions,
                    integer_threshold=1e-9
                )
                gapfiller.model.solver.configuration.tolerances.feasibility = 1e-9
                gapfiller.model.solver.configuration.tolerances.integrality = 1e-9
                
                solutions = gapfiller.fill()
                
                if solutions:
                    return solutions[0]
                return []
                
            except Exception as e:
                logger.error(f"Gap-filling failed for gene {gene_id}: {e}")
                return []
    
    def export_results_to_csv(
        self,
        results: Dict[str, GeneEssentialityResult],
        output_path: Union[str, Path]
    ) -> None:
        """
        Export essential gene results to a CSV file.
        
        Args:
            results: Dictionary of essentiality results
            output_path: Path for the output CSV file
        """
        # Find maximum number of essential genes
        max_genes = max(len(r.essential_genes) for r in results.values())
        
        # Create DataFrame
        data = {}
        for substrate, result in results.items():
            genes = result.essential_genes + [''] * (max_genes - len(result.essential_genes))
            data[substrate] = genes
        
        df = pd.DataFrame(data)
        df.to_csv(output_path, index=False)
        logger.info(f"Results exported to {output_path}")
    
    def export_results_to_json(
        self,
        results: Dict[str, GeneEssentialityResult],
        output_path: Union[str, Path]
    ) -> None:
        """
        Export essential gene results to a JSON file.
        
        Args:
            results: Dictionary of essentiality results
            output_path: Path for the output JSON file
        """
        import json
        
        export_data = {
            substrate: {
                'essential_genes': result.essential_genes,
                'growth_rate': result.growth_rate
            }
            for substrate, result in results.items()
        }
        
        with open(output_path, 'w') as f:
            json.dump(export_data, f, indent=2)
        
        logger.info(f"Results exported to {output_path}")
    
    def generate_comparison_matrix(
        self,
        results: Dict[str, GeneEssentialityResult]
    ) -> pd.DataFrame:
        """
        Generate a matrix comparing gene essentiality across substrates.
        
        Args:
            results: Dictionary of essentiality results
            
        Returns:
            DataFrame with genes as rows and substrates as columns
        """
        # Collect all unique genes
        all_genes = set()
        for result in results.values():
            all_genes.update(result.essential_genes)
        
        # Create comparison matrix
        matrix_data = []
        for gene in sorted(all_genes):
            row = {'gene': gene}
            for substrate, result in results.items():
                row[substrate] = 1 if gene in result.essential_genes else 0
            matrix_data.append(row)
        
        return pd.DataFrame(matrix_data).set_index('gene')


def create_bar_chart(
    data: Dict[str, int],
    output_path: Union[str, Path],
    title: str = "Conditionally Essential Genes",
    xlabel: str = "Number of Genes",
    figsize: Tuple[float, float] = (10, 8)
) -> None:
    """
    Create a horizontal bar chart of conditionally essential genes.
    
    Args:
        data: Dictionary mapping substrate names to gene counts
        output_path: Path for the output figure
        title: Chart title
        xlabel: X-axis label
        figsize: Figure size (width, height)
    """
    import matplotlib.pyplot as plt
    
    # Sort by values
    sorted_data = dict(sorted(data.items(), key=lambda x: x[1], reverse=True))
    
    fig, ax = plt.subplots(figsize=figsize)
    
    y_pos = np.arange(len(sorted_data))
    ax.barh(y_pos, list(sorted_data.values()), align='center')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(list(sorted_data.keys()))
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    logger.info(f"Figure saved to {output_path}")
