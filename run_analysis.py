#!/usr/bin/env python3
"""
Main Script for Essential Gene Analysis

This script runs the complete workflow for identifying conditionally
essential genes and finding rescue reactions through gap-filling.

Usage:
    python run_analysis.py --config config/config.yaml
    python run_analysis.py  # Uses default config path
"""

import argparse
import json
import logging
import sys
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import Optional

import yaml

from src.essential_genes import EssentialGeneAnalyzer, CarbonSubstrate
from src.visualization import (
    create_horizontal_bar_chart,
    create_vertical_bar_chart,
    create_heatmap,
    generate_summary_statistics,
)


# =============================================================================
# Configuration
# =============================================================================

@dataclass
class Config:
    """Configuration for the analysis pipeline."""
    # Model settings
    model_path: str
    
    # Solver settings
    solver: str = 'glpk'
    processes: int = 4
    
    # Medium conditions
    carbon_uptake_rate: float = -10.0
    oxygen_uptake_rate: float = -20.0
    
    # Analysis settings
    growth_threshold: float = 1e-6
    reference_substrate: str = 'α-D-Glucose'
    
    # Input files
    carbon_substrates_file: Optional[str] = None
    off_reactions_file: Optional[str] = None
    universal_model_path: Optional[str] = None
    
    # Output settings
    output_dir: str = './outputs'
    export_format: str = 'both'  # 'csv', 'json', or 'both'
    generate_figures: bool = True
    
    # Column mappings for input files
    substrate_name_column: str = 'cmpd name'
    substrate_reaction_column: str = 'EX_ rxn ?'
    substrate_filter_column: str = 'finding'
    substrate_filter_value: str = 'growth'
    off_reactions_column: str = 'reactions'
    
    # Logging
    verbose: bool = False
    
    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'Config':
        """Load configuration from a YAML file."""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        
        # Handle None values and convert to proper types
        if config_dict is None:
            config_dict = {}
        
        return cls(**{k: v for k, v in config_dict.items() if v is not None or k == 'off_reactions_file'})
    
    @classmethod
    def from_json(cls, json_path: str) -> 'Config':
        """Load configuration from a JSON file."""
        with open(json_path, 'r') as f:
            config_dict = json.load(f)
        return cls(**config_dict)


def load_config(config_path: str) -> Config:
    """Load configuration from file (YAML or JSON)."""
    config_path = Path(config_path)
    
    if not config_path.exists():
        raise FileNotFoundError(f"Config file not found: {config_path}")
    
    if config_path.suffix in ('.yaml', '.yml'):
        return Config.from_yaml(str(config_path))
    elif config_path.suffix == '.json':
        return Config.from_json(str(config_path))
    else:
        raise ValueError(f"Unsupported config format: {config_path.suffix}")


# =============================================================================
# Logging Setup
# =============================================================================

def setup_logging(output_dir: Path, verbose: bool = False) -> logging.Logger:
    """Set up logging to both file and console."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / f"analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"
    
    # Clear any existing handlers
    root_logger = logging.getLogger()
    root_logger.handlers = []
    
    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler(sys.stdout)
        ]
    )
    return logging.getLogger(__name__)


# =============================================================================
# Main Analysis Pipeline
# =============================================================================

def run_analysis(config: Config):
    """
    Run the complete essential gene analysis pipeline.
    
    This workflow:
    1. Loads the metabolic model
    2. Identifies essential genes for each carbon substrate
    3. Finds conditionally essential genes (essential on specific substrates but not reference)
    4. Exports results and generates visualizations
    
    Args:
        config: Configuration object with all settings
    """
    # Setup output directory and logging
    output_dir = Path(config.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = setup_logging(output_dir, config.verbose)
    
    logger.info("=" * 60)
    logger.info("ESSENTIAL GENE ANALYSIS PIPELINE")
    logger.info("=" * 60)
    logger.info(f"Configuration loaded successfully")
    logger.info(f"Model: {config.model_path}")
    logger.info(f"Output directory: {output_dir.absolute()}")
    
    # =========================================================================
    # Step 1: Initialize the analyzer
    # =========================================================================
    logger.info("\n[Step 1] Loading metabolic model...")
    
    try:
        analyzer = EssentialGeneAnalyzer(
            model_path=config.model_path,
            solver=config.solver,
            processes=config.processes
        )
    except FileNotFoundError as e:
        logger.error(f"Model file not found: {config.model_path}")
        logger.error("Please check the 'model_path' in your config file.")
        raise
    
    # =========================================================================
    # Step 2: Load carbon substrates
    # =========================================================================
    logger.info("\n[Step 2] Loading carbon substrates...")
    
    if not config.carbon_substrates_file:
        logger.error("No carbon substrates file specified in config!")
        raise ValueError("carbon_substrates_file is required in config")
    
    try:
        analyzer.load_carbon_substrates(
            config.carbon_substrates_file,
            name_column=config.substrate_name_column,
            reaction_column=config.substrate_reaction_column,
            filter_column=config.substrate_filter_column,
            filter_value=config.substrate_filter_value
        )
    except FileNotFoundError:
        logger.error(f"Carbon substrates file not found: {config.carbon_substrates_file}")
        raise
    
    # =========================================================================
    # Step 3: Load reactions to turn off (optional)
    # =========================================================================
    if config.off_reactions_file:
        logger.info("\n[Step 3] Loading off reactions...")
        try:
            analyzer.load_off_reactions(
                config.off_reactions_file,
                column=config.off_reactions_column
            )
        except FileNotFoundError:
            logger.warning(f"Off reactions file not found: {config.off_reactions_file}")
            logger.warning("Continuing without turning off reactions...")
    else:
        logger.info("\n[Step 3] No off reactions file specified, skipping...")
    
    # =========================================================================
    # Step 4: Find essential genes for all substrates
    # =========================================================================
    logger.info("\n[Step 4] Finding essential genes for each carbon substrate...")
    logger.info(f"This may take a while depending on the number of substrates...")
    
    def progress_callback(current, total):
        logger.info(f"Progress: {current}/{total} substrates analyzed")
    
    results = analyzer.analyze_all_substrates(progress_callback=progress_callback)
    
    # =========================================================================
    # Step 5: Export raw results
    # =========================================================================
    logger.info("\n[Step 5] Exporting results...")
    
    if config.export_format in ('csv', 'both'):
        csv_path = output_dir / 'essential_genes_matrix.csv'
        analyzer.export_results_to_csv(results, csv_path)
        logger.info(f"  Saved: {csv_path}")
    
    if config.export_format in ('json', 'both'):
        json_path = output_dir / 'essential_genes.json'
        analyzer.export_results_to_json(results, json_path)
        logger.info(f"  Saved: {json_path}")
    
    # =========================================================================
    # Step 6: Find conditionally essential genes
    # =========================================================================
    logger.info(f"\n[Step 6] Finding conditionally essential genes...")
    logger.info(f"  Reference substrate: {config.reference_substrate}")
    
    try:
        cond_essential = analyzer.find_conditionally_essential_genes(
            reference_substrate=config.reference_substrate
        )
    except ValueError as e:
        logger.error(f"Reference substrate '{config.reference_substrate}' not found in results.")
        logger.error("Available substrates: " + ", ".join(results.keys()))
        raise
    
    # Export conditionally essential genes
    cond_path = output_dir / 'conditionally_essential.json'
    cond_export = {k: list(v) for k, v in cond_essential.items()}
    with open(cond_path, 'w') as f:
        json.dump(cond_export, f, indent=2)
    logger.info(f"  Saved: {cond_path}")
    
    # Create detailed report
    cond_report_path = output_dir / 'conditionally_essential_detailed.txt'
    with open(cond_report_path, 'w') as f:
        f.write("CONDITIONALLY ESSENTIAL GENES REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Reference substrate: {config.reference_substrate}\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for substrate, genes in sorted(cond_essential.items()):
            if genes:
                f.write(f"\n{substrate} ({len(genes)} genes):\n")
                f.write("-" * 40 + "\n")
                for gene in sorted(genes):
                    try:
                        gene_obj = analyzer.model.genes.get_by_id(gene)
                        gene_name = gene_obj.name if gene_obj.name else gene
                    except:
                        gene_name = gene
                    f.write(f"  {gene}: {gene_name}\n")
    logger.info(f"  Saved: {cond_report_path}")
    
    # =========================================================================
    # Step 7: Generate visualizations
    # =========================================================================
    if config.generate_figures:
        logger.info("\n[Step 7] Generating visualizations...")
        
        # Count conditionally essential genes per substrate
        cond_counts = {k: len(v) for k, v in cond_essential.items() if v}
        
        if cond_counts:
            try:
                # Horizontal bar chart
                create_horizontal_bar_chart(
                    cond_counts,
                    output_dir / 'conditionally_essential_horizontal.png',
                    title='Conditionally Essential Genes by Carbon Substrate',
                    xlabel='Number of Conditionally Essential Genes'
                )
                logger.info(f"  Saved: conditionally_essential_horizontal.png")
                
                # Vertical bar chart
                create_vertical_bar_chart(
                    cond_counts,
                    output_dir / 'conditionally_essential_vertical.png',
                    title='Conditionally Essential Genes by Carbon Substrate',
                    ylabel='Number of Conditionally Essential Genes'
                )
                logger.info(f"  Saved: conditionally_essential_vertical.png")
            except Exception as e:
                logger.warning(f"Could not create conditionally essential plots: {e}")
        
        # Essential genes count per substrate
        essential_counts = {s: len(r.essential_genes) for s, r in results.items()}
        try:
            create_horizontal_bar_chart(
                essential_counts,
                output_dir / 'essential_genes_count.png',
                title='Essential Genes by Carbon Substrate',
                xlabel='Number of Essential Genes',
                color='#3498DB'
            )
            logger.info(f"  Saved: essential_genes_count.png")
        except Exception as e:
            logger.warning(f"Could not create essential genes plot: {e}")
        
        # Comparison matrix
        matrix = analyzer.generate_comparison_matrix(results)
        matrix_path = output_dir / 'gene_essentiality_matrix.csv'
        matrix.to_csv(matrix_path)
        logger.info(f"  Saved: {matrix_path}")
        
        # Heatmap (if not too large)
        if len(matrix) <= 200:
            try:
                create_heatmap(
                    matrix,
                    output_dir / 'gene_essentiality_heatmap.png',
                    title='Gene Essentiality Across Carbon Substrates'
                )
                logger.info(f"  Saved: gene_essentiality_heatmap.png")
            except Exception as e:
                logger.warning(f"Could not create heatmap: {e}")
        else:
            logger.info(f"  Skipping heatmap (too many genes: {len(matrix)})")
        
        # Summary statistics
        try:
            summary_path = output_dir / 'summary_statistics.csv'
            generate_summary_statistics(
                {s: r.essential_genes for s, r in results.items()},
                summary_path
            )
            logger.info(f"  Saved: {summary_path}")
        except Exception as e:
            logger.warning(f"Could not generate summary statistics: {e}")
    else:
        logger.info("\n[Step 7] Skipping visualizations (disabled in config)")
    
    # =========================================================================
    # Step 8: Generate final report
    # =========================================================================
    logger.info("\n[Step 8] Generating final report...")
    
    cond_counts = {k: len(v) for k, v in cond_essential.items() if v}
    
    report_path = output_dir / 'analysis_report.txt'
    with open(report_path, 'w') as f:
        f.write("ESSENTIAL GENE ANALYSIS REPORT\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Model: {config.model_path}\n")
        f.write(f"Analysis date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Solver: {config.solver}\n")
        f.write(f"Reference substrate: {config.reference_substrate}\n\n")
        
        f.write("SUMMARY\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total substrates analyzed: {len(results)}\n")
        
        all_genes = set()
        for r in results.values():
            all_genes.update(r.essential_genes)
        f.write(f"Total unique essential genes: {len(all_genes)}\n")
        f.write(f"Substrates with conditionally essential genes: {len(cond_counts)}\n\n")
        
        f.write("ESSENTIAL GENES BY SUBSTRATE\n")
        f.write("-" * 40 + "\n")
        for substrate, result in sorted(results.items(), key=lambda x: -len(x[1].essential_genes)):
            f.write(f"{substrate}: {len(result.essential_genes)} genes (growth: {result.growth_rate:.4f})\n")
        
        if cond_counts:
            f.write("\n\nCONDITIONALLY ESSENTIAL GENES SUMMARY\n")
            f.write("-" * 40 + "\n")
            for substrate, count in sorted(cond_counts.items(), key=lambda x: -x[1]):
                f.write(f"{substrate}: {count} genes\n")
    
    logger.info(f"  Saved: {report_path}")
    
    # =========================================================================
    # Done
    # =========================================================================
    logger.info("\n" + "=" * 60)
    logger.info("ANALYSIS COMPLETE")
    logger.info("=" * 60)
    logger.info(f"\nOutput files saved to: {output_dir.absolute()}")
    logger.info("\nGenerated files:")
    for f in sorted(output_dir.glob('*')):
        logger.info(f"  - {f.name}")
    
    return results, cond_essential


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    """Main entry point with argument parsing."""
    parser = argparse.ArgumentParser(
        description='Essential Gene Analysis Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
    # Run with config file
    python run_analysis.py --config config/config.yaml
    
    # Run with default config location
    python run_analysis.py
    
    # Run with verbose output
    python run_analysis.py --config config/config.yaml --verbose
        '''
    )
    
    parser.add_argument(
        '--config', '-c',
        default='config/config.yaml',
        help='Path to configuration file (YAML or JSON). Default: config/config.yaml'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output (overrides config setting)'
    )
    
    args = parser.parse_args()
    
    # Load configuration
    try:
        config = load_config(args.config)
        print(f"Loaded configuration from: {args.config}")
    except FileNotFoundError:
        print(f"ERROR: Config file not found: {args.config}")
        print("\nPlease create a config file or specify a valid path.")
        print("Example config file location: config/config.yaml")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to load config: {e}")
        sys.exit(1)
    
    # Override verbose if specified on command line
    if args.verbose:
        config.verbose = True
    
    # Run analysis
    try:
        run_analysis(config)
    except FileNotFoundError as e:
        print(f"\nERROR: {e}")
        print("Please check the file paths in your configuration.")
        sys.exit(1)
    except Exception as e:
        print(f"\nERROR: Analysis failed: {e}")
        if config.verbose:
            import traceback
            traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
