"""
Configuration module for Essential Gene Analysis.

This module handles configuration loading and validation for the analysis pipeline.
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional
import yaml
import json


@dataclass
class AnalysisConfig:
    """Configuration for essential gene analysis."""
    
    # Model settings
    model_path: str
    model_format: str = 'sbml'  # 'sbml' or 'json'
    
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
    export_format: str = 'csv'  # 'csv', 'json', or 'both'
    generate_figures: bool = True
    
    # Column mappings for input files
    substrate_name_column: str = 'cmpd name'
    substrate_reaction_column: str = 'EX_ rxn ?'
    substrate_filter_column: str = 'finding'
    substrate_filter_value: str = 'growth'
    off_reactions_column: str = 'reactions'
    
    # Gap-filling settings
    gapfill_demand_reactions: bool = False
    gapfill_exchange_reactions: bool = True
    gapfill_integer_threshold: float = 1e-9
    
    def __post_init__(self):
        """Validate configuration after initialization."""
        self.output_dir = Path(self.output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        if self.model_path:
            model_path = Path(self.model_path)
            if not model_path.exists():
                raise FileNotFoundError(f"Model file not found: {self.model_path}")
    
    @classmethod
    def from_yaml(cls, yaml_path: str) -> 'AnalysisConfig':
        """Load configuration from a YAML file."""
        with open(yaml_path, 'r') as f:
            config_dict = yaml.safe_load(f)
        return cls(**config_dict)
    
    @classmethod
    def from_json(cls, json_path: str) -> 'AnalysisConfig':
        """Load configuration from a JSON file."""
        with open(json_path, 'r') as f:
            config_dict = json.load(f)
        return cls(**config_dict)
    
    def to_yaml(self, yaml_path: str) -> None:
        """Save configuration to a YAML file."""
        config_dict = self._to_dict()
        with open(yaml_path, 'w') as f:
            yaml.dump(config_dict, f, default_flow_style=False)
    
    def to_json(self, json_path: str) -> None:
        """Save configuration to a JSON file."""
        config_dict = self._to_dict()
        with open(json_path, 'w') as f:
            json.dump(config_dict, f, indent=2)
    
    def _to_dict(self) -> Dict:
        """Convert configuration to dictionary."""
        return {
            'model_path': str(self.model_path),
            'model_format': self.model_format,
            'solver': self.solver,
            'processes': self.processes,
            'carbon_uptake_rate': self.carbon_uptake_rate,
            'oxygen_uptake_rate': self.oxygen_uptake_rate,
            'growth_threshold': self.growth_threshold,
            'reference_substrate': self.reference_substrate,
            'carbon_substrates_file': self.carbon_substrates_file,
            'off_reactions_file': self.off_reactions_file,
            'universal_model_path': self.universal_model_path,
            'output_dir': str(self.output_dir),
            'export_format': self.export_format,
            'generate_figures': self.generate_figures,
            'substrate_name_column': self.substrate_name_column,
            'substrate_reaction_column': self.substrate_reaction_column,
            'substrate_filter_column': self.substrate_filter_column,
            'substrate_filter_value': self.substrate_filter_value,
            'off_reactions_column': self.off_reactions_column,
            'gapfill_demand_reactions': self.gapfill_demand_reactions,
            'gapfill_exchange_reactions': self.gapfill_exchange_reactions,
            'gapfill_integer_threshold': self.gapfill_integer_threshold,
        }


def create_default_config(output_path: str = 'config.yaml') -> AnalysisConfig:
    """
    Create a default configuration file.
    
    Args:
        output_path: Path for the configuration file
        
    Returns:
        The created AnalysisConfig object
    """
    config = AnalysisConfig(
        model_path='models/iML1515.xml',
        carbon_substrates_file='data/carbon_substrate.csv',
        off_reactions_file='data/off_reactions_aerobic.csv',
    )
    
    if output_path.endswith('.yaml') or output_path.endswith('.yml'):
        config.to_yaml(output_path)
    else:
        config.to_json(output_path)
    
    return config
