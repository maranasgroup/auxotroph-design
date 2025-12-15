"""
Essential Gene Analysis Package

A tool for identifying conditionally essential genes under different carbon
substrate conditions and finding potential rescue reactions through gap-filling.

Author: Roghaye Mohammadbeygi
"""

from .essential_genes import (
    EssentialGeneAnalyzer,
    CarbonSubstrate,
    GeneEssentialityResult,
    ConditionalEssentialityResult,
    create_bar_chart,
)
from .config import AnalysisConfig, create_default_config
from .visualization import (
    create_horizontal_bar_chart,
    create_vertical_bar_chart,
    create_heatmap,
    create_venn_diagram,
    create_upset_plot,
    generate_summary_statistics,
)

__version__ = '1.0.0'
__author__ = 'Roghaye Mohammadbeygi'
__all__ = [
    'EssentialGeneAnalyzer',
    'CarbonSubstrate',
    'GeneEssentialityResult',
    'ConditionalEssentialityResult',
    'AnalysisConfig',
    'create_default_config',
    'create_bar_chart',
    'create_horizontal_bar_chart',
    'create_vertical_bar_chart',
    'create_heatmap',
    'create_venn_diagram',
    'create_upset_plot',
    'generate_summary_statistics',
]
