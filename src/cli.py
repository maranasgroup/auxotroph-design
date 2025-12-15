#!/usr/bin/env python3
"""
Command-line interface for Essential Gene Analysis.

This module provides a CLI for running essential gene analysis from the command line.

Usage:
    python -m auxotroph-design.cli analyze --config config.yaml
    python -m auxotroph-design.cli gapfill --config config.yaml
"""

import argparse
import logging
import sys
from pathlib import Path
from typing import Optional

from .essential_genes import EssentialGeneAnalyzer, create_bar_chart
from .config import AnalysisConfig, create_default_config


def setup_logging(verbose: bool = False) -> None:
    """Configure logging based on verbosity level."""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(sys.stdout),
        ]
    )


def run_analysis(config: AnalysisConfig) -> None:
    """
    Run the complete essential gene analysis pipeline.
    
    Args:
        config: Analysis configuration object
    """
    logger = logging.getLogger(__name__)
    
    # Initialize analyzer
    logger.info("Initializing analyzer...")
    analyzer = EssentialGeneAnalyzer(
        model_path=config.model_path,
        solver=config.solver,
        processes=config.processes
    )
    
    # Load carbon substrates
    if config.carbon_substrates_file:
        analyzer.load_carbon_substrates(
            config.carbon_substrates_file,
            name_column=config.substrate_name_column,
            reaction_column=config.substrate_reaction_column,
            filter_column=config.substrate_filter_column,
            filter_value=config.substrate_filter_value
        )
    
    # Load off reactions
    if config.off_reactions_file:
        analyzer.load_off_reactions(
            config.off_reactions_file,
            column=config.off_reactions_column
        )
    
    # Run analysis
    logger.info("Running essential gene analysis...")
    results = analyzer.analyze_all_substrates()
    
    # Export results
    output_dir = Path(config.output_dir)
    
    if config.export_format in ('csv', 'both'):
        csv_path = output_dir / 'essential_genes_matrix.csv'
        analyzer.export_results_to_csv(results, csv_path)
    
    if config.export_format in ('json', 'both'):
        json_path = output_dir / 'essential_genes.json'
        analyzer.export_results_to_json(results, json_path)
    
    # Find conditionally essential genes
    logger.info("Identifying conditionally essential genes...")
    cond_essential = analyzer.find_conditionally_essential_genes(
        reference_substrate=config.reference_substrate
    )
    
    # Export conditional essentiality
    cond_counts = {k: len(v) for k, v in cond_essential.items() if v}
    
    import json
    cond_path = output_dir / 'conditionally_essential.json'
    with open(cond_path, 'w') as f:
        json.dump({k: list(v) for k, v in cond_essential.items()}, f, indent=2)
    logger.info(f"Conditionally essential genes saved to {cond_path}")
    
    # Generate figures
    if config.generate_figures and cond_counts:
        fig_path = output_dir / 'conditionally_essential_genes.png'
        create_bar_chart(
            cond_counts,
            fig_path,
            title='Conditionally Essential Genes by Carbon Substrate'
        )
    
    # Generate comparison matrix
    matrix = analyzer.generate_comparison_matrix(results)
    matrix_path = output_dir / 'gene_essentiality_matrix.csv'
    matrix.to_csv(matrix_path)
    logger.info(f"Comparison matrix saved to {matrix_path}")
    
    logger.info("Analysis complete!")


def run_gapfilling(config: AnalysisConfig, gene_id: str, substrate: str) -> None:
    """
    Run gap-filling analysis for a specific gene knockout.
    
    Args:
        config: Analysis configuration object
        gene_id: Gene ID to knock out
        substrate: Carbon substrate for growth condition
    """
    logger = logging.getLogger(__name__)
    
    # Initialize analyzer
    analyzer = EssentialGeneAnalyzer(
        model_path=config.model_path,
        solver=config.solver,
        processes=config.processes
    )
    
    # Load carbon substrates
    if config.carbon_substrates_file:
        analyzer.load_carbon_substrates(
            config.carbon_substrates_file,
            name_column=config.substrate_name_column,
            reaction_column=config.substrate_reaction_column,
        )
    
    # Find substrate
    substrate_obj = None
    for sub in analyzer.carbon_substrates:
        if sub.name == substrate:
            substrate_obj = sub
            break
    
    if not substrate_obj:
        logger.error(f"Substrate not found: {substrate}")
        return
    
    # Get rescue reactions
    rescue_rxns = analyzer.identify_rescue_reactions(gene_id, substrate_obj)
    
    if rescue_rxns:
        logger.info(f"Rescue reactions for {gene_id} on {substrate}:")
        for rxn in rescue_rxns:
            logger.info(f"  - {rxn}")
    else:
        logger.info(f"No rescue reactions found for {gene_id}")


def main():
    """Main entry point for the CLI."""
    parser = argparse.ArgumentParser(
        description='Essential Gene Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
Examples:
  # Generate default configuration file
  python -m src.cli init-config --output config.yaml
  
  # Run full analysis
  python -m src.cli analyze --config config.yaml
  
  # Run gap-filling for a specific gene
  python -m src.cli gapfill --config config.yaml --gene b1234 --substrate "D-Glucose"
        '''
    )
    
    parser.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Enable verbose output'
    )
    
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # init-config command
    init_parser = subparsers.add_parser(
        'init-config',
        help='Create a default configuration file'
    )
    init_parser.add_argument(
        '--output', '-o',
        default='config.yaml',
        help='Output path for configuration file'
    )
    
    # analyze command
    analyze_parser = subparsers.add_parser(
        'analyze',
        help='Run essential gene analysis'
    )
    analyze_parser.add_argument(
        '--config', '-c',
        required=True,
        help='Path to configuration file'
    )
    
    # gapfill command
    gapfill_parser = subparsers.add_parser(
        'gapfill',
        help='Run gap-filling analysis'
    )
    gapfill_parser.add_argument(
        '--config', '-c',
        required=True,
        help='Path to configuration file'
    )
    gapfill_parser.add_argument(
        '--gene', '-g',
        required=True,
        help='Gene ID to analyze'
    )
    gapfill_parser.add_argument(
        '--substrate', '-s',
        required=True,
        help='Carbon substrate for growth condition'
    )
    
    args = parser.parse_args()
    
    setup_logging(args.verbose)
    
    if args.command == 'init-config':
        create_default_config(args.output)
        print(f"Configuration file created: {args.output}")
    
    elif args.command == 'analyze':
        if args.config.endswith('.yaml') or args.config.endswith('.yml'):
            config = AnalysisConfig.from_yaml(args.config)
        else:
            config = AnalysisConfig.from_json(args.config)
        run_analysis(config)
    
    elif args.command == 'gapfill':
        if args.config.endswith('.yaml') or args.config.endswith('.yml'):
            config = AnalysisConfig.from_yaml(args.config)
        else:
            config = AnalysisConfig.from_json(args.config)
        run_gapfilling(config, args.gene, args.substrate)
    
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
