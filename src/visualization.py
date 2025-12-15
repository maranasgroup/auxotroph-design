"""
Visualization module for Essential Gene Analysis.

This module provides publication-quality figure generation for gene essentiality
analysis results.
"""

from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd


def create_horizontal_bar_chart(
    data: Dict[str, int],
    output_path: Union[str, Path],
    title: str = "Conditionally Essential Genes",
    xlabel: str = "Number of Genes",
    figsize: Tuple[float, float] = (8, 10),
    color: str = '#2C5F2D',
    dpi: int = 300,
    font_size: int = 10
) -> None:
    """
    Create a publication-quality horizontal bar chart.
    
    Follows PLOS Computational Biology style guidelines.
    
    Args:
        data: Dictionary mapping category names to values
        output_path: Path for the output figure
        title: Chart title
        xlabel: X-axis label
        figsize: Figure size (width, height) in inches
        color: Bar color
        dpi: Resolution for output file
        font_size: Base font size
    """
    import matplotlib.pyplot as plt
    
    # Sort by values in descending order
    sorted_data = dict(sorted(data.items(), key=lambda x: x[1], reverse=True))
    
    # Filter out zeros
    sorted_data = {k: v for k, v in sorted_data.items() if v > 0}
    
    if not sorted_data:
        print("No data to plot")
        return
    
    # Set up the figure with publication-quality settings
    plt.rcParams.update({
        'font.size': font_size,
        'font.family': 'sans-serif',
        'axes.labelsize': font_size,
        'axes.titlesize': font_size + 2,
        'xtick.labelsize': font_size - 1,
        'ytick.labelsize': font_size - 1,
        'axes.linewidth': 0.8,
        'xtick.major.width': 0.8,
        'ytick.major.width': 0.8,
    })
    
    fig, ax = plt.subplots(figsize=figsize)
    
    y_pos = np.arange(len(sorted_data))
    bars = ax.barh(y_pos, list(sorted_data.values()), align='center', 
                   color=color, edgecolor='black', linewidth=0.5)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(list(sorted_data.keys()))
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontweight='bold')
    
    # Add value labels on bars
    for bar, value in zip(bars, sorted_data.values()):
        ax.text(bar.get_width() + 0.5, bar.get_y() + bar.get_height()/2,
                str(value), va='center', ha='left', fontsize=font_size - 2)
    
    # Clean up spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight', 
                facecolor='white', edgecolor='none')
    plt.close()


def create_vertical_bar_chart(
    data: Dict[str, int],
    output_path: Union[str, Path],
    title: str = "Conditionally Essential Genes",
    ylabel: str = "Number of Genes",
    figsize: Tuple[float, float] = (12, 6),
    color: str = '#2C5F2D',
    dpi: int = 300,
    font_size: int = 10,
    rotate_labels: bool = True
) -> None:
    """
    Create a publication-quality vertical bar chart.
    
    Args:
        data: Dictionary mapping category names to values
        output_path: Path for the output figure
        title: Chart title
        ylabel: Y-axis label
        figsize: Figure size (width, height) in inches
        color: Bar color
        dpi: Resolution for output file
        font_size: Base font size
        rotate_labels: Whether to rotate x-axis labels
    """
    import matplotlib.pyplot as plt
    
    # Sort by values in descending order
    sorted_data = dict(sorted(data.items(), key=lambda x: x[1], reverse=True))
    sorted_data = {k: v for k, v in sorted_data.items() if v > 0}
    
    if not sorted_data:
        print("No data to plot")
        return
    
    plt.rcParams.update({
        'font.size': font_size,
        'font.family': 'sans-serif',
    })
    
    fig, ax = plt.subplots(figsize=figsize)
    
    x_pos = np.arange(len(sorted_data))
    bars = ax.bar(x_pos, list(sorted_data.values()), align='center',
                  color=color, edgecolor='black', linewidth=0.5)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(list(sorted_data.keys()), 
                       rotation=45 if rotate_labels else 0,
                       ha='right' if rotate_labels else 'center')
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontweight='bold')
    
    # Add value labels on bars
    for bar, value in zip(bars, sorted_data.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                str(value), ha='center', va='bottom', fontsize=font_size - 2)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()


def create_heatmap(
    matrix: pd.DataFrame,
    output_path: Union[str, Path],
    title: str = "Gene Essentiality Matrix",
    figsize: Optional[Tuple[float, float]] = None,
    cmap: str = 'RdYlGn_r',
    dpi: int = 300,
    show_values: bool = False
) -> None:
    """
    Create a heatmap of gene essentiality across conditions.
    
    Args:
        matrix: DataFrame with genes as rows and conditions as columns
        output_path: Path for the output figure
        title: Chart title
        figsize: Figure size (auto-calculated if None)
        cmap: Colormap name
        dpi: Resolution for output file
        show_values: Whether to show values in cells
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    if figsize is None:
        width = max(8, len(matrix.columns) * 0.5)
        height = max(6, len(matrix) * 0.15)
        figsize = (width, height)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    sns.heatmap(
        matrix,
        cmap=cmap,
        ax=ax,
        cbar_kws={'label': 'Essential (1) / Non-essential (0)'},
        linewidths=0.5 if len(matrix) < 50 else 0,
        annot=show_values and len(matrix) < 30,
        fmt='d' if show_values else ''
    )
    
    ax.set_title(title, fontweight='bold')
    ax.set_xlabel('Carbon Substrate')
    ax.set_ylabel('Gene')
    
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()


def create_venn_diagram(
    sets: Dict[str, set],
    output_path: Union[str, Path],
    title: str = "Gene Set Overlap",
    figsize: Tuple[float, float] = (8, 8),
    dpi: int = 300
) -> None:
    """
    Create a Venn diagram showing overlap between gene sets.
    
    Note: Works best with 2-3 sets.
    
    Args:
        sets: Dictionary mapping set names to gene sets
        output_path: Path for the output figure
        title: Chart title
        figsize: Figure size (width, height)
        dpi: Resolution for output file
    """
    import matplotlib.pyplot as plt
    
    try:
        from matplotlib_venn import venn2, venn3
    except ImportError:
        print("matplotlib_venn not installed. Skipping Venn diagram.")
        return
    
    if len(sets) < 2 or len(sets) > 3:
        print("Venn diagram requires 2-3 sets")
        return
    
    fig, ax = plt.subplots(figsize=figsize)
    
    set_list = list(sets.values())
    set_labels = list(sets.keys())
    
    if len(sets) == 2:
        venn2(set_list, set_labels=set_labels, ax=ax)
    else:
        venn3(set_list, set_labels=set_labels, ax=ax)
    
    ax.set_title(title, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()


def create_upset_plot(
    sets: Dict[str, set],
    output_path: Union[str, Path],
    title: str = "Gene Set Intersections",
    figsize: Tuple[float, float] = (12, 8),
    dpi: int = 300
) -> None:
    """
    Create an UpSet plot for visualizing set intersections.
    
    Better than Venn diagrams for more than 3 sets.
    
    Args:
        sets: Dictionary mapping set names to gene sets
        output_path: Path for the output figure
        title: Chart title
        figsize: Figure size (width, height)
        dpi: Resolution for output file
    """
    import matplotlib.pyplot as plt
    
    try:
        from upsetplot import UpSet, from_contents
    except ImportError:
        print("upsetplot not installed. Skipping UpSet plot.")
        return
    
    data = from_contents(sets)
    
    fig = plt.figure(figsize=figsize)
    upset = UpSet(data, show_counts=True, sort_by='cardinality')
    upset.plot(fig=fig)
    
    plt.suptitle(title, fontweight='bold', y=1.02)
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()


def create_rescue_network(
    rescue_data: Dict[str, Dict[str, List[str]]],
    output_path: Union[str, Path],
    title: str = "Gene-Reaction Rescue Network",
    figsize: Tuple[float, float] = (14, 10),
    dpi: int = 300
) -> None:
    """
    Create a network visualization of genes and their rescue reactions.
    
    Args:
        rescue_data: Nested dict: {substrate: {gene_id: [rescue_reactions]}}
        output_path: Path for the output figure
        title: Chart title
        figsize: Figure size (width, height)
        dpi: Resolution for output file
    """
    import matplotlib.pyplot as plt
    
    try:
        import networkx as nx
    except ImportError:
        print("networkx not installed. Skipping network plot.")
        return
    
    G = nx.Graph()
    
    # Add nodes and edges
    for substrate, genes in rescue_data.items():
        for gene, reactions in genes.items():
            G.add_node(gene, node_type='gene', substrate=substrate)
            for rxn in reactions:
                G.add_node(rxn, node_type='reaction')
                G.add_edge(gene, rxn)
    
    if not G.nodes():
        print("No data for network plot")
        return
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Set up node colors
    node_colors = []
    for node in G.nodes():
        if G.nodes[node].get('node_type') == 'gene':
            node_colors.append('#E74C3C')  # Red for genes
        else:
            node_colors.append('#3498DB')  # Blue for reactions
    
    pos = nx.spring_layout(G, k=2, iterations=50)
    
    nx.draw(
        G, pos, ax=ax,
        node_color=node_colors,
        node_size=500,
        with_labels=True,
        font_size=8,
        font_weight='bold',
        edge_color='gray',
        alpha=0.7
    )
    
    ax.set_title(title, fontweight='bold')
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#E74C3C', label='Gene'),
        Patch(facecolor='#3498DB', label='Reaction')
    ]
    ax.legend(handles=legend_elements, loc='upper left')
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=dpi, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()


def generate_summary_statistics(
    results: Dict[str, List[str]],
    output_path: Union[str, Path]
) -> pd.DataFrame:
    """
    Generate summary statistics for essential gene analysis.
    
    Args:
        results: Dictionary mapping substrates to lists of essential genes
        output_path: Path for the output CSV file
        
    Returns:
        DataFrame with summary statistics
    """
    stats = []
    
    all_genes = set()
    for genes in results.values():
        all_genes.update(genes)
    
    for substrate, genes in results.items():
        gene_set = set(genes)
        stats.append({
            'Substrate': substrate,
            'Essential Genes': len(genes),
            'Unique Genes': len(gene_set - set().union(*[
                set(g) for s, g in results.items() if s != substrate
            ])),
            'Shared with All': len(gene_set.intersection(*[
                set(g) for g in results.values()
            ]))
        })
    
    df = pd.DataFrame(stats)
    df.to_csv(output_path, index=False)
    
    return df
