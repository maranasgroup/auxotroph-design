"""
Iterate through all carbon substrates and save OFF genes and OFF reactions.
"""
import pandas as pd
import json
from datetime import datetime
from boolean_gem_simple import BooleanGEMNetwork


def analyze_all_substrates(gene_rules_file, tf_activity_file, env_template_file,
                           derived_conditions_file, gem_path,
                           output_dir='substrate_analysis',
                           oxygen_conditions=[True, False]):
    """
    Run simulations for all carbon substrates and save OFF genes/reactions.
    
    Parameters:
    -----------
    gene_rules_file : str
    tf_activity_file : str
    env_template_file : str
    derived_conditions_file : str
    gem_path : str
    output_dir : str
        Directory to save results
    oxygen_conditions : list
        List of oxygen conditions to test [True] or [True, False]
    """
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    # All available carbon substrates
    carbon_substrates = [
        'glucose', 'arabinose', 'galactose', 'lactose', 'maltose',
        'xylose', 'glycerol', 'ribose', 'fucose', 'rhamnose',
        'trehalose', 'acetate', 'succinate', 'fumarate', 'malate',
        'lactate', 'formate', 'mannose', 'glycolate',
        'galacturonate', 'glucuronate', 'N-acetylglucosamine', 'deoxyribose'
    ]
    
    # Create network
    print("="*80)
    print("INITIALIZING NETWORK")
    print("="*80)
    network = BooleanGEMNetwork(
        gene_rules_file, tf_activity_file,
        env_template_file, derived_conditions_file,
        gem_path
    )
    
    # Store results
    all_results = []
    condition_summaries = []
    
    # Iterate through all conditions
    for oxygen in oxygen_conditions:
        oxygen_str = 'aerobic' if oxygen else 'anaerobic'
        
        for carbon_source in carbon_substrates:
            condition_name = f"{oxygen_str}_{carbon_source}"
            print("\n" + "="*80)
            print(f"ANALYZING: {condition_name}")
            print("="*80)
            
            try:
                # Run simulation
                gene_states, solution, iters = network.run_coupled_simulation(
                    carbon_source=carbon_source,
                    oxygen=oxygen,
                    identify_essential=True,  # Identify essential genes
                    force_essential_on=False,
                    debug_mode=False
                )
                
                # Get OFF genes
                off_genes = get_off_genes(network)
                
                # Get OFF reactions (blocked in GEM)
                off_reactions = get_off_reactions(network)
                
                # Get ON genes and reactions for reference
                on_genes = get_on_genes(network)
                on_reactions = get_on_reactions(network, solution)
                
                # Summary statistics
                summary = {
                    'condition': condition_name,
                    'carbon_source': carbon_source,
                    'oxygen': oxygen_str,
                    'converged': iters < 50,
                    'iterations': iters,
                    'growth_rate': solution.objective_value if solution.status == 'optimal' else 0,
                    'fba_status': solution.status,
                    'total_genes': len(network.gene_state),
                    'on_genes': len(on_genes),
                    'off_genes': len(off_genes),
                    'essential_genes': len(network.essential_genes),
                    'total_reactions': len(network.gem.reactions),
                    'on_reactions': len(on_reactions),
                    'off_reactions': len(off_reactions),
                    'active_tfs': len(network.get_active_tfs())
                }
                
                condition_summaries.append(summary)
                
                # Detailed results
                result = {
                    'condition': condition_name,
                    'carbon_source': carbon_source,
                    'oxygen': oxygen_str,
                    'growth_rate': summary['growth_rate'],
                    'off_genes': off_genes,
                    'off_reactions': off_reactions,
                    'on_genes': on_genes[:100],  # Limit to first 100 for file size
                    'on_reactions': on_reactions[:100],
                    'essential_genes': list(network.essential_genes),
                    'active_tfs': network.get_active_tfs()
                }
                
                all_results.append(result)
                
                # Save individual condition files
                save_condition_results(
                    condition_name, off_genes, off_reactions,
                    on_genes, on_reactions, network, solution,
                    output_dir
                )
                
                print(f"\n✓ {condition_name} completed:")
                print(f"  Growth: {summary['growth_rate']:.4f}")
                print(f"  OFF genes: {len(off_genes)}")
                print(f"  OFF reactions: {len(off_reactions)}")
                
            except Exception as e:
                print(f"\n✗ {condition_name} FAILED: {e}")
                condition_summaries.append({
                    'condition': condition_name,
                    'carbon_source': carbon_source,
                    'oxygen': oxygen_str,
                    'error': str(e),
                    'growth_rate': 0
                })
    
    # Save comprehensive summary
    save_comprehensive_summary(condition_summaries, all_results, output_dir)
    
    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print(f"Results saved to: {output_dir}/")
    print(f"  - summary.csv: Overview of all conditions")
    print(f"  - all_results.json: Detailed results")
    print(f"  - Individual files: {condition_name}_*.csv")


def get_off_genes(network):
    """Get list of OFF genes with details."""
    off_genes = []
    for _, row in network.gene_rules.iterrows():
        bnum = row['bNum']
        gene_name = row['Gene']
        if network.gene_state.get(bnum, 0) == 0:
            off_genes.append({
                'bNum': bnum,
                'gene_name': gene_name if pd.notna(gene_name) else '',
                'rule': row['Rule'],
                'is_essential': bnum in network.essential_genes
            })
    return off_genes


def get_on_genes(network):
    """Get list of ON genes with details."""
    on_genes = []
    for _, row in network.gene_rules.iterrows():
        bnum = row['bNum']
        gene_name = row['Gene']
        if network.gene_state.get(bnum, 0) == 1:
            on_genes.append({
                'bNum': bnum,
                'gene_name': gene_name if pd.notna(gene_name) else '',
                'rule': row['Rule'],
                'is_essential': bnum in network.essential_genes
            })
    return on_genes


def get_off_reactions(network):
    """Get list of OFF (blocked) reactions in GEM."""
    off_reactions = []
    for rxn in network.gem.reactions:
        # Check if reaction is blocked (both bounds are 0)
        if rxn.lower_bound == 0 and rxn.upper_bound == 0:
            off_reactions.append({
                'reaction_id': rxn.id,
                'reaction_name': rxn.name,
                'subsystem': rxn.subsystem,
                'gene_reaction_rule': rxn.gene_reaction_rule,
                'reversible': rxn.reversibility
            })
    return off_reactions


def get_on_reactions(network, solution):
    """Get list of ON (active) reactions with flux values."""
    on_reactions = []
    if solution.status == 'optimal':
        for rxn in network.gem.reactions:
            flux = solution.fluxes.get(rxn.id, 0)
            if abs(flux) > 1e-6:  # Has significant flux
                on_reactions.append({
                    'reaction_id': rxn.id,
                    'reaction_name': rxn.name,
                    'flux': flux,
                    'subsystem': rxn.subsystem,
                    'gene_reaction_rule': rxn.gene_reaction_rule
                })
    return on_reactions


def save_condition_results(condition_name, off_genes, off_reactions,
                          on_genes, on_reactions, network, solution,
                          output_dir):
    """Save results for a single condition to CSV files."""
    
    # Save OFF genes
    if off_genes:
        df_off_genes = pd.DataFrame(off_genes)
        df_off_genes.to_csv(f'{output_dir}/{condition_name}_OFF_genes.csv', index=False)
    
    # Save OFF reactions
    if off_reactions:
        df_off_reactions = pd.DataFrame(off_reactions)
        df_off_reactions.to_csv(f'{output_dir}/{condition_name}_OFF_reactions.csv', index=False)
    
    # Save ON genes
    if on_genes:
        df_on_genes = pd.DataFrame(on_genes)
        df_on_genes.to_csv(f'{output_dir}/{condition_name}_ON_genes.csv', index=False)
    
    # Save ON reactions
    if on_reactions:
        df_on_reactions = pd.DataFrame(on_reactions)
        df_on_reactions.to_csv(f'{output_dir}/{condition_name}_ON_reactions.csv', index=False)
    
    # Save active TFs
    active_tfs = network.get_active_tfs()
    if active_tfs:
        df_tfs = pd.DataFrame({'tf_name': active_tfs})
        df_tfs.to_csv(f'{output_dir}/{condition_name}_active_TFs.csv', index=False)


def save_comprehensive_summary(condition_summaries, all_results, output_dir):
    """Save comprehensive summary files."""
    
    # Summary CSV
    df_summary = pd.DataFrame(condition_summaries)
    df_summary.to_csv(f'{output_dir}/summary.csv', index=False)
    
    # Detailed JSON
    with open(f'{output_dir}/all_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    
    # Create comparison matrix: OFF genes across conditions
    create_off_genes_matrix(all_results, output_dir)
    
    # Create comparison matrix: OFF reactions across conditions
    create_off_reactions_matrix(all_results, output_dir)


def create_off_genes_matrix(all_results, output_dir):
    """Create a matrix showing which genes are OFF in which conditions."""
    
    # Collect all unique genes
    all_genes = set()
    for result in all_results:
        if 'off_genes' in result:
            for gene in result['off_genes']:
                all_genes.add(gene['bNum'])
    
    if not all_genes:
        return
    
    # Create matrix
    matrix_data = []
    for gene in sorted(all_genes):
        row = {'gene': gene}
        for result in all_results:
            condition = result['condition']
            off_gene_ids = [g['bNum'] for g in result.get('off_genes', [])]
            row[condition] = 1 if gene in off_gene_ids else 0
        matrix_data.append(row)
    
    df_matrix = pd.DataFrame(matrix_data)
    df_matrix.to_csv(f'{output_dir}/OFF_genes_matrix.csv', index=False)
    
    print(f"\n✓ Created OFF genes matrix: {len(all_genes)} genes × {len(all_results)} conditions")


def create_off_reactions_matrix(all_results, output_dir):
    """Create a matrix showing which reactions are OFF in which conditions."""
    
    # Collect all unique reactions
    all_reactions = set()
    for result in all_results:
        if 'off_reactions' in result:
            for rxn in result['off_reactions']:
                all_reactions.add(rxn['reaction_id'])
    
    if not all_reactions:
        return
    
    # Create matrix
    matrix_data = []
    for rxn_id in sorted(all_reactions):
        row = {'reaction': rxn_id}
        for result in all_results:
            condition = result['condition']
            off_rxn_ids = [r['reaction_id'] for r in result.get('off_reactions', [])]
            row[condition] = 1 if rxn_id in off_rxn_ids else 0
        matrix_data.append(row)
    
    df_matrix = pd.DataFrame(matrix_data)
    df_matrix.to_csv(f'{output_dir}/OFF_reactions_matrix.csv', index=False)
    
    print(f"✓ Created OFF reactions matrix: {len(all_reactions)} reactions × {len(all_results)} conditions")


def quick_analysis(gene_rules_file, tf_activity_file, env_template_file,
                  derived_conditions_file, gem_path,
                  substrates=['glucose', 'arabinose', 'galacturonate'],
                  output_dir='quick_analysis'):
    """
    Quick analysis of a few substrates (for testing).
    """
    print("Running quick analysis on selected substrates...")
    
    import os
    os.makedirs(output_dir, exist_ok=True)
    
    network = BooleanGEMNetwork(
        gene_rules_file, tf_activity_file,
        env_template_file, derived_conditions_file,
        gem_path
    )
    
    for carbon_source in substrates:
        condition_name = f"aerobic_{carbon_source}"
        print(f"\n{'='*60}")
        print(f"Analyzing: {condition_name}")
        print('='*60)
        
        gene_states, solution, iters = network.run_coupled_simulation(
            carbon_source=carbon_source,
            oxygen=True,
            identify_essential=True,
            debug_mode=False
        )
        
        # Get and save results
        off_genes = get_off_genes(network)
        off_reactions = get_off_reactions(network)
        on_genes = get_on_genes(network)
        on_reactions = get_on_reactions(network, solution)
        
        save_condition_results(
            condition_name, off_genes, off_reactions,
            on_genes, on_reactions, network, solution,
            output_dir
        )
        
        print(f"✓ Results saved:")
        print(f"  OFF genes: {len(off_genes)}")
        print(f"  OFF reactions: {len(off_reactions)}")
        print(f"  ON genes: {len(on_genes)}")
        print(f"  ON reactions: {len(on_reactions)}")


if __name__ == "__main__":
    import sys
    
    # Configuration
    gene_rules_file = 'gene_rules.csv'
    tf_activity_file = 'tf_activity.csv'
    env_template_file = 'environment_template.csv'
    derived_conditions_file = 'derived_conditions.csv'
    gem_path = 'iAF1260b.xml'
    
    # Choose mode
    if len(sys.argv) > 1 and sys.argv[1] == 'quick':
        # Quick mode - test with a few substrates
        quick_analysis(
            gene_rules_file, tf_activity_file,
            env_template_file, derived_conditions_file,
            gem_path,
            substrates=['glucose', 'arabinose', 'galacturonate'],
            output_dir='quick_analysis'
        )
    else:
        # Full mode - all substrates, aerobic and anaerobic
        analyze_all_substrates(
            gene_rules_file, tf_activity_file,
            env_template_file, derived_conditions_file,
            gem_path,
            output_dir='substrate_analysis',
            oxygen_conditions=[True, False]  # Both aerobic and anaerobic
        )
    
    print("\n" + "="*80)
    print("DONE!")
    print("="*80)
    print("\nCheck the output directory for results:")
    print("  - summary.csv: Quick overview")
    print("  - *_OFF_genes.csv: OFF genes per condition")
    print("  - *_OFF_reactions.csv: OFF reactions per condition")
    print("  - *_ON_genes.csv: ON genes per condition")
    print("  - *_ON_reactions.csv: ON reactions per condition")
    print("  - OFF_genes_matrix.csv: Gene OFF/ON across all conditions")
    print("  - OFF_reactions_matrix.csv: Reaction OFF/ON across all conditions")
