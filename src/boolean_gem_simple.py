import pandas as pd
import re
import cobra
from cobra.exceptions import OptimizationError


class BooleanGEMNetwork:
    def __init__(self, gene_rules_file, tf_activity_file, env_template_file,
                 derived_conditions_file, gem_path):
        """Initialize the integrated Boolean-GEM network."""
        self.gene_rules = pd.read_csv(gene_rules_file, sep='\t')
        self.tf_activity = pd.read_csv(tf_activity_file, encoding='utf-8-sig')
        self.env_template = pd.read_csv(env_template_file)
        self.derived_conditions = pd.read_csv(derived_conditions_file, sep='\t')

        # Load GEM
        print(f"Loading GEM from {gem_path}...")
        self.gem = cobra.io.read_sbml_model(gem_path)
        print(f"  Loaded model: {self.gem.id}")
        print(f"  Reactions: {len(self.gem.reactions)}")
        print(f"  Metabolites: {len(self.gem.metabolites)}")
        print(f"  Genes: {len(self.gem.genes)}")

        self.gene_state = {}
        self.tf_state = {}
        self.env_state = {}

        # Essential genes (never knocked out)
        self.essential_genes = set()

        # Metabolite mapping: Boolean ID -> GEM metabolite ID
        self.metabolite_map = {
            'AGDC': 'acgam6p_c',  # N-Acetyl-D-glucosamine-6-P
            'ALTRH': 'altr_c',  # Altronate
            'MNNH': 'manna_c',  # Mannonate
            'GUI1': 'guln__D_c',  # D-Gulonate
            'GUI2': 'guln__L_c',  # L-Gulonate
            'MANAO': 'manna_c',  # Mannonate (duplicate of MNNH)
            'PPM2': 'drib_c',  # Deoxyribose
            'TAGURr': 'tagur_c'  # Tagaturonate
        }

        # Carbon source mapping
        self.carbon_exchange_map = {
            'glucose': 'EX_glc__D_e',
            'arabinose': 'EX_arab__L_e',
            'galactose': 'EX_gal_e',
            'lactose': 'EX_lcts_e',
            'maltose': 'EX_malt_e',
            'xylose': 'EX_xyl__D_e',
            'glycerol': 'EX_glyc_e',
            'ribose': 'EX_rib__D_e',
            'fucose': 'EX_fuc__L_e',
            'rhamnose': 'EX_rmn_e',
            'trehalose': 'EX_tre_e',
            'acetate': 'EX_ac_e',
            'succinate': 'EX_succ_e',
            'fumarate': 'EX_fum_e',
            'malate': 'EX_mal__L_e',
            'lactate': 'EX_lac__D_e',
            'formate': 'EX_for_e',
            'mannose': 'EX_man_e',
            'glycolate': 'EX_glyclt_e',
            'galacturonate': 'EX_galur_e',
            'glucuronate': 'EX_glcur_e',
            'N-acetylglucosamine': 'EX_acgam_e',
            'deoxyribose': 'EX_drib_e',
        }

    def identify_essential_genes(self, carbon_source='glucose', oxygen=True,
                                 growth_threshold=0.001, max_test_genes=None):
        """
        Identify essential genes by single knockout analysis.

        Parameters:
        -----------
        carbon_source : str
            Carbon source for testing
        oxygen : bool
            Aerobic or anaerobic
        growth_threshold : float
            Minimum growth to consider viable (default: 0.001)
        max_test_genes : int or None
            Maximum genes to test (for faster testing). None = test all

        Returns:
        --------
        set : Essential gene IDs
        """
        print(f"\nIdentifying essential genes on {carbon_source} "
              f"({'aerobic' if oxygen else 'anaerobic'})...")

        # Setup medium
        self._setup_gem_medium(carbon_source, oxygen)

        # Get wildtype growth
        try:
            wt_solution = self.gem.optimize()
            wt_growth = wt_solution.objective_value
            print(f"  Wildtype growth: {wt_growth:.4f}")

            if wt_growth < growth_threshold:
                print(f"  WARNING: Wildtype growth is below threshold!")

        except OptimizationError:
            print("  ERROR: Wildtype not viable!")
            return set()

        # essential = set()
        total_genes = len(self.gem.genes)
        essential_genes = cobra.flux_analysis.variability.find_essential_genes(self.gem)
        essential = [gene.id for gene in essential_genes]
        # Limit number of genes to test if requested
        genes_to_test = list(self.gem.genes)
        # if max_test_genes and max_test_genes < total_genes:
        #     genes_to_test = genes_to_test[:max_test_genes]
        #     print(f"  Testing {max_test_genes} of {total_genes} genes (fast mode)...")
        # else:
        #     print(f"  Testing all {total_genes} genes (this may take a while)...")
        #
        # # Test each gene
        # for i, gene in enumerate(genes_to_test):
        #     if (i + 1) % 50 == 0:
        #         print(f"    Progress: {i + 1}/{len(genes_to_test)}, "
        #               f"{len(essential)} essential")
        #
        #     # Test knockout
        #     with self.gem as model:
        #         gene.knock_out()
        #         try:
        #             solution = model.optimize()
        #             growth = solution.objective_value
        #
        #             if growth < growth_threshold:
        #                 essential.add(gene.id)
        #
        #         except OptimizationError:
        #             essential.add(gene.id)

        self.essential_genes = essential

        pct = 100 * len(essential) / len(genes_to_test)
        print(f"  ✓ Found {len(essential)} essential genes ({pct:.1f}%)")

        if 0 < len(essential) <= 20:
            print(f"  Essential genes: {', '.join(sorted(essential))}")
        elif len(essential) > 20:
            examples = sorted(essential)[:10]
            print(f"  Examples: {', '.join(examples)}...")

        return essential

    def initialize_conditions(self, carbon_source='glucose', oxygen=True,
                              identify_essential=True, force_essential_on=False):
        """
        Initialize conditions.

        Parameters:
        -----------
        carbon_source : str
        oxygen : bool
        identify_essential : bool
            If True, identify essential genes (recommended, but slow)
        force_essential_on : bool
            If True, force essential genes ON in Boolean network
            If False, just protect them from knockout in GEM
        """
        carbon_map = {
            'glucose': 'glc-D(e) > 0',
            'arabinose': 'arab-L(e) > 0',
            'galactose': 'gal(e) > 0',
            'lactose': 'lcts(e) > 0',
            'maltose': 'malt(e) > 0',
            'xylose': 'xyl-D(e) > 0',
            'glycerol': 'glyc(e) > 0',
            'ribose': 'rib-D(e) > 0',
            'fucose': 'fuc-L(e) > 0',
            'rhamnose': 'rmn(e) > 0',
            'trehalose': 'tre(e) > 0',
            'acetate': 'ac(e) > 0',
            'succinate': 'succ(e) > 0',
            'fumarate': 'fum(e) > 0',
            'malate': 'mal-L(e) > 0',
            'lactate': 'lac-D(e) > 0',
            'formate': 'for(e) > 0',
            'mannose': 'man(e) > 0',
            'glycolate': 'glyclt(e) > 0',
            'galacturonate': 'galur(e) > 0',
            'glucuronate': 'glcur(e) > 0',
            'N-acetylglucosamine': 'acgam(e) > 0',
            'deoxyribose': 'cytd(e) > 0',
        }

        print(f"\nInitializing {'aerobic' if oxygen else 'anaerobic'}_{carbon_source}...")

        # Initialize genes to OFF
        for _, row in self.gene_rules.iterrows():
            self.gene_state[row['bNum']] = 0
            if pd.notna(row['Gene']) and row['Gene']:
                self.gene_state[row['Gene']] = 0

        # Initialize environment
        for _, row in self.env_template.iterrows():
            self.env_state[row['condition']] = False

        # Set oxygen
        self.env_state['o2(e) > 0'] = oxygen

        # Set carbon source
        if carbon_source in carbon_map:
            self.env_state[carbon_map[carbon_source]] = True
        else:
            raise ValueError(f"Unknown carbon source: {carbon_source}")

        # Growth conditions
        self.env_state['Growth > 0'] = True
        self.env_state['Growth>0'] = True
        self.env_state['NOT (Growth > 0)'] = False
        self.env_state['NOT (Growth>0)'] = False
        self.env_state['pH < 4'] = False
        self.env_state['pH < 7'] = False
        self.env_state['pi(e)<0.004E-6 M'] = False

        # Special conditions
        for cond in ['Stringent', 'Rich Medium', 'Stress', 'Oxidative Stress',
                     'Salicylate', 'dipyridyl', 'high NAD', 'high osmolarity',
                     'Surplus PYR', 'Surplus FDP', 'heat_shock']:
            self.env_state[cond] = False
            self.env_state[f'"{cond}"'] = False

        # Initialize internal metabolites
        for metab in ['AGDC', 'ALTRH', 'GUI1', 'GUI2', 'MANAO', 'MNNH', 'PPM2', 'TAGURr']:
            for op in ['>0', '> 0', '<0', '< 0']:
                self.env_state[f'{metab}{op}'] = False

        self._update_derived_conditions()

        # Initialize TFs
        for _, row in self.tf_activity.iterrows():
            self.tf_state[row['tf_name']] = 0
        self.tf_state['ON'] = 1
        self.tf_state['OFF'] = 0

        # Setup GEM medium
        self._setup_gem_medium(carbon_source, oxygen)

        # Identify essential genes (slow!)
        if identify_essential and len(self.essential_genes) == 0:
            self.identify_essential_genes(carbon_source, oxygen)

        # Optionally force essential genes ON
        if force_essential_on:
            n_forced = 0
            for gene_id in self.essential_genes:
                if gene_id in self.gene_state:
                    self.gene_state[gene_id] = 1
                    n_forced += 1
            print(f"  Forced {n_forced} essential genes to ON in Boolean network")

    def _setup_gem_medium(self, carbon_source, oxygen):
        """Configure GEM medium."""
        # self.gem.medium = {}

        # Carbon source
        if carbon_source in self.carbon_exchange_map:
            ex_rxn = self.carbon_exchange_map[carbon_source]
            if ex_rxn in self.gem.reactions:
                self.gem.medium[ex_rxn] = 10.0

        # Oxygen
        if oxygen and 'EX_o2_e' in self.gem.reactions:
            self.gem.medium['EX_o2_e'] = 20.0

        # Essential nutrients
        essential = {
            'EX_pi_e': 10.0, 'EX_nh4_e': 10.0, 'EX_so4_e': 10.0,
            'EX_k_e': 10.0, 'EX_na1_e': 10.0, 'EX_mg2_e': 10.0,
            'EX_ca2_e': 10.0, 'EX_cl_e': 10.0, 'EX_fe2_e': 10.0,
            'EX_h2o_e': 1000.0, 'EX_h_e': 1000.0,
        }

        for rxn_id, uptake in essential.items():
            if rxn_id in self.gem.reactions:
                self.gem.medium[rxn_id] = uptake

    def _apply_gene_states_to_gem(self):
        """
        Apply Boolean gene states to GEM.
        ESSENTIAL GENES ARE NEVER KNOCKED OUT.
        """
        knocked_out = 0
        protected = 0

        for gene_id, state in self.gene_state.items():
            if gene_id not in self.gem.genes:
                continue

            gene = self.gem.genes.get_by_id(gene_id)

            if state == 0:  # Gene OFF in Boolean network
                # Check if essential
                if gene_id in self.essential_genes:
                    protected += 1
                    # DON'T knockout - keep gene functional
                    continue
                else:
                    with self.gem as model:
                        model.genes.get_by_id(gene_id)
                        optim = model.optimize()
                        if optim.objective_value < 0.2:
                            self.gene_state[gene_id] = 1
                            continue
                        else:
                            # Knockout non-essential gene
                            gene.knock_out()

                            knocked_out += 1
            else:
                # Gene ON - make sure it's functional
                gene.functional = True

        print(f"    GEM: {knocked_out} genes knocked out, {protected} essential genes protected")

    def _extract_metabolite_states(self, solution, flux_threshold=0.001):
        """Extract metabolite states from FBA solution."""
        if solution.status != 'optimal':
            print(f"    Warning: FBA status = {solution.status}")
            return

        active_metabolites = []

        for bool_id, gem_id in self.metabolite_map.items():
            if gem_id not in self.gem.metabolites:
                continue

            metabolite = self.gem.metabolites.get_by_id(gem_id)

            # Sum fluxes through this metabolite
            total_flux = 0.0
            for rxn in metabolite.reactions:
                if rxn.id in solution.fluxes:
                    flux = abs(solution.fluxes[rxn.id])
                    coeff = rxn.metabolites.get(metabolite, 0)
                    if coeff != 0:
                        total_flux += flux * abs(coeff)

            # Set Boolean state
            is_active = total_flux > flux_threshold
            self.env_state[f'{bool_id}>0'] = is_active
            self.env_state[f'{bool_id} > 0'] = is_active

            if is_active:
                active_metabolites.append(f"{bool_id}={total_flux:.4f}")

        if active_metabolites:
            print(f"    Active metabolites: {', '.join(active_metabolites)}")

    def _update_derived_conditions(self):
        """Update derived CRP conditions."""
        for _, row in self.derived_conditions.iterrows():
            cond_name = row['condition_name']
            rule = row['rule']
            try:
                result = self._simple_evaluate(rule)
                self.env_state[cond_name] = result
                self.env_state[f'"{cond_name}"'] = result
            except:
                self.env_state[cond_name] = False
                self.env_state[f'"{cond_name}"'] = False

    def _simple_evaluate(self, rule, verbose_errors=False):
        """
        Simple rule evaluation with robust error handling.

        Parameters:
        -----------
        rule : str
            Boolean rule to evaluate
        verbose_errors : bool
            If True, print detailed error messages
        """
        if pd.isna(rule) or rule == '' or rule == 'ON':
            return True
        if rule == 'OFF' or rule == '(OFF)':
            return False

        rule_str = str(rule).strip()
        namespace = {'__builtins__': {}}

        # Add all gene states to namespace
        for key, val in self.gene_state.items():
            safe_key = str(key).replace('-', '_').replace(' ', '_')
            namespace[safe_key] = bool(val)
            namespace[str(key)] = bool(val)  # Keep original too

        # Add all TF states to namespace
        for key, val in self.tf_state.items():
            safe_key = str(key).replace('-', '_').replace(' ', '_')
            namespace[safe_key] = bool(val)
            namespace[str(key)] = bool(val)

        # Add all environment states to namespace
        for key, val in self.env_state.items():
            safe_key = str(key).replace('-', '_').replace(' ', '_')
            namespace[safe_key] = bool(val)
            namespace[str(key)] = bool(val)

        # Replace logical operators
        rule_str = rule_str.replace(' AND ', ' and ')
        rule_str = rule_str.replace(' OR ', ' or ')
        rule_str = rule_str.replace('NOT ', 'not ')
        rule_str = rule_str.replace('NOT(', 'not (')

        # Handle quoted strings
        def replace_quote(match):
            content = match.group(1)
            clean = content.replace(' ', '_').replace('-', '_')
            namespace[clean] = self.env_state.get(f'"{content}"', self.env_state.get(content, False))
            return clean

        rule_str = re.sub(r'"([^"]+)"', replace_quote, rule_str)

        # Handle metabolite comparisons like "glc-D(e) > 0"
        pattern = r'([a-zA-Z0-9_\-]+)\(e\)\s*([><])\s*(\d+(?:\.\d+)?(?:E-?\d+)?(?:\s*M)?)'

        def replace_comparison(match):
            var, op, val = match.groups()
            full_cond = f'{var}(e) {op} {val}'
            if full_cond in self.env_state:
                return str(self.env_state[full_cond])
            full_cond_no_space = f'{var}(e){op}{val}'
            if full_cond_no_space in self.env_state:
                return str(self.env_state[full_cond_no_space])
            return 'False'

        rule_str = re.sub(pattern, replace_comparison, rule_str)

        # Handle internal metabolite comparisons like "AGDC > 0"
        pattern2 = r'([A-Z0-9]+)\s*([><])\s*(\d+(?:\.\d+)?)'

        def replace_comparison2(match):
            var, op, val = match.groups()
            full_cond = f'{var} {op} {val}'
            if full_cond in self.env_state:
                return str(self.env_state[full_cond])
            full_cond_no_space = f'{var}{op}{val}'
            if full_cond_no_space in self.env_state:
                return str(self.env_state[full_cond_no_space])
            return 'False'

        rule_str = re.sub(pattern2, replace_comparison2, rule_str)

        # Replace any remaining gene/TF names that have hyphens
        # Find all potential identifiers with hyphens
        identifier_pattern = r'\b([a-zA-Z][a-zA-Z0-9\-]*)\b'

        def replace_identifier(match):
            identifier = match.group(1)
            # Check if it's a keyword
            if identifier.lower() in ['and', 'or', 'not', 'true', 'false']:
                return identifier
            # Check if it exists in our states
            if identifier in namespace:
                return identifier
            # Try with underscores instead of hyphens
            safe_identifier = identifier.replace('-', '_')
            if safe_identifier in namespace:
                return safe_identifier
            # If not found anywhere, assume it's False
            namespace[identifier] = False
            return identifier

        rule_str = re.sub(identifier_pattern, replace_identifier, rule_str)

        try:
            result = eval(rule_str, namespace)
            return bool(result)
        except Exception as e:
            if verbose_errors:
                print(f"    Error evaluating rule: {rule}")
                print(f"    Processed as: {rule_str}")
                print(f"    Error: {e}")
            return False

    def update_tf_activity(self, verbose_errors=True):
        """Update TF activities."""
        changed = False
        for _, row in self.tf_activity.iterrows():
            tf_name = row['tf_name']
            rule = row['rule']
            new_state = 1 if self._simple_evaluate(rule, verbose_errors) else 0
            if self.tf_state[tf_name] != new_state:
                self.tf_state[tf_name] = new_state
                changed = True
        return changed

    def update_gene_expression(self, verbose_errors=False):
        """Update gene expression."""
        changed = False
        for _, row in self.gene_rules.iterrows():
            bnum = row['bNum']
            rule = row['Rule']
            new_state = 1 if self._simple_evaluate(rule, verbose_errors) else 0
            if self.gene_state[bnum] != new_state:
                self.gene_state[bnum] = new_state
                if pd.notna(row['Gene']) and row['Gene']:
                    self.gene_state[row['Gene']] = new_state
                changed = True
        return changed

    def run_coupled_simulation(self, carbon_source='glucose', oxygen=True,
                               max_iterations=50, flux_threshold=0.001,
                               identify_essential=True, force_essential_on=False,
                               debug_mode=False):
        """
        Run coupled Boolean-GEM simulation.

        Parameters:
        -----------
        carbon_source : str
        oxygen : bool
        max_iterations : int
        flux_threshold : float
            Minimum flux to consider metabolite active
        identify_essential : bool
            If True, identify essential genes (SLOW but prevents crashes)
        force_essential_on : bool
            If True, force essential genes ON in Boolean network
        debug_mode : bool
            If True, print detailed error messages for rule evaluation failures
        """
        self.initialize_conditions(carbon_source, oxygen,
                                   identify_essential, force_essential_on)

        for iteration in range(max_iterations):
            print(f"\n--- Iteration {iteration + 1} ---")

            # 1. Update Boolean network
            tf_changed = self.update_tf_activity(verbose_errors=debug_mode)
            gene_changed = self.update_gene_expression(verbose_errors=debug_mode)

            # 2. Apply gene states to GEM (protecting essential genes)
            self._apply_gene_states_to_gem()

            # 3. Run FBA
            try:
                solution = self.gem.optimize()
                print(f"    FBA: growth={solution.objective_value:.4f}, status={solution.status}")
            except OptimizationError as e:
                print(f"    FBA FAILED: {e}")
                print(f"    Too many genes knocked out - model infeasible!")
                break

            # 4. Extract metabolite states
            self._extract_metabolite_states(solution, flux_threshold)

            # 5. Update derived conditions
            self._update_derived_conditions()

            expressed = sum(1 for v in self.gene_state.values() if v == 1)
            active_tfs = sum(1 for k, v in self.tf_state.items()
                             if v == 1 and k not in ['ON', 'OFF'])
            print(f"    Boolean: {expressed} genes, {active_tfs} TFs active")

            if not tf_changed and not gene_changed:
                print(f"\n✓ Converged after {iteration + 1} iterations")
                return self.gene_state, solution, iteration + 1

        print(f"\n⚠ Did not converge within {max_iterations} iterations")
        return self.gene_state, solution, max_iterations

    def get_expressed_genes(self):
        """Get list of expressed genes."""
        expressed = []
        for _, row in self.gene_rules.iterrows():
            if self.gene_state.get(row['bNum'], 0) == 1:
                expressed.append((row['bNum'], row['Gene']))
        return expressed

    def get_active_tfs(self):
        """Get list of active TFs."""
        return [tf for tf, state in self.tf_state.items()
                if state == 1 and tf not in ['ON', 'OFF']]


# Example usage
if __name__ == "__main__":
    network = BooleanGEMNetwork(
        gene_rules_file='gene_rules.csv',
        tf_activity_file='tf_activity.csv',
        env_template_file='environment_template.csv',
        derived_conditions_file='derived_conditions.csv',
        gem_path='iAF1260b.xml'
    )

    # Test with glucose
    print("\n" + "=" * 80)
    print("GLUCOSE (aerobic)")
    print("=" * 80)
    gene_states, solution, iters = network.run_coupled_simulation(
        'glucose', oxygen=True,
        identify_essential=True,  # Find essential genes (slow first time)
        force_essential_on=True,
        debug_mode=True  # Don't force them ON, just protect
    )
    print(f"\nExpressed genes: {len(network.get_expressed_genes())}")
    print(f"Active TFs: {', '.join(sorted(network.get_active_tfs()))}")
    print(f"Essential genes protected: {len(network.essential_genes)}")
