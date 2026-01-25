import ROOT
import numpy as np
from collections import defaultdict

def load_data(file_path, tree_name, variable_name):
    """Load a single variable from a ROOT file's TTree as a 2D array (events x objects)."""
    file = ROOT.TFile(file_path)
    tree = file.Get(tree_name)
    data = []
    for event in tree:
        value = getattr(event, variable_name)
        try:
            arr = list(value)
        except TypeError:
            arr = [value]
        data.append(arr)
    return np.array(data, dtype=object)

def apply_pt_eta_cuts(pts, etas, charges, pt_cut=15000):
    """Apply pT cut and return filtered arrays."""
    mask = pts > pt_cut
    return pts[mask], etas[mask], charges[mask]

def check_quality(quality_arr, n_required, indices=None):
    """Check if leptons pass quality requirements."""
    if len(quality_arr) < n_required:
        return False
    if indices is None:
        indices = range(n_required)
    return all(quality_arr[i] for i in indices if i < len(quality_arr))

def count_jets(jet_pts, threshold=30000):
    """Count jets above pT threshold."""
    return np.sum(np.array(jet_pts) > threshold)

class ChannelSelector:
    """Handles event selection and data collection for analysis channels."""
    
    def __init__(self):
        self.channels = {
            'dielectron': {'counts': [0, 0, 0, 0], 'data': defaultdict(list)},
            'electronmuon': {'counts': [0, 0, 0, 0], 'data': defaultdict(list)},
            'dimuon': {'counts': [0, 0, 0, 0], 'data': defaultdict(list)},
            'mu1jet': {'counts': [0, 0, 0, 0], 'data': defaultdict(list)},
            'el1jet': {'counts': [0, 0, 0, 0], 'data': defaultdict(list)},
        }
    def calculate_met(self, nu_pt, nu_phi, nu_eta):
        """Calculate generator-level MET from neutrinos."""
        nu_pt = np.array(nu_pt)
        nu_phi = np.array(nu_phi)
        nu_eta = np.array(nu_eta)
        nu_pt = nu_pt * np.sin(2*np.arctan(np.exp(-nu_eta)))
        nu_px = nu_pt * np.cos(nu_phi)
        nu_py = nu_pt * np.sin(nu_phi)
        met_x = np.sum(nu_px)
        met_y = np.sum(nu_py)
        met = np.sqrt(met_x**2 + met_y**2)
        met_phi = np.arctan2(met_y, met_x)
        return met, met_phi
    
    def add_lepton_data(self, channel, pts, etas, centrality_flag):
        """Add lepton kinematic data."""
        if isinstance(pts, (list, np.ndarray)):
            self.channels[channel]['data']['pts'].extend(pts if isinstance(pts, list) else pts.tolist())
            self.channels[channel]['data']['etas'].extend(etas if isinstance(etas, list) else etas.tolist())
        else:
            self.channels[channel]['data']['pts'].append(pts)
            self.channels[channel]['data']['etas'].append(etas)
        self.channels[channel]['data']['centrality'].append(centrality_flag)
    
    def add_event_data(self, channel, jet_mult, met_tots, met_phis, centrality_flag):
        """Add per-event data (jets, MET, centrality)."""
        self.channels[channel]['data']['jet_mults'].append(jet_mult)
        self.channels[channel]['data']['centrality'].append(centrality_flag)
        for i in range(1, 6):
            self.channels[channel]['data'][f'met{i}_pts'].append(met_tots[i-1])
            self.channels[channel]['data'][f'met{i}_etas'].append(met_phis[i-1])
    
    def process_dielectron(self, el_data, jet_pts, met_data, centrality_flag, nu_data):
        """Process dielectron channel."""
        el_pts, el_etas, el_charges = apply_pt_eta_cuts(
            np.array(el_data['pts']), np.array(el_data['etas']), np.array(el_data['charges'])
        )
        
        if len(el_pts) == 2 and np.prod(el_charges) < 0:
            jet_mult = count_jets(jet_pts)
            if jet_mult > 0:
                self.channels['dielectron']['counts'][0] += 1
                
                # Check quality requirements
                for op_idx, op_key in enumerate(['loose', 'medium', 'tight'], 1):
                    if check_quality(el_data[op_key], 2):
                        self.channels['dielectron']['counts'][op_idx] += 1
                
                # Store data only for loose events
                if check_quality(el_data['loose'], 2):
                    self.add_lepton_data('dielectron', el_pts, el_etas, centrality_flag)
                    self.add_event_data('dielectron', jet_mult, met_data['tots'], met_data['phis'], centrality_flag)
                    met_tot, met_phi = self.calculate_met(
                        nu_data['pt'], nu_data['phi'], nu_data['eta']
                    )
                    self.channels['dielectron']['data']['gen_met'].append(met_tot)
                    self.channels['dielectron']['data']['gen_met_phi'].append(met_phi)

    def process_electronmuon(self, el_data, mu_data, jet_pts, met_data, centrality_flag, nu_data):
        """Process electron-muon channel."""
        el_pts, el_etas, el_charges = apply_pt_eta_cuts(
            np.array(el_data['pts']), np.array(el_data['etas']), np.array(el_data['charges'])
        )
        mu_pts, mu_etas, mu_charges = apply_pt_eta_cuts(
            np.array(mu_data['pts']), np.array(mu_data['etas']), np.array(mu_data['charges'])
        )
        
        if len(el_pts) == 1 and len(mu_pts) == 1 and el_charges[0] * mu_charges[0] < 0:
            jet_mult = count_jets(jet_pts)
            if jet_mult > 0:
                self.channels['electronmuon']['counts'][0] += 1
                
                for op_idx, op_key in enumerate(['loose', 'medium', 'tight'], 1):
                    if check_quality(el_data[op_key], 1) and check_quality(mu_data[op_key], 1):
                        self.channels['electronmuon']['counts'][op_idx] += 1
                
                if check_quality(el_data['loose'], 1) and check_quality(mu_data['loose'], 1):
                    self.add_lepton_data('electronmuon', [el_pts[0], mu_pts[0]], [el_etas[0], mu_etas[0]], centrality_flag)
                    self.add_event_data('electronmuon', jet_mult, met_data['tots'], met_data['phis'], centrality_flag)
                    
                    # Calculate generator-level MET from neutrinos
                    met_tot, met_phi = self.calculate_met(
                        nu_data['pt'], nu_data['phi'], nu_data['eta']
                    )
                    self.channels['electronmuon']['data']['gen_met'].append(met_tot)
                    self.channels['electronmuon']['data']['gen_met_phi'].append(met_phi)
    
    def process_dimuon(self, mu_data, jet_pts, met_data, centrality_flag, nu_data):
        """Process dimuon channel."""
        mu_pts, mu_etas, mu_charges = apply_pt_eta_cuts(
            np.array(mu_data['pts']), np.array(mu_data['etas']), np.array(mu_data['charges'])
        )
        
        if len(mu_pts) == 2 and np.prod(mu_charges) < 0:
            jet_mult = count_jets(jet_pts)
            if jet_mult > 0:
                self.channels['dimuon']['counts'][0] += 1
                
                for op_idx, op_key in enumerate(['loose', 'medium', 'tight'], 1):
                    if check_quality(mu_data[op_key], 2):
                        self.channels['dimuon']['counts'][op_idx] += 1
                
                if check_quality(mu_data['loose'], 2):
                    self.add_lepton_data('dimuon', mu_pts, mu_etas, centrality_flag)
                    self.add_event_data('dimuon', jet_mult, met_data['tots'], met_data['phis'], centrality_flag)
                    
                    # Calculate generator-level MET from neutrinos
                    met_tot, met_phi = self.calculate_met(
                        nu_data['pt'], nu_data['phi'], nu_data['eta']
                    )
                    self.channels['dimuon']['data']['gen_met'].append(met_tot)
                    self.channels['dimuon']['data']['gen_met_phi'].append(met_phi)
    
    def process_mu1jet(self, mu_data, jet_pts, met_data, centrality_flag, nu_data):
        """Process single muon + jets channel."""
        mu_pts, mu_etas, _ = apply_pt_eta_cuts(
            np.array(mu_data['pts']), np.array(mu_data['etas']), np.array(mu_data['charges'])
        )
        
        if len(mu_pts) == 1:
            jet_mult = count_jets(jet_pts)
            if jet_mult >= 2:
                self.channels['mu1jet']['counts'][0] += 1
                
                for op_idx, op_key in enumerate(['loose', 'medium', 'tight'], 1):
                    if check_quality(mu_data[op_key], 1):
                        self.channels['mu1jet']['counts'][op_idx] += 1
                
                if check_quality(mu_data['loose'], 1):
                    self.add_lepton_data('mu1jet', mu_pts[0], mu_etas[0], centrality_flag)
                    self.add_event_data('mu1jet', jet_mult, met_data['tots'], met_data['phis'], centrality_flag)
                    
                    # Calculate generator-level MET from neutrinos
                    met_tot, met_phi = self.calculate_met(
                        nu_data['pt'], nu_data['phi'], nu_data['eta']
                    )
                    self.channels['mu1jet']['data']['gen_met'].append(met_tot)
                    self.channels['mu1jet']['data']['gen_met_phi'].append(met_phi)
    
    def process_el1jet(self, el_data, jet_pts, met_data, centrality_flag, nu_data):
        """Process single electron + jets channel."""
        el_pts, el_etas, _ = apply_pt_eta_cuts(
            np.array(el_data['pts']), np.array(el_data['etas']), np.array(el_data['charges'])
        )
        
        if len(el_pts) == 1:
            jet_mult = count_jets(jet_pts)
            if jet_mult >= 2:
                self.channels['el1jet']['counts'][0] += 1
                
                for op_idx, op_key in enumerate(['loose', 'medium', 'tight'], 1):
                    if check_quality(el_data[op_key], 1):
                        self.channels['el1jet']['counts'][op_idx] += 1
                
                if check_quality(el_data['loose'], 1):
                    self.add_lepton_data('el1jet', el_pts[0], el_etas[0], centrality_flag)
                    self.add_event_data('el1jet', jet_mult, met_data['tots'], met_data['phis'], centrality_flag)
                    
                    # Calculate generator-level MET from neutrinos
                    met_tot, met_phi = self.calculate_met(
                        nu_data['pt'], nu_data['phi'], nu_data['eta']
                    )
                    self.channels['el1jet']['data']['gen_met'].append(met_tot)
                    self.channels['el1jet']['data']['gen_met_phi'].append(met_phi)
    
    def get_results(self):
        """Convert collected data to numpy arrays and return."""
        data_to_save = {}
        
        # Channel name mappings for output
        output_names = {
            'dielectron': 'el',
            'electronmuon': 'elmu',
            'dimuon': 'mu',
            'mu1jet': 'mu1jet_mu',
            'el1jet': 'mu1jet_el',
        }
        
        for channel, output_prefix in output_names.items():
            channel_data = self.channels[channel]['data']
            
            # Convert to numpy arrays
            data_to_save[f'{output_prefix}_pts'] = np.array(channel_data['pts']) / 1000.0 if channel_data['pts'] else np.array([])
            data_to_save[f'{output_prefix}_etas'] = np.array(channel_data['etas']) if channel_data['etas'] else np.array([])
            data_to_save[f'{output_prefix}_jet_mults'] = np.array(channel_data['jet_mults']) if channel_data['jet_mults'] else np.array([])
            data_to_save[f'{output_prefix}_centrality'] = np.array(channel_data['centrality']) if channel_data['centrality'] else np.array([])
            
            # Add generator-level MET if available
            if 'gen_met' in channel_data and channel_data['gen_met']:
                data_to_save[f'{output_prefix}_gen_met'] = np.array(channel_data['gen_met'])
                data_to_save[f'{output_prefix}_gen_met_phi'] = np.array(channel_data['gen_met_phi'])
            
            for i in range(1, 6):
                data_to_save[f'{output_prefix}_met{i}_pts'] = np.array(channel_data[f'met{i}_pts']) if channel_data[f'met{i}_pts'] else np.array([])
                data_to_save[f'{output_prefix}_met{i}_etas'] = np.array(channel_data[f'met{i}_etas']) if channel_data[f'met{i}_etas'] else np.array([])
        
        return data_to_save
    
    def print_summary(self):
        """Print event count summary."""
        channel_names = {
            'dielectron': 'Dielectron',
            'electronmuon': 'Electron-Muon',
            'dimuon': 'Dimuon',
            'mu1jet': 'Muon + >=2 jets',
            'el1jet': 'Electron + >=2 jets',
        }
        
        for channel, name in channel_names.items():
            counts = self.channels[channel]['counts']
            print(f"{name} channel: Selected events = {counts[0]}, loose = {counts[1]}, medium = {counts[2]}, tight = {counts[3]}")
        
        # Calculate totals
        total_loose = sum(ch['counts'][1] for ch in self.channels.values())
        total_medium = sum(ch['counts'][2] for ch in self.channels.values())
        total_tight = sum(ch['counts'][3] for ch in self.channels.values())
        
        print(f"\nTotal selected events across all channels (loose) = {total_loose}")
        print(f"Total selected events across all channels (medium) = {total_medium}")
        print(f"Total selected events across all channels (tight) = {total_tight}")
        
        # Print percentages
        def pct(part, whole):
            return (part / whole * 100) if whole > 0 else 0.0
        
        print("\nPercentages of totals (per OP):")
        for channel, name in channel_names.items():
            counts = self.channels[channel]['counts']
            print(f"  {name}: loose% = {pct(counts[1], total_loose):.2f}%, "
                  f"medium% = {pct(counts[2], total_medium):.2f}%, "
                  f"tight% = {pct(counts[3], total_tight):.2f}%")

def process_central_peripheral(data_to_save):
    """Process central/peripheral data for MET distributions."""
    channel_prefixes = ['el', 'elmu', 'mu', 'mu1jet_mu', 'mu1jet_el']
    
    for i in range(1, 6):
        try:
            vals_list = []
            flags_list = []
            
            for prefix in channel_prefixes:
                vals = np.asarray(data_to_save.get(f'{prefix}_met{i}_pts', [])).ravel()
                flags = np.asarray(data_to_save.get(f'{prefix}_centrality', [])).ravel()
                
                if vals.size == 0 and flags.size == 0:
                    continue
                
                # Align lengths
                if vals.size == flags.size:
                    vals_list.append(vals)
                    flags_list.append(flags)
                elif flags.size == 1 and vals.size > 0:
                    flags_list.append(np.full(vals.size, flags[0], dtype=int))
                    vals_list.append(vals)
                elif vals.size == 1 and flags.size > 0:
                    vals_list.append(np.full(flags.size, vals[0]))
                    flags_list.append(flags)
            
            if len(vals_list) == 0:
                data_to_save[f'met{i}_central_vals'] = np.array([])
                data_to_save[f'met{i}_peripheral_vals'] = np.array([])
                continue
            
            all_vals = np.concatenate(vals_list)
            all_flags = np.concatenate(flags_list)
            
            if all_vals.size == 0 or all_flags.size == 0 or all_vals.size != all_flags.size:
                data_to_save[f'met{i}_central_vals'] = np.array([])
                data_to_save[f'met{i}_peripheral_vals'] = np.array([])
                continue
            
            data_to_save[f'met{i}_central_vals'] = all_vals[all_flags == 1]
            data_to_save[f'met{i}_peripheral_vals'] = all_vals[all_flags == 0]
            
        except Exception as e:
            print(f"Warning: Failed to process central/peripheral data for met{i}: {e}")
            data_to_save[f'met{i}_central_vals'] = np.array([])
            data_to_save[f'met{i}_peripheral_vals'] = np.array([])

# ========== MAIN DIFFERENCES FROM ORIGINAL CODE ==========

# BEFORE: 5+ separate loops over events (once per channel)
# for i in range(el_charge.shape[0]):  # dielectron loop
# for i in range(el_charge.shape[0]):  # electron-muon loop  
# for i in range(mu_charge.shape[0]):  # dimuon loop
# etc...

# AFTER: Single loop processes all channels
# for i in range(n_events):
#     selector.process_dielectron(...)
#     selector.process_electronmuon(...)
#     selector.process_dimuon(...)
#     selector.process_mu1jet(...)
#     selector.process_el1jet(...)

# BEFORE: Repeated code like this 5+ times:
# high_pt_mask = el_pts > 15000
# high_pt_el_pts = el_pts[high_pt_mask]
# high_pt_el_etas = el_etas[high_pt_mask]
# high_pt_el_charges = el_charges[high_pt_mask]

# AFTER: Single reusable function:
# apply_pt_eta_cuts(pts, etas, charges, pt_cut=15000)

# BEFORE: Manual checking like:
# if len(el_loose_i) >= 2 and el_loose_i[0] and el_loose_i[1]:

# AFTER: Reusable function:
# check_quality(quality_arr, n_required)

def main():
    """Main analysis function."""
    file_mc = "output_mcv4.root"
    tree = "reco"
    tree_nu = "particleLevel"
    
    # Load all data once
    print("Loading data from ROOT file...")
    variables = {
        'el': ['charge', 'pt_NOSYS', 'eta', 'phi', 'select_loose_HI_NOSYS', 'select_medium_HI_NOSYS', 'select_tight_HI_NOSYS'],
        'mu': ['charge', 'pt_NOSYS', 'eta', 'phi', 'select_loose_NOSYS', 'select_medium_NOSYS', 'select_tight_NOSYS'],
        'jet': ['pt_NOSYS', 'eta', 'phi'],
        'met': ['1_phi', '1_tot', '2_phi', '2_tot', '3_phi', '3_tot', '4_phi', '4_tot', '5_phi', '5_tot'],
        'other': ['fcal_sum_et'],
        'PL_nu': ['pt', 'phi', 'eta']
    }
    
    data = {}
    for particle, vars in variables.items():
        if particle == 'met':
            for var in vars:
                data[f'met{var}'] = load_data(file_mc, tree, f'met{var}')
        elif particle == 'other':
            for var in vars:
                data[var] = load_data(file_mc, tree, var)
        elif particle == 'PL_nu':
            for var in vars:
                data[f'PL_nu_{var}'] = load_data(file_mc, tree_nu, f'PL_nu_{var}')
        else:
            for var in vars:
                key = f'{particle}_{var}'
                data[key] = load_data(file_mc, tree, key)
    
    # Calculate centrality once
    fcal = data['fcal_sum_et']
    centrality_flag = np.where(fcal > np.median(fcal), 1, 0)
    
    # Initialize selector
    selector = ChannelSelector()
    
    # *** SINGLE LOOP OVER ALL EVENTS ***
    print("Processing events...")
    n_events = len(data['el_charge'])
    
    for i in range(n_events):
        # Prepare electron data
        el_data = {
            'pts': data['el_pt_NOSYS'][i],
            'etas': data['el_eta'][i],
            'charges': data['el_charge'][i],
            'loose': data['el_select_loose_HI_NOSYS'][i],
            'medium': data['el_select_medium_HI_NOSYS'][i],
            'tight': data['el_select_tight_HI_NOSYS'][i],
        }
        
        # Prepare muon data
        mu_data = {
            'pts': data['mu_pt_NOSYS'][i],
            'etas': data['mu_eta'][i],
            'charges': data['mu_charge'][i],
            'loose': data['mu_select_loose_NOSYS'][i],
            'medium': data['mu_select_medium_NOSYS'][i],
            'tight': data['mu_select_tight_NOSYS'][i],
        }
        
        # Prepare MET data
        met_data = {
            'tots': [data[f'met{j}_tot'][i] for j in range(1, 6)],
            'phis': [data[f'met{j}_phi'][i] for j in range(1, 6)],
        }
        nu_data = {
            'pt': data['PL_nu_pt'][i],
            'phi': data['PL_nu_phi'][i],
            'eta': data['PL_nu_eta'][i],
        }
        
        jet_pts = data['jet_pt_NOSYS'][i]
        cent_flag = centrality_flag[i]
        
        # *** Process all 5 channels in single pass ***
        selector.process_dielectron(el_data, jet_pts, met_data, cent_flag, nu_data)
        selector.process_electronmuon(el_data, mu_data, jet_pts, met_data, cent_flag, nu_data)
        selector.process_dimuon(mu_data, jet_pts, met_data, cent_flag, nu_data)
        selector.process_mu1jet(mu_data, jet_pts, met_data, cent_flag, nu_data)
        selector.process_el1jet(el_data, jet_pts, met_data, cent_flag, nu_data)
    
    # Get results
    data_to_save = selector.get_results()
    
    # Process central/peripheral data
    process_central_peripheral(data_to_save)
    
    # Print summary
    selector.print_summary()
    
    # Save results
    output_filename = 'analysis_output.npz'
    try:
        np.savez_compressed(output_filename, **data_to_save)
        print(f"\nSuccessfully saved analysis data to {output_filename}")
    except Exception as e:
        print(f"\nError: Failed to save data to {output_filename}: {e}")

if __name__ == "__main__":
    main()