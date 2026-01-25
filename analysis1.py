import ROOT
import numpy as np

def load_data(file_path, tree_name, variable_name):
    """Load a single variable from a ROOT file's TTree as a 2D array (events x objects)."""
    file = ROOT.TFile(file_path)
    tree = file.Get(tree_name)
    data = []
    for event in tree:
        value = getattr(event, variable_name)
        # If value is a collection (e.g., std::vector), convert to list
        try:
            arr = list(value)
        except TypeError:
            arr = [value]
        data.append(arr)
    return np.array(data, dtype=object)  # dtype=object for jagged arrays

# --- 1. Load Data ---
file_mc = "output_mcv3.root"
tree = "reco"

el_charge = load_data(file_mc, tree, "el_charge")
el_pt = load_data(file_mc, tree, "el_pt_NOSYS")
el_eta = load_data(file_mc, tree, "el_eta")
el_phi = load_data(file_mc, tree, "el_phi")
mu_charge = load_data(file_mc, tree, "mu_charge")
mu_pt = load_data(file_mc, tree, "mu_pt_NOSYS")
mu_eta = load_data(file_mc, tree, "mu_eta")
mu_phi = load_data(file_mc, tree, "mu_phi")
jet_pt = load_data(file_mc, tree, "jet_pt_NOSYS")
jet_eta = load_data(file_mc, tree, "jet_eta")
jet_phi = load_data(file_mc, tree, "jet_phi")
mu_loose = load_data(file_mc, tree, "mu_select_loose_NOSYS")
el_loose = load_data(file_mc, tree, "el_select_loose_HI_NOSYS")
met1_phi = load_data(file_mc, tree, "met1_phi")
met1_tot = load_data(file_mc, tree, "met1_tot")
met2_phi = load_data(file_mc, tree, "met2_phi")
met2_tot = load_data(file_mc, tree, "met2_tot")
met3_phi = load_data(file_mc, tree, "met3_phi")
met3_tot = load_data(file_mc, tree, "met3_tot")
met4_phi = load_data(file_mc, tree, "met4_phi")
met4_tot = load_data(file_mc, tree, "met4_tot")
met5_phi = load_data(file_mc, tree, "met5_phi")
met5_tot = load_data(file_mc, tree, "met5_tot")
fcal = load_data(file_mc, tree, "fcal_sum_et")
el_medium = load_data(file_mc, tree, "el_select_medium_HI_NOSYS")
el_tight = load_data(file_mc, tree, "el_select_tight_HI_NOSYS")
mu_medium = load_data(file_mc, tree, "mu_select_medium_NOSYS")
mu_tight = load_data(file_mc, tree, "mu_select_tight_NOSYS")
centrality_flag = np.where(fcal > np.median(fcal), 1, 0)

data_to_save = {}

selected_el_pts = []
selected_el_etas = []
selected_el_jet_mults = []
selected_el_met1_pts = []
selected_el_met2_pts = []
selected_el_met3_pts = []
selected_el_met4_pts = []
selected_el_met5_pts = []
selected_el_met1_etas = []
selected_el_met2_etas = []
selected_el_met3_etas = []
selected_el_met4_etas = []
selected_el_met5_etas = []
n_events_dielectron = [0, 0, 0, 0]  # [total passing kinematic+jet cuts, both loose, both medium, both tight]
selected_el_centrality = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    el_loose_i = np.array(el_loose[i])
    el_medium_i = np.array(el_medium[i])
    el_tight_i = np.array(el_tight[i])


    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_etas = el_etas[high_pt_mask]
    high_pt_el_charges = el_charges[high_pt_mask]

    if len(high_pt_el_pts) == 2 and np.prod(high_pt_el_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult > 0:
            n_events_dielectron[0] += 1
            # If both electrons pass the loose operating point (counting per-event)
            if len(el_loose_i) >= 2 and el_loose_i[0] and el_loose_i[1]:
                n_events_dielectron[1] += 1
                selected_el_pts.extend(high_pt_el_pts.tolist())
                selected_el_etas.extend(high_pt_el_etas.tolist())
                selected_el_jet_mults.append(jet_mult)
                selected_el_met1_pts.append(met1_tot[i])
                selected_el_met2_pts.append(met2_tot[i])
                selected_el_met3_pts.append(met3_tot[i])
                selected_el_met4_pts.append(met4_tot[i])
                selected_el_met5_pts.append(met5_tot[i])
                selected_el_met1_etas.append(met1_phi[i])
                selected_el_met2_etas.append(met2_phi[i])
                selected_el_met3_etas.append(met3_phi[i])
                selected_el_met4_etas.append(met4_phi[i])
                selected_el_met5_etas.append(met5_phi[i])
                selected_el_centrality.append(centrality_flag[i])
            if len(el_medium_i) >= 2 and el_medium_i[0] and el_medium_i[1]:
                n_events_dielectron[2] += 1
            if len(el_tight_i) >= 2 and el_tight_i[0] and el_tight_i[1]:
                n_events_dielectron[3] += 1

data_to_save['el_pts'] = np.array(selected_el_pts) / 1000.0
data_to_save['el_etas'] = np.array(selected_el_etas)
data_to_save['el_jet_mults'] = np.array(selected_el_jet_mults)
data_to_save['el_met1_pts'] = np.array(selected_el_met1_pts)
data_to_save['el_met2_pts'] = np.array(selected_el_met2_pts)
data_to_save['el_met3_pts'] = np.array(selected_el_met3_pts)
data_to_save['el_met4_pts'] = np.array(selected_el_met4_pts)
data_to_save['el_met5_pts'] = np.array(selected_el_met5_pts)
data_to_save['el_met1_etas'] = np.array(selected_el_met1_etas)
data_to_save['el_met2_etas'] = np.array(selected_el_met2_etas)
data_to_save['el_met3_etas'] = np.array(selected_el_met3_etas)
data_to_save['el_met4_etas'] = np.array(selected_el_met4_etas)
data_to_save['el_met5_etas'] = np.array(selected_el_met5_etas)
data_to_save['el_centrality'] = np.array(selected_el_centrality)

# --- 3. Electron-muon channel ---
selected_elmu_pts = []
selected_elmu_etas = []
selected_elmu_jet_mults = []
selected_elmu_met1_pts = []
selected_elmu_met2_pts = []
selected_elmu_met3_pts = []
selected_elmu_met4_pts = []
selected_elmu_met5_pts = []
selected_elmu_met1_etas = []
selected_elmu_met2_etas = []
selected_elmu_met3_etas = []
selected_elmu_met4_etas = []
selected_elmu_met5_etas = []
n_events_electronmuon = [0, 0, 0, 0]  # [total, both loose, both medium, both tight]
selected_elmu_centrality = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    el_loose_i = np.array(el_loose[i])
    mu_loose_i = np.array(mu_loose[i])
    el_medium_i = np.array(el_medium[i])
    el_tight_i = np.array(el_tight[i])
    mu_medium_i = np.array(mu_medium[i])
    mu_tight_i = np.array(mu_tight[i])

    high_pt_el_mask = el_pts > 15000
    high_pt_mu_mask = mu_pts > 15000
    high_pt_el_pts = el_pts[high_pt_el_mask]
    high_pt_el_etas = el_etas[high_pt_el_mask]
    high_pt_el_charges = el_charges[high_pt_el_mask]
    high_pt_mu_pts = mu_pts[high_pt_mu_mask]
    high_pt_mu_etas = mu_etas[high_pt_mu_mask]
    high_pt_mu_charges = mu_charges[high_pt_mu_mask]

    if len(high_pt_el_pts) == 1 and len(high_pt_mu_pts) == 1:
        if high_pt_el_charges[0] * high_pt_mu_charges[0] < 0:
            jet_pts = np.array(jet_pt[i])
            jet_mult = np.sum(jet_pts > 30000)
            if jet_mult > 0:
                n_events_electronmuon[0] += 1
                # require both lepton objects to pass given operating points
                if len(el_loose_i) >= 1 and len(mu_loose_i) >= 1 and el_loose_i[0] and mu_loose_i[0]:
                    n_events_electronmuon[1] += 1
                    selected_elmu_pts.extend([high_pt_el_pts[0], high_pt_mu_pts[0]])
                    selected_elmu_etas.extend([high_pt_el_etas[0], high_pt_mu_etas[0]])
                    selected_elmu_jet_mults.append(jet_mult)
                    selected_elmu_met1_pts.append(met1_tot[i])
                    selected_elmu_met2_pts.append(met2_tot[i])
                    selected_elmu_met3_pts.append(met3_tot[i])
                    selected_elmu_met4_pts.append(met4_tot[i])
                    selected_elmu_met5_pts.append(met5_tot[i])
                    selected_elmu_met1_etas.append(met1_phi[i])
                    selected_elmu_met2_etas.append(met2_phi[i])
                    selected_elmu_met3_etas.append(met3_phi[i])
                    selected_elmu_met4_etas.append(met4_phi[i])
                    selected_elmu_met5_etas.append(met5_phi[i])
                    selected_elmu_centrality.append(centrality_flag[i])
                if len(el_medium_i) >= 1 and len(mu_medium_i) >= 1 and el_medium_i[0] and mu_medium_i[0]:
                    n_events_electronmuon[2] += 1
                if len(el_tight_i) >= 1 and len(mu_tight_i) >= 1 and el_tight_i[0] and mu_tight_i[0]:
                    n_events_electronmuon[3] += 1

data_to_save['elmu_pts'] = np.array(selected_elmu_pts) / 1000.0
data_to_save['elmu_etas'] = np.array(selected_elmu_etas)
data_to_save['elmu_jet_mults'] = np.array(selected_elmu_jet_mults)
data_to_save['elmu_met1_pts'] = np.array(selected_elmu_met1_pts)
data_to_save['elmu_met2_pts'] = np.array(selected_elmu_met2_pts)
data_to_save['elmu_met3_pts'] = np.array(selected_elmu_met3_pts)
data_to_save['elmu_met4_pts'] = np.array(selected_elmu_met4_pts)
data_to_save['elmu_met5_pts'] = np.array(selected_elmu_met5_pts)
data_to_save['elmu_met1_etas'] = np.array(selected_elmu_met1_etas)
data_to_save['elmu_met2_etas'] = np.array(selected_elmu_met2_etas)
data_to_save['elmu_met3_etas'] = np.array(selected_elmu_met3_etas)
data_to_save['elmu_met4_etas'] = np.array(selected_elmu_met4_etas)
data_to_save['elmu_met5_etas'] = np.array(selected_elmu_met5_etas)
data_to_save['elmu_centrality'] = np.array(selected_elmu_centrality)

# --- 4. Dimuon channel ---
selected_mu_pts = []
selected_mu_etas = []
selected_mu_jet_mults = []
selected_mu_met1_pts = []
selected_mu_met2_pts = []
selected_mu_met3_pts = []
selected_mu_met4_pts = []
selected_mu_met5_pts = []
selected_mu_met1_etas = []
selected_mu_met2_etas = []
selected_mu_met3_etas = []
selected_mu_met4_etas = []
selected_mu_met5_etas = []
n_events_dimuon = [0, 0, 0, 0]  # [total, both loose, both medium, both tight]
selected_mu_centrality = []

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    mu_loose_i = np.array(mu_loose[i])
    mu_medium_i = np.array(mu_medium[i])
    mu_tight_i = np.array(mu_tight[i])

    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    if len(high_pt_mu_pts) == 2 and np.prod(high_pt_mu_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult > 0:
            n_events_dimuon[0] += 1
            if len(mu_loose_i) >= 2 and mu_loose_i[0] and mu_loose_i[1]:
                n_events_dimuon[1] += 1
                selected_mu_pts.extend(high_pt_mu_pts.tolist())
                selected_mu_etas.extend(high_pt_mu_etas.tolist())
                selected_mu_jet_mults.append(jet_mult)
                selected_mu_met1_pts.append(met1_tot[i])
                selected_mu_met2_pts.append(met2_tot[i])
                selected_mu_met3_pts.append(met3_tot[i])
                selected_mu_met4_pts.append(met4_tot[i])
                selected_mu_met5_pts.append(met5_tot[i])
                selected_mu_met1_etas.append(met1_phi[i])
                selected_mu_met2_etas.append(met2_phi[i])
                selected_mu_met3_etas.append(met3_phi[i])
                selected_mu_met4_etas.append(met4_phi[i])
                selected_mu_met5_etas.append(met5_phi[i])
                selected_mu_centrality.append(centrality_flag[i])
            if len(mu_medium_i) >= 2 and mu_medium_i[0] and mu_medium_i[1]:
                n_events_dimuon[2] += 1
            if len(mu_tight_i) >= 2 and mu_tight_i[0] and mu_tight_i[1]:
                n_events_dimuon[3] += 1

data_to_save['mu_pts'] = np.array(selected_mu_pts) / 1000.0
data_to_save['mu_etas'] = np.array(selected_mu_etas)
data_to_save['mu_jet_mults'] = np.array(selected_mu_jet_mults)
data_to_save['mu_met1_pts'] = np.array(selected_mu_met1_pts)
data_to_save['mu_met2_pts'] = np.array(selected_mu_met2_pts)
data_to_save['mu_met3_pts'] = np.array(selected_mu_met3_pts)
data_to_save['mu_met4_pts'] = np.array(selected_mu_met4_pts)
data_to_save['mu_met5_pts'] = np.array(selected_mu_met5_pts)
data_to_save['mu_met1_etas'] = np.array(selected_mu_met1_etas)
data_to_save['mu_met2_etas'] = np.array(selected_mu_met2_etas)
data_to_save['mu_met3_etas'] = np.array(selected_mu_met3_etas)
data_to_save['mu_met4_etas'] = np.array(selected_mu_met4_etas)
data_to_save['mu_met5_etas'] = np.array(selected_mu_met5_etas)
data_to_save['mu_centrality'] = np.array(selected_mu_centrality)

# --- 5. Muon + >=2 jets channel ---
selected_mu1jet_mu_pts = []
selected_mu1jet_mu_etas = []
selected_mu1jet_mu_jet_mults = []
selected_mu1jet_met1_pts = []
selected_mu1jet_met2_pts = []
selected_mu1jet_met3_pts = []
selected_mu1jet_met4_pts = []
selected_mu1jet_met5_pts = []
selected_mu1jet_met1_etas = []
selected_mu1jet_met2_etas = []
selected_mu1jet_met3_etas = []
selected_mu1jet_met4_etas = []
selected_mu1jet_met5_etas = []
n_events_mul1jet = [0, 0, 0, 0]  # [total, loose, medium, tight]
selected_mu1jet_centrality = [] # Renamed for clarity

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    mu_loose_i = np.array(mu_loose[i])
    mu_medium_i = np.array(mu_medium[i])
    mu_tight_i = np.array(mu_tight[i])

    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]

    if len(high_pt_mu_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            n_events_mul1jet[0] += 1
            if len(mu_loose_i) >= 1 and mu_loose_i[0]:
                n_events_mul1jet[1] += 1
                selected_mu1jet_mu_pts.append(high_pt_mu_pts[0])
                selected_mu1jet_mu_etas.append(high_pt_mu_etas[0])
                selected_mu1jet_mu_jet_mults.append(jet_mult)
                selected_mu1jet_met1_pts.append(met1_tot[i])
                selected_mu1jet_met2_pts.append(met2_tot[i])
                selected_mu1jet_met3_pts.append(met3_tot[i])
                selected_mu1jet_met4_pts.append(met4_tot[i])
                selected_mu1jet_met5_pts.append(met5_tot[i])
                selected_mu1jet_met1_etas.append(met1_phi[i])
                selected_mu1jet_met2_etas.append(met2_phi[i])
                selected_mu1jet_met3_etas.append(met3_phi[i])
                selected_mu1jet_met4_etas.append(met4_phi[i])
                selected_mu1jet_met5_etas.append(met5_phi[i])
                selected_mu1jet_centrality.append(centrality_flag[i])
            if len(mu_medium_i) >= 1 and mu_medium_i[0]:
                n_events_mul1jet[2] += 1
            if len(mu_tight_i) >= 1 and mu_tight_i[0]:
                n_events_mul1jet[3] += 1

data_to_save['mu1jet_mu_pts'] = np.array(selected_mu1jet_mu_pts) / 1000.0
data_to_save['mu1jet_mu_etas'] = np.array(selected_mu1jet_mu_etas)
data_to_save['mu1jet_mu_jet_mults'] = np.array(selected_mu1jet_mu_jet_mults)
data_to_save['mu1jet_met1_pts'] = np.array(selected_mu1jet_met1_pts)
data_to_save['mu1jet_met2_pts'] = np.array(selected_mu1jet_met2_pts)
data_to_save['mu1jet_met3_pts'] = np.array(selected_mu1jet_met3_pts)
data_to_save['mu1jet_met4_pts'] = np.array(selected_mu1jet_met4_pts)
data_to_save['mu1jet_met5_pts'] = np.array(selected_mu1jet_met5_pts)
data_to_save['mu1jet_met1_etas'] = np.array(selected_mu1jet_met1_etas)
data_to_save['mu1jet_met2_etas'] = np.array(selected_mu1jet_met2_etas)
data_to_save['mu1jet_met3_etas'] = np.array(selected_mu1jet_met3_etas)
data_to_save['mu1jet_met4_etas'] = np.array(selected_mu1jet_met4_etas)
data_to_save['mu1jet_met5_etas'] = np.array(selected_mu1jet_met5_etas)
data_to_save['mu1jet_centrality'] = np.array(selected_mu1jet_centrality)

# --- 6. Electron + >=2 jets channel ---
selected_mu1jet_el_pts = []
selected_mu1jet_el_etas = []
selected_mu1jet_el_jet_mults = []
selected_mu1jet_el_met1_pts = []
selected_mu1jet_el_met2_pts = []
selected_mu1jet_el_met3_pts = []
selected_mu1jet_el_met4_pts = []
selected_mu1jet_el_met5_pts = []
selected_mu1jet_el_met1_etas = []
selected_mu1jet_el_met2_etas = []
selected_mu1jet_el_met3_etas = []
selected_mu1jet_el_met4_etas = []
selected_mu1jet_el_met5_etas = []
n_events_el1jet = [0, 0, 0, 0]  # [total, loose, medium, tight]
selected_mu1jet_el_centrality = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    el_loose_i = np.array(el_loose[i])
    el_medium_i = np.array(el_medium[i])
    el_tight_i = np.array(el_tight[i])

    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_etas = el_etas[high_pt_mask]

    if len(high_pt_el_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            n_events_el1jet[0] += 1
            if len(el_loose_i) >= 1 and el_loose_i[0]:
                n_events_el1jet[1] += 1
                selected_mu1jet_el_pts.append(high_pt_el_pts[0])
                selected_mu1jet_el_etas.append(high_pt_el_etas[0])
                selected_mu1jet_el_jet_mults.append(jet_mult)
                selected_mu1jet_el_met1_pts.append(met1_tot[i])
                selected_mu1jet_el_met2_pts.append(met2_tot[i])
                selected_mu1jet_el_met3_pts.append(met3_tot[i])
                selected_mu1jet_el_met4_pts.append(met4_tot[i])
                selected_mu1jet_el_met5_pts.append(met5_tot[i])
                selected_mu1jet_el_met1_etas.append(met1_phi[i])
                selected_mu1jet_el_met2_etas.append(met2_phi[i])
                selected_mu1jet_el_met3_etas.append(met3_phi[i])
                selected_mu1jet_el_met4_etas.append(met4_phi[i])
                selected_mu1jet_el_met5_etas.append(met5_phi[i])
                selected_mu1jet_el_centrality.append(centrality_flag[i])
            if len(el_medium_i) >= 1 and el_medium_i[0]:
                n_events_el1jet[2] += 1
            if len(el_tight_i) >= 1 and el_tight_i[0]:
                n_events_el1jet[3] += 1

data_to_save['mu1jet_el_pts'] = np.array(selected_mu1jet_el_pts) / 1000.0
data_to_save['mu1jet_el_etas'] = np.array(selected_mu1jet_el_etas)
data_to_save['mu1jet_el_jet_mults'] = np.array(selected_mu1jet_el_jet_mults)
data_to_save['mu1jet_el_met1_pts'] = np.array(selected_mu1jet_el_met1_pts)
data_to_save['mu1jet_el_met2_pts'] = np.array(selected_mu1jet_el_met2_pts)
data_to_save['mu1jet_el_met3_pts'] = np.array(selected_mu1jet_el_met3_pts)
data_to_save['mu1jet_el_met4_pts'] = np.array(selected_mu1jet_el_met4_pts)
data_to_save['mu1jet_el_met5_pts'] = np.array(selected_mu1jet_el_met5_pts)
data_to_save['mu1jet_el_met1_etas'] = np.array(selected_mu1jet_el_met1_etas)
data_to_save['mu1jet_el_met2_etas'] = np.array(selected_mu1jet_el_met2_etas)
data_to_save['mu1jet_el_met3_etas'] = np.array(selected_mu1jet_el_met3_etas)
data_to_save['mu1jet_el_met4_etas'] = np.array(selected_mu1jet_el_met4_etas)
data_to_save['mu1jet_el_met5_etas'] = np.array(selected_mu1jet_el_met5_etas)
data_to_save['mu1jet_el_centrality'] = np.array(selected_mu1jet_el_centrality)


# --- 7. Print event counts including operating points ---
print(f"Dielectron channel: Selected events = {n_events_dielectron[0]}, loose = {n_events_dielectron[1]}, medium = {n_events_dielectron[2]}, tight = {n_events_dielectron[3]}")
print(f"Electron-Muon channel: Selected events = {n_events_electronmuon[0]}, loose = {n_events_electronmuon[1]}, medium = {n_events_electronmuon[2]}, tight = {n_events_electronmuon[3]}")
print(f"Dimuon channel: Selected events = {n_events_dimuon[0]}, loose = {n_events_dimuon[1]}, medium = {n_events_dimuon[2]}, tight = {n_events_dimuon[3]}")
print(f"Muon + >=2 jets channel: Selected events = {n_events_mul1jet[0]}, loose = {n_events_mul1jet[1]}, medium = {n_events_mul1jet[2]}, tight = {n_events_mul1jet[3]}")
print(f"Electron + >=2 jets channel: Selected events = {n_events_el1jet[0]}, loose = {n_events_el1jet[1]}, medium = {n_events_el1jet[2]}, tight = {n_events_el1jet[3]}")

total_loose = (n_events_dielectron[1] + n_events_electronmuon[1] +
               n_events_dimuon[1] + n_events_mul1jet[1] +
               n_events_el1jet[1])
total_medium = (n_events_dielectron[2] + n_events_electronmuon[2] +
                n_events_dimuon[2] + n_events_mul1jet[2] +
                n_events_el1jet[2])
total_tight = (n_events_dielectron[3] + n_events_electronmuon[3] +
               n_events_dimuon[3] + n_events_mul1jet[3] +
               n_events_el1jet[3])

print(f"Total selected events across all channels (loose) = {total_loose}")
print(f"Total selected events across all channels (medium) = {total_medium}")
print(f"Total selected events across all channels (tight) = {total_tight}")

def pct(part, whole):
    return (part / whole * 100) if whole > 0 else 0.0

print("Percentages of totals (per OP):")
print(f" Dielectron: loose% = {pct(n_events_dielectron[1], total_loose):.2f}%, medium% = {pct(n_events_dielectron[2], total_medium):.2f}%, tight% = {pct(n_events_dielectron[3], total_tight):.2f}%")
print(f" Electron-Muon: loose% = {pct(n_events_electronmuon[1], total_loose):.2f}%, medium% = {pct(n_events_electronmuon[2], total_medium):.2f}%, tight% = {pct(n_events_electronmuon[3], total_tight):.2f}%")
print(f" Dimuon: loose% = {pct(n_events_dimuon[1], total_loose):.2f}%, medium% = {pct(n_events_dimuon[2], total_medium):.2f}%, tight% = {pct(n_events_dimuon[3], total_tight):.2f}%")
print(f" Mu+>=2jets: loose% = {pct(n_events_mul1jet[1], total_loose):.2f}%, medium% = {pct(n_events_mul1jet[2], total_medium):.2f}%, tight% = {pct(n_events_mul1jet[3], total_tight):.2f}%")
print(f" El+>=2jets: loose% = {pct(n_events_el1jet[1], total_loose):.2f}%, medium% = {pct(n_events_el1jet[2], total_medium):.2f}%, tight% = {pct(n_events_el1jet[3], total_tight):.2f}%")

# --- 8. Combine Central/Peripheral data for saving ---
# This logic is copied from your original script, but instead of plotting,
# we save the combined arrays.
for i in range(1, 6):
    try:
        vals_list = []
        flags_list = []

        channel_keys = [
            (data_to_save[f'el_met{i}_pts'], data_to_save['el_centrality']),
            (data_to_save[f'elmu_met{i}_pts'], data_to_save['elmu_centrality']),
            (data_to_save[f'mu_met{i}_pts'], data_to_save['mu_centrality']),
            (data_to_save[f'mu1jet_met{i}_pts'], data_to_save['mu1jet_centrality']),
            (data_to_save[f'mu1jet_el_met{i}_pts'], data_to_save['mu1jet_el_centrality']),
        ]
        
        for vals, flags in channel_keys:
            vals = np.asarray(vals).ravel()
            flags = np.asarray(flags).ravel()

            if vals.size == 0 and flags.size == 0:
                continue
            
            # try to align lengths if one is a single value
            if vals.size == flags.size:
                vals_list.append(vals)
                flags_list.append(flags)
            else:
                if flags.size == 1 and vals.size > 0:
                    flags_list.append(np.full(vals.size, flags[0], dtype=int))
                    vals_list.append(vals)
                elif vals.size == 1 and flags.size > 0:
                    vals_list.append(np.full(flags.size, vals[0]))
                    flags_list.append(flags)
                else:
                    continue # lengths incompatible

        if len(vals_list) == 0:
            continue

        all_vals = np.concatenate(vals_list) if any(v.size > 0 for v in vals_list) else np.array([])
        all_flags = np.concatenate(flags_list) if any(f.size > 0 for f in flags_list) else np.array([])

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


# --- 9. Save all data to a compressed file ---
output_filename = 'analysis_output.npz'
try:
    np.savez_compressed(output_filename, **data_to_save)
    print(f"\nSuccessfully saved analysis data to {output_filename}")
except Exception as e:
    print(f"\nError: Failed to save data to {output_filename}: {e}")

