import ROOT
import numpy as np
import matplotlib.pyplot as plt


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

file_mc = "output_mc.root"
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


# Dielectron channel (already present)
selected_el_pts = []

for i in range(el_charge.shape[0]):
    # Find electrons with pt > 15000
    el_pts = np.array(el_pt[i])
    el_charges = np.array(el_charge[i])
    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_charges = el_charges[high_pt_mask]

    # Check for exactly two electrons with pt > 15000 and opposite charge
    if len(high_pt_el_pts) == 2 and np.prod(high_pt_el_charges) < 0:
        # Check for jets with pt > 30000
        jet_pts = np.array(jet_pt[i])
        if np.sum(jet_pts > 30000) > 0:
            # Add both electron pts to the list
            selected_el_pts.extend(high_pt_el_pts.tolist())

selected_el_pts = np.array(selected_el_pts)/1000.0  # Convert to GeV

# Electron-muon channel
selected_elmu_pts = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_charges = np.array(el_charge[i])
    mu_pts = np.array(mu_pt[i])
    mu_charges = np.array(mu_charge[i])

    # Find electrons and muons with pt > 15000
    high_pt_el_mask = el_pts > 15000
    high_pt_mu_mask = mu_pts > 15000
    high_pt_el_pts = el_pts[high_pt_el_mask]
    high_pt_el_charges = el_charges[high_pt_el_mask]
    high_pt_mu_pts = mu_pts[high_pt_mu_mask]
    high_pt_mu_charges = mu_charges[high_pt_mu_mask]

    # Check for exactly one electron and one muon with opposite charge
    if len(high_pt_el_pts) == 1 and len(high_pt_mu_pts) == 1:
        if high_pt_el_charges[0] * high_pt_mu_charges[0] < 0:
            jet_pts = np.array(jet_pt[i])
            if np.sum(jet_pts > 30000) > 0:
                # Add both lepton pts to the list
                selected_elmu_pts.extend([high_pt_el_pts[0], high_pt_mu_pts[0]])

selected_elmu_pts = np.array(selected_elmu_pts)/1000.0  # Convert to GeV

# Dimuon channel
selected_mu_pts = []

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_charges = np.array(mu_charge[i])
    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    # Check for exactly two muons with pt > 15000 and opposite charge
    if len(high_pt_mu_pts) == 2 and np.prod(high_pt_mu_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        if np.sum(jet_pts > 30000) > 0:
            selected_mu_pts.extend(high_pt_mu_pts.tolist())

selected_mu_pts = np.array(selected_mu_pts)/1000.0  # Convert to GeV

selected_mu1jet_pts = []

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_charges = np.array(mu_charge[i])
    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    # Check for exactly one muon with pt > 15000
    if len(high_pt_mu_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        # Require at least 2 jets with pt > 30000
        if np.sum(jet_pts > 30000) >= 2:
            selected_mu1jet_pts.append(high_pt_mu_pts[0])

selected_mu1jet_pts = np.array(selected_mu1jet_pts)/1000.0  # Convert to GeV

# Plot histograms separately

plt.figure()
plt.hist(selected_el_pts, bins=50, range=(0, max(selected_el_pts) * 1.1), histtype='step', color='blue')
plt.xlabel('Electron $p_T$ [GeV]')
plt.ylabel('Events')
plt.title('MC: Dielectron channel: Electron $p_T$')
plt.savefig('MC_dielectron_pt.png')
plt.show()

plt.figure()
plt.hist(selected_elmu_pts, bins=50, range=(0, max(selected_elmu_pts) * 1.1), histtype='step', color='green')
plt.xlabel('Lepton $p_T$ [GeV]')
plt.ylabel('Events')
plt.title('MC: Electron-Muon channel: Lepton $p_T$')
plt.savefig('MC_elmu_pt.png')
plt.show()

plt.figure()
plt.hist(selected_mu_pts, bins=50, range=(0, max(selected_mu_pts) * 1.1), histtype='step', color='red')
plt.xlabel('Muon $p_T$ [GeV]')
plt.ylabel('Events')
plt.title('MC: Dimuon channel: Muon $p_T$')
plt.savefig('MC_dimuon_pt.png')
plt.show()

plt.figure()
plt.hist(selected_mu1jet_pts, bins=50, range=(0, max(selected_mu_pts) * 1.1), histtype='step', color='red')
plt.xlabel('Lepton $p_T$ [GeV]')
plt.ylabel('Events')
plt.title('MC: Lep + jet channel: Lepton $p_T$')
plt.savefig('MC_lep+jet.png')
plt.show()