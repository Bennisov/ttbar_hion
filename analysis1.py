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
def plot_histogram(data, bins, xlabel, ylabel, title, filename, color='blue'):
    plt.figure()
    plt.hist(data, bins=bins, range=(0, max(data) * 1.1 if len(data) > 0 else 1), histtype='step', color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(filename)
    plt.show()

def plot_histogram_eta(data, bins, xlabel, ylabel, title, filename, color='blue'):
    plt.figure()
    plt.hist(data, bins=bins, range=(-3, 3), histtype='step', color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(filename)
    plt.show()

def plot_histogram_int(data, bins, xlabel, ylabel, title, filename, color='blue'):
    plt.figure()
    plt.hist(data, bins=bins, range=(0, max(data) + 1 if len(data) > 0 else 1), histtype='step', color=color)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.savefig(filename)
    plt.show()

# Dielectron channel (already present)
selected_el_pts = []
selected_el_etas = []
selected_el_jet_mults = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_etas = el_etas[high_pt_mask]
    high_pt_el_charges = el_charges[high_pt_mask]

    if len(high_pt_el_pts) == 2 and np.prod(high_pt_el_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult > 0:
            selected_el_pts.extend(high_pt_el_pts.tolist())
            selected_el_etas.extend(high_pt_el_etas.tolist())
            selected_el_jet_mults.extend([jet_mult])

selected_el_pts = np.array(selected_el_pts)/1000.0
selected_el_etas = np.array(selected_el_etas)
selected_el_jet_mults = np.array(selected_el_jet_mults)

# Electron-muon channel
selected_elmu_pts = []
selected_elmu_etas = []
selected_elmu_jet_mults = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])

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
                selected_elmu_pts.extend([high_pt_el_pts[0], high_pt_mu_pts[0]])
                selected_elmu_etas.extend([high_pt_el_etas[0], high_pt_mu_etas[0]])
                selected_elmu_jet_mults.extend([jet_mult])

selected_elmu_pts = np.array(selected_elmu_pts)/1000.0
selected_elmu_etas = np.array(selected_elmu_etas)
selected_elmu_jet_mults = np.array(selected_elmu_jet_mults)

# Dimuon channel
selected_mu_pts = []
selected_mu_etas = []
selected_mu_jet_mults = []

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    if len(high_pt_mu_pts) == 2 and np.prod(high_pt_mu_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult > 0:
            selected_mu_pts.extend(high_pt_mu_pts.tolist())
            selected_mu_etas.extend(high_pt_mu_etas.tolist())
            selected_mu_jet_mults.extend([jet_mult])

selected_mu_pts = np.array(selected_mu_pts)/1000.0
selected_mu_etas = np.array(selected_mu_etas)
selected_mu_jet_mults = np.array(selected_mu_jet_mults)

# Muon + >=2 jets channel
selected_mu1jet_mu_pts = []
selected_mu1jet_mu_etas = []
selected_mu1jet_mu_jet_mults = []

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    if len(high_pt_mu_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            selected_mu1jet_mu_pts.append(high_pt_mu_pts[0])
            selected_mu1jet_mu_etas.append(high_pt_mu_etas[0])
            selected_mu1jet_mu_jet_mults.append(jet_mult)

selected_mu1jet_mu_pts = np.array(selected_mu1jet_mu_pts)/1000.0
selected_mu1jet_mu_etas = np.array(selected_mu1jet_mu_etas)
selected_mu1jet_mu_jet_mults = np.array(selected_mu1jet_mu_jet_mults)

# Electron + >=2 jets channel
selected_mu1jet_el_pts = []
selected_mu1jet_el_etas = []
selected_mu1jet_el_jet_mults = []

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_etas = el_etas[high_pt_mask]
    high_pt_el_charges = el_charges[high_pt_mask]

    if len(high_pt_el_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            selected_mu1jet_el_pts.append(high_pt_el_pts[0])
            selected_mu1jet_el_etas.append(high_pt_el_etas[0])
            selected_mu1jet_el_jet_mults.append(jet_mult)

selected_mu1jet_el_pts = np.array(selected_mu1jet_el_pts)/1000.0
selected_mu1jet_el_etas = np.array(selected_mu1jet_el_etas)
selected_mu1jet_el_jet_mults = np.array(selected_mu1jet_el_jet_mults)

# Plotting
plot_histogram(selected_el_pts, bins=50, xlabel='Electron $p_T$ [GeV]', ylabel='Events',
               title='MC: Dielectron channel: Electron $p_T$', filename='MC_dielectron_pt.png', color='blue')
plot_histogram_eta(selected_el_etas, bins=50, xlabel='Electron $\eta$', ylabel='Events',
                   title='MC: Dielectron channel: Electron $\eta$', filename='MC_dielectron_eta.png', color='blue')
plot_histogram_int(selected_el_jet_mults, bins=10, xlabel='Jet multiplicity', ylabel='Events',
                   title='MC: Dielectron channel: Jet multiplicity', filename='MC_dielectron_jetmult.png', color='blue')

plot_histogram(selected_elmu_pts, bins=50, xlabel='Lepton $p_T$ [GeV]', ylabel='Events',
               title='MC: Electron-Muon channel: Lepton $p_T$', filename='MC_elmu_pt.png', color='green')
plot_histogram_eta(selected_elmu_etas, bins=50, xlabel='Lepton $\eta$', ylabel='Events',
                   title='MC: Electron-Muon channel: Lepton $\eta$', filename='MC_elmu_eta.png', color='green')
plot_histogram_int(selected_elmu_jet_mults, bins=10, xlabel='Jet multiplicity', ylabel='Events',
                   title='MC: Electron-Muon channel: Jet multiplicity', filename='MC_elmu_jetmult.png', color='green')

plot_histogram(selected_mu_pts, bins=50, xlabel='Muon $p_T$ [GeV]', ylabel='Events',
               title='MC: Dimuon channel: Muon $p_T$', filename='MC_dimuon_pt.png', color='red')
plot_histogram_eta(selected_mu_etas, bins=50, xlabel='Muon $\eta$', ylabel='Events',
                   title='MC: Dimuon channel: Muon $\eta$', filename='MC_dimuon_eta.png', color='red')
plot_histogram_int(selected_mu_jet_mults, bins=10, xlabel='Jet multiplicity', ylabel='Events',
                   title='MC: Dimuon channel: Jet multiplicity', filename='MC_dimuon_jetmult.png', color='red')

plot_histogram(selected_mu1jet_mu_pts, bins=50, xlabel='Muon $p_T$ [GeV]', ylabel='Events',
               title='MC: Muon + jet channel: Muon $p_T$', filename='MC_mu+jet_pt.png', color='red')
plot_histogram_eta(selected_mu1jet_mu_etas, bins=50, xlabel='Muon $\eta$', ylabel='Events',
                   title='MC: Muon + jet channel: Muon $\eta$', filename='MC_mu+jet_eta.png', color='red')
plot_histogram_int(selected_mu1jet_mu_jet_mults, bins=10, xlabel='Jet multiplicity', ylabel='Events',
                   title='MC: Muon + jet channel: Jet multiplicity', filename='MC_mu+jet_jetmult.png', color='red')

plot_histogram(selected_mu1jet_el_pts, bins=50, xlabel='Electron $p_T$ [GeV]', ylabel='Events',
               title='MC: Electron + jet channel: Electron $p_T$', filename='MC_el+jet_pt.png', color='red')
plot_histogram_eta(selected_mu1jet_el_etas, bins=50, xlabel='Electron $\eta$', ylabel='Events',
                   title='MC: Electron + jet channel: Electron $\eta$', filename='MC_el+jet_eta.png', color='red')
plot_histogram_int(selected_mu1jet_el_jet_mults, bins=10, xlabel='Jet multiplicity', ylabel='Events',
                   title='MC: Electron + jet channel: Jet multiplicity', filename='MC_el+jet_jetmult.png', color='red')