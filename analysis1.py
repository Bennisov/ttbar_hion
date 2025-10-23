import ROOT
import numpy as np
import atexit
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

file_mc = "output_mcv2.root"
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
el_loose = load_data(file_mc, tree, "el_select_loose_NOSYS")
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



# Dielectron channel (already present)
def plot_five_channel_pts(el_pts, elmu_pts, mu_pts, mu1jet_mu_pts, mu1jet_el_pts,
                          bins=50, xlabel='Lepton $p_T$ [GeV]', ylabel='Events',
                          title='MC: Lepton pT — 5 channels', filename='combined_pt.png'):
    """Overlay 5 channel pT histograms on a single plot."""
    datasets = [
        (np.asarray(el_pts), 'Dielectron', 'blue'),
        (np.asarray(elmu_pts), 'Electron-Muon', 'green'),
        (np.asarray(mu_pts), 'Dimuon', 'red'),
        (np.asarray(mu1jet_mu_pts), 'Muon + >=2 jets', 'purple'),
        (np.asarray(mu1jet_el_pts), 'Electron + >=2 jets', 'orange'),
    ]

    # Determine common range (handle empty datasets)
    max_vals = [d.max() if d.size > 0 else 0 for d, _, _ in datasets]
    max_val = max(max_vals) if max_vals else 1
    if max_val == 0:
        max_val = 1.0
    hist_range = (0.0, max_val * 1.1)

    plt.figure()
    for data, label, color in datasets:
        if data.size == 0:
            continue
        plt.hist(data, bins=bins, range=hist_range, histtype='step', label=label, color=color, linewidth=1.5, density=True)


    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

def _auto_plot_at_exit():
    # Try to call plot_five_channel_pts after the script has populated the selected_* arrays.
    try:
        el = globals().get('selected_el_pts', None)
        elmu = globals().get('selected_elmu_pts', None)
        mu = globals().get('selected_mu_pts', None)
        mu1jet_mu = globals().get('selected_mu1jet_mu_pts', None)
        mu1jet_el = globals().get('selected_mu1jet_el_pts', None)
        if any(v is not None for v in [el, elmu, mu, mu1jet_mu, mu1jet_el]):
            # Ensure we pass arrays (empty if None)
            plot_five_channel_pts(
                el_pts = el if el is not None else np.array([]),
                elmu_pts = elmu if elmu is not None else np.array([]),
                mu_pts = mu if mu is not None else np.array([]),
                mu1jet_mu_pts = mu1jet_mu if mu1jet_mu is not None else np.array([]),
                mu1jet_el_pts = mu1jet_el if mu1jet_el is not None else np.array([])
            )
    except Exception:
        # Fail silently so script can continue even if plotting fails.
        pass

# Register for execution at program exit so the plot is produced after selections are filled.
atexit.register(_auto_plot_at_exit)

# Dielectron channel (updated to use per-event operating point and MET storage)
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
n_events_dielectron = [0, 0]  # [total passing kinematic+jet cuts, passing additionally loose operating point]

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    el_loose_i = np.array(el_loose[i])  # per-event operating-point flags

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
            if el_loose_i[0] and el_loose_i[1]:
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

selected_el_pts = np.array(selected_el_pts) / 1000.0
selected_el_etas = np.array(selected_el_etas)
selected_el_jet_mults = np.array(selected_el_jet_mults)
selected_el_met1_pts = np.array(selected_el_met1_pts)
selected_el_met2_pts = np.array(selected_el_met2_pts)
selected_el_met3_pts = np.array(selected_el_met3_pts)
selected_el_met4_pts = np.array(selected_el_met4_pts)
selected_el_met5_pts = np.array(selected_el_met5_pts)
selected_el_met1_etas = np.array(selected_el_met1_etas)
selected_el_met2_etas = np.array(selected_el_met2_etas)
selected_el_met3_etas = np.array(selected_el_met3_etas)
selected_el_met4_etas = np.array(selected_el_met4_etas)
selected_el_met5_etas = np.array(selected_el_met5_etas)

# Electron-muon channel (now mirrors dielectron: operating-point cut and MET storage)
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
n_events_electronmuon = [0, 0]  # [total passing kinematic+jet cuts, passing additionally loose OP for both]

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    el_loose_i = np.array(el_loose[i])
    mu_loose_i = np.array(mu_loose[i])

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
                # require both lepton objects to pass loose operating point (one electron, one muon)
                if el_loose_i[0] and mu_loose_i[0]:
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

selected_elmu_pts = np.array(selected_elmu_pts) / 1000.0
selected_elmu_etas = np.array(selected_elmu_etas)
selected_elmu_jet_mults = np.array(selected_elmu_jet_mults)
selected_elmu_met1_pts = np.array(selected_elmu_met1_pts)
selected_elmu_met2_pts = np.array(selected_elmu_met2_pts)
selected_elmu_met3_pts = np.array(selected_elmu_met3_pts)
selected_elmu_met4_pts = np.array(selected_elmu_met4_pts)
selected_elmu_met5_pts = np.array(selected_elmu_met5_pts)
selected_elmu_met1_etas = np.array(selected_elmu_met1_etas)
selected_elmu_met2_etas = np.array(selected_elmu_met2_etas)
selected_elmu_met3_etas = np.array(selected_elmu_met3_etas)
selected_elmu_met4_etas = np.array(selected_elmu_met4_etas)
selected_elmu_met5_etas = np.array(selected_elmu_met5_etas)

# Dimuon channel (with operating-point cut and MET storage)
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
n_events_dimuon = [0, 0]  # [total passing kinematic+jet cuts, passing additionally loose OP]

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    mu_loose_i = np.array(mu_loose[i])

    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]
    high_pt_mu_charges = mu_charges[high_pt_mask]

    if len(high_pt_mu_pts) == 2 and np.prod(high_pt_mu_charges) < 0:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult > 0:
            n_events_dimuon[0] += 1
            if mu_loose_i[0] and mu_loose_i[1]:
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

selected_mu_pts = np.array(selected_mu_pts) / 1000.0
selected_mu_etas = np.array(selected_mu_etas)
selected_mu_jet_mults = np.array(selected_mu_jet_mults)
selected_mu_met1_pts = np.array(selected_mu_met1_pts)
selected_mu_met2_pts = np.array(selected_mu_met2_pts)
selected_mu_met3_pts = np.array(selected_mu_met3_pts)
selected_mu_met4_pts = np.array(selected_mu_met4_pts)
selected_mu_met5_pts = np.array(selected_mu_met5_pts)
selected_mu_met1_etas = np.array(selected_mu_met1_etas)
selected_mu_met2_etas = np.array(selected_mu_met2_etas)
selected_mu_met3_etas = np.array(selected_mu_met3_etas)
selected_mu_met4_etas = np.array(selected_mu_met4_etas)
selected_mu_met5_etas = np.array(selected_mu_met5_etas)

# Muon + >=2 jets channel (with operating-point cut and MET storage)
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
n_events_mul1jet = [0, 0]  # [total, passing loose OP]

for i in range(mu_charge.shape[0]):
    mu_pts = np.array(mu_pt[i])
    mu_etas = np.array(mu_eta[i])
    mu_charges = np.array(mu_charge[i])
    mu_loose_i = np.array(mu_loose[i])

    high_pt_mask = mu_pts > 15000
    high_pt_mu_pts = mu_pts[high_pt_mask]
    high_pt_mu_etas = mu_etas[high_pt_mask]

    if len(high_pt_mu_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            n_events_mul1jet[0] += 1
            if mu_loose_i[0]:
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

selected_mu1jet_mu_pts = np.array(selected_mu1jet_mu_pts) / 1000.0
selected_mu1jet_mu_etas = np.array(selected_mu1jet_mu_etas)
selected_mu1jet_mu_jet_mults = np.array(selected_mu1jet_mu_jet_mults)
selected_mu1jet_met1_pts = np.array(selected_mu1jet_met1_pts)
selected_mu1jet_met2_pts = np.array(selected_mu1jet_met2_pts)
selected_mu1jet_met3_pts = np.array(selected_mu1jet_met3_pts)
selected_mu1jet_met4_pts = np.array(selected_mu1jet_met4_pts)
selected_mu1jet_met5_pts = np.array(selected_mu1jet_met5_pts)
selected_mu1jet_met1_etas = np.array(selected_mu1jet_met1_etas)
selected_mu1jet_met2_etas = np.array(selected_mu1jet_met2_etas)
selected_mu1jet_met3_etas = np.array(selected_mu1jet_met3_etas)
selected_mu1jet_met4_etas = np.array(selected_mu1jet_met4_etas)
selected_mu1jet_met5_etas = np.array(selected_mu1jet_met5_etas)

# Electron + >=2 jets channel (with operating-point cut and MET storage)
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
n_events_el1jet = [0, 0]  # [total, passing loose OP]

for i in range(el_charge.shape[0]):
    el_pts = np.array(el_pt[i])
    el_etas = np.array(el_eta[i])
    el_charges = np.array(el_charge[i])
    el_loose_i = np.array(el_loose[i])

    high_pt_mask = el_pts > 15000
    high_pt_el_pts = el_pts[high_pt_mask]
    high_pt_el_etas = el_etas[high_pt_mask]

    if len(high_pt_el_pts) == 1:
        jet_pts = np.array(jet_pt[i])
        jet_mult = np.sum(jet_pts > 30000)
        if jet_mult >= 2:
            n_events_el1jet[0] += 1
            if el_loose_i[0]:
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

selected_mu1jet_el_pts = np.array(selected_mu1jet_el_pts) / 1000.0
selected_mu1jet_el_etas = np.array(selected_mu1jet_el_etas)
selected_mu1jet_el_jet_mults = np.array(selected_mu1jet_el_jet_mults)
selected_mu1jet_el_met1_pts = np.array(selected_mu1jet_el_met1_pts)
selected_mu1jet_el_met2_pts = np.array(selected_mu1jet_el_met2_pts)
selected_mu1jet_el_met3_pts = np.array(selected_mu1jet_el_met3_pts)
selected_mu1jet_el_met4_pts = np.array(selected_mu1jet_el_met4_pts)
selected_mu1jet_el_met5_pts = np.array(selected_mu1jet_el_met5_pts)
selected_mu1jet_el_met1_etas = np.array(selected_mu1jet_el_met1_etas)
selected_mu1jet_el_met2_etas = np.array(selected_mu1jet_el_met2_etas)
selected_mu1jet_el_met3_etas = np.array(selected_mu1jet_el_met3_etas)
selected_mu1jet_el_met4_etas = np.array(selected_mu1jet_el_met4_etas)
selected_mu1jet_el_met5_etas = np.array(selected_mu1jet_el_met5_etas)

print(f"Dielectron channel: Selected events = {n_events_dielectron[0]}, after cuts = {n_events_dielectron[1]}")
print(f"Electron-Muon channel: Selected events = {n_events_electronmuon[0]}, after cuts = {n_events_electronmuon[1]}")
print(f"Dimuon channel: Selected events = {n_events_dimuon[0]}, after cuts = {n_events_dimuon[1]}")
print(f"Muon + >=2 jets channel: Selected events = {n_events_mul1jet[0]}, after cuts = {n_events_mul1jet[1]}")
print(f"Electron + >=2 jets channel: Selected events = {n_events_el1jet[0]}, after cuts = {n_events_el1jet[1]}")

# Combined plotting utilities (use same style as plot_five_channel_pts)
def plot_five_channel_generic(el, elmu, mu, mu1jet_mu, mu1jet_el,
                              bins=50, xlabel='Variable', ylabel='Events',
                              title='MC: Variable — 5 channels', filename='combined.png',
                              hist_range=None, colors=None):
    datasets = [
        (np.asarray(el), 'Dielectron', 'blue'),
        (np.asarray(elmu), 'Electron-Muon', 'green'),
        (np.asarray(mu), 'Dimuon', 'red'),
        (np.asarray(mu1jet_mu), 'Muon + >=2 jets', 'purple'),
        (np.asarray(mu1jet_el), 'Electron + >=2 jets', 'orange'),
    ]
    if colors is not None:
        datasets = [(d, lab, colors[i]) for i, (d, lab, _) in enumerate(datasets)]

    # Determine automatic range if not provided
    if hist_range is None:
        max_vals = [d.max() if d.size > 0 else 0 for d, _, _ in datasets]
        max_val = max(max_vals) if max_vals else 1
        if max_val == 0:
            max_val = 1.0
        hist_range = (0.0, max_val * 1.1)

    plt.figure()
    any_plotted = False
    for data, label, color in datasets:
        if data.size == 0:
            continue
        any_plotted = True
        plt.hist(data, bins=bins, range=hist_range, histtype='step', label=label, color=color, linewidth=1.5, density=True)

    if not any_plotted:
        return  # nothing to plot

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig(filename)
    plt.show()

# 1) Lepton pT (reuse existing function for style consistency)
try:
    plot_five_channel_pts(
        el_pts=selected_el_pts,
        elmu_pts=selected_elmu_pts,
        mu_pts=selected_mu_pts,
        mu1jet_mu_pts=selected_mu1jet_mu_pts,
        mu1jet_el_pts=selected_mu1jet_el_pts,
        bins=50,
        xlabel='Lepton $p_T$ [GeV]',
        title='MC: Lepton $p_T$ — 5 channels',
        filename='combined_pt.png'
    )
except Exception:
    pass

# 2) Lepton eta
plot_five_channel_generic(
    selected_el_etas, selected_elmu_etas, selected_mu_etas,
    selected_mu1jet_mu_etas, selected_mu1jet_el_etas,
    bins=50,
    xlabel='Lepton $\\eta$',
    title='MC: Lepton $\\eta$ — 5 channels',
    filename='combined_eta.png',
    hist_range=(-3.0, 3.0)
)

# 3) Jet multiplicity (integer bins)
plot_five_channel_generic(
    selected_el_jet_mults, selected_elmu_jet_mults, selected_mu_jet_mults,
    selected_mu1jet_mu_jet_mults, selected_mu1jet_el_jet_mults,
    bins=np.arange(0, 11) - 0.5,  # integer bins 0..10
    xlabel='Jet multiplicity',
    title='MC: Jet multiplicity — 5 channels',
    filename='combined_jetmult.png',
    hist_range=( -0.5, 9.5 )
)

# 4) MET totals for met1..met5 (plot totals)
met_labels = ['met1', 'met2', 'met3', 'met4', 'met5']
for i, met_label in enumerate(met_labels, start=1):
    try:
        el_met = globals().get(f'selected_el_met{i}_pts', np.array([]))
        elmu_met = globals().get(f'selected_elmu_met{i}_pts', np.array([]))
        mu_met = globals().get(f'selected_mu_met{i}_pts', np.array([]))
        mu1jet_mu_met = globals().get(f'selected_mu1jet_met{i}_pts', np.array([]))
        mu1jet_el_met = globals().get(f'selected_mu1jet_el_met{i}_pts', np.array([]))

        plot_five_channel_generic(
            el_met, elmu_met, mu_met, mu1jet_mu_met, mu1jet_el_met,
            bins=50,
            xlabel=f'{met_label} total [units]',
            title=f'MC: {met_label} total — 5 channels',
            filename=f'combined_{met_label}_tot.png'
        )
    except Exception:
        pass

# 5) MET phi (angles) for met1..met5
for i, met_label in enumerate(met_labels, start=1):
    try:
        el_met_phi = globals().get(f'selected_el_met{i}_etas', np.array([]))
        elmu_met_phi = globals().get(f'selected_elmu_met{i}_etas', np.array([]))
        mu_met_phi = globals().get(f'selected_mu_met{i}_etas', np.array([]))
        mu1jet_mu_met_phi = globals().get(f'selected_mu1jet_met{i}_etas', np.array([]))
        mu1jet_el_met_phi = globals().get(f'selected_mu1jet_el_met{i}_etas', np.array([]))

        plot_five_channel_generic(
            el_met_phi, elmu_met_phi, mu_met_phi, mu1jet_mu_met_phi, mu1jet_el_met_phi,
            bins=64,
            xlabel=f'{met_label} $\\phi$ [rad]',
            title=f'MC: {met_label} $\\phi$ — 5 channels',
            filename=f'combined_{met_label}_phi.png',
            hist_range=(-np.pi, np.pi)
        )
    except Exception:
        pass