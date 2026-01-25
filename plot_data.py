import numpy as np
import matplotlib.pyplot as plt

INPUT_FILENAME = 'data_histograms.npz'

def load_data():
    try:
        return np.load(INPUT_FILENAME, allow_pickle=True)
    except FileNotFoundError:
        print(f"Error: {INPUT_FILENAME} not found.")
        exit()

def plot_prebinned(data, keys, labels, xlabel, title, filename, log_y=True):
    plt.figure(figsize=(8, 6))
    colors = ['blue', 'red', 'green', 'purple', 'orange']
    
    for i, key in enumerate(keys):
        try:
            counts = data[f"{key}_counts"]
            edges = data[f"{key}_edges"]
            if np.sum(counts) > 0:
                plt.stairs(counts, edges, label=labels[i], color=colors[i], linewidth=1.5)
        except KeyError:
            pass

    plt.xlabel(xlabel)
    plt.ylabel("Events")
    plt.title(title)
    plt.legend()
    
    if log_y:
        plt.yscale('log')
        plt.ylim(bottom=0.1)
    else:
        plt.ylim(bottom=0)

    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(filename)
    print(f"Saved {filename}")
    plt.close()

def main():
    data = load_data()
    print("Data loaded. Generating plots...")

    keys = ["el_pt", "elmu_pt", "mu_pt", "mu1jet_pt", "el1jet_pt"]
    labels = ["Dielectron", "Electron-Muon", "Dimuon", "Muon + >=2 jets", "Electron + >=2 jets"]

    # 1. Pt
    plot_prebinned(data, keys, labels, 
                   xlabel="Lepton $p_T$ [GeV]", title="Data: Lepton $p_T$", filename="data_lepton_pt.png")

    # 2. Eta
    plot_prebinned(data, [k.replace("_pt", "_eta") for k in keys], labels,
                   xlabel="Lepton $\eta$", title="Data: Lepton $\eta$", filename="data_lepton_eta.png", log_y=False)

    # 3. Jet Multiplicity
    plot_prebinned(data, [k.replace("_pt", "_njet") for k in keys], labels,
                   xlabel="Jet Multiplicity", title="Data: Jet Multiplicity", filename="data_jet_multiplicity.png")

    # 4. FCal (Rescaled 0-10)
    if 'fcal_counts' in data:
        plt.figure(figsize=(8, 6))
        plt.stairs(data['fcal_counts'], data['fcal_edges'], color='gray', fill=True, alpha=0.5)
        plt.xlabel("FCal $\Sigma E_T$ (0-10)")
        plt.ylabel("Events")
        plt.title("Data: FCal Sum Et Distribution")
        plt.yscale('log')
        plt.grid(alpha=0.3)
        plt.tight_layout()
        plt.savefig("data_fcal_sum_et.png")
        print("Saved data_fcal_sum_et.png")

if __name__ == "__main__":
    main()