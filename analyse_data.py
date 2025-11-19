import ROOT
import numpy as np

# --- Configuration ---
FILE_PATH = "output_data.root"  # Your data file
TREE_NAME = "reco"
OUTPUT_FILE = "data_histograms.npz"

# Enable multi-threading for speed
ROOT.EnableImplicitMT()

def analyze_data():
    print(f"Initializing RDataFrame on {FILE_PATH}...")
    df = ROOT.RDataFrame(TREE_NAME, FILE_PATH)

    # --- 1. Definitions & Cuts (Lazy Evaluation) ---
    
    # Define Masks: Pt > 15 GeV AND Loose Operating Point
    # We use the variable names provided: el_select_loose_HI_NOSYS and mu_select_loose_NOSYS
    df = df.Define("el_mask_good", "el_pt_NOSYS > 15000 && el_select_loose_HI_NOSYS") \
           .Define("mu_mask_good", "mu_pt_NOSYS > 15000 && mu_select_loose_NOSYS") \
           .Define("jet_mask", "jet_pt_NOSYS > 30000")
           
    # Define Good Objects (Apply the masks above)
    df = df.Define("n_good_jets", "Sum(jet_mask)") \
           .Define("good_el_pt_GeV", "el_pt_NOSYS[el_mask_good] / 1000.0") \
           .Define("good_el_eta", "el_eta[el_mask_good]") \
           .Define("good_el_charge", "el_charge[el_mask_good]") \
           .Define("good_mu_pt_GeV", "mu_pt_NOSYS[mu_mask_good] / 1000.0") \
           .Define("good_mu_eta", "mu_eta[mu_mask_good]") \
           .Define("good_mu_charge", "mu_charge[mu_mask_good]")

    # --- 2. Channel Filters ---
    
    # Dielectron: 2 Good Els, Opp Charge, >= 1 Jet
    df_el = df.Filter("Sum(el_mask_good) == 2", "2 Good Electrons") \
              .Filter("good_el_charge[0] * good_el_charge[1] < 0", "El Opp Charge") \
              .Filter("n_good_jets > 0", "Jets > 0 (El)")

    # Dimuon: 2 Good Mus, Opp Charge, >= 1 Jet
    df_mu = df.Filter("Sum(mu_mask_good) == 2", "2 Good Muons") \
              .Filter("good_mu_charge[0] * good_mu_charge[1] < 0", "Mu Opp Charge") \
              .Filter("n_good_jets > 0", "Jets > 0 (Mu)")

    # El-Mu: 1 Good El, 1 Good Mu, Opp Charge, >= 1 Jet
    df_elmu = df.Filter("Sum(el_mask_good) == 1 && Sum(mu_mask_good) == 1", "1 El 1 Mu") \
                .Filter("good_el_charge[0] * good_mu_charge[0] < 0", "ElMu Opp Charge") \
                .Filter("n_good_jets > 0", "Jets > 0 (ElMu)")
                
    # Combine vectors for generic plotting for ElMu channel
    df_elmu = df_elmu.Define("elmu_pt_combined", "Concatenate(good_el_pt_GeV, good_mu_pt_GeV)") \
                     .Define("elmu_eta_combined", "Concatenate(good_el_eta, good_mu_eta)")

    # Mu + Jets: 1 Good Mu, >= 2 Jets
    df_mu1jet = df.Filter("Sum(mu_mask_good) == 1", "1 Muon") \
                  .Filter("n_good_jets >= 2", "Jets >= 2 (Mu)")

    # El + Jets: 1 Good El, >= 2 Jets
    df_el1jet = df.Filter("Sum(el_mask_good) == 1", "1 Electron") \
                  .Filter("n_good_jets >= 2", "Jets >= 2 (El)")

    # --- 3. Book Event Counts ---
    count_total = df.Count()
    count_el = df_el.Count()
    count_mu = df_mu.Count()
    count_elmu = df_elmu.Count()
    count_mu1jet = df_mu1jet.Count()
    count_el1jet = df_el1jet.Count()

    # --- 4. Define Histogram Models ---
    h_pt_model = ("", "", 50, 0, 200)      # Pt: 0-200 GeV
    h_eta_model = ("", "", 30, -3.0, 3.0)  # Eta: -3 to 3
    h_njet_model = ("", "", 10, 0, 10)     # Jets: 0-10
    h_fcal_model = ("", "", 100, 0, 10)    # FCal: 0-10

    print("Booking histograms...")

    # Book Histograms (Pointers)
    h_pt = {
        "el": df_el.Histo1D(h_pt_model, "good_el_pt_GeV"),
        "mu": df_mu.Histo1D(h_pt_model, "good_mu_pt_GeV"),
        "elmu": df_elmu.Histo1D(h_pt_model, "elmu_pt_combined"),
        "mu1jet": df_mu1jet.Histo1D(h_pt_model, "good_mu_pt_GeV"),
        "el1jet": df_el1jet.Histo1D(h_pt_model, "good_el_pt_GeV")
    }

    h_eta = {
        "el": df_el.Histo1D(h_eta_model, "good_el_eta"),
        "mu": df_mu.Histo1D(h_eta_model, "good_mu_eta"),
        "elmu": df_elmu.Histo1D(h_eta_model, "elmu_eta_combined"),
        "mu1jet": df_mu1jet.Histo1D(h_eta_model, "good_mu_eta"),
        "el1jet": df_el1jet.Histo1D(h_eta_model, "good_el_eta")
    }

    h_njet = {
        "el": df_el.Histo1D(h_njet_model, "n_good_jets"),
        "mu": df_mu.Histo1D(h_njet_model, "n_good_jets"),
        "elmu": df_elmu.Histo1D(h_njet_model, "n_good_jets"),
        "mu1jet": df_mu1jet.Histo1D(h_njet_model, "n_good_jets"),
        "el1jet": df_el1jet.Histo1D(h_njet_model, "n_good_jets")
    }

    h_fcal = df.Histo1D(h_fcal_model, "fcal_sum_et")

    # --- 5. Execution ---
    print("Processing events... (This performs the main loop)")
    
    # Force execution by asking for the total count
    total_events = count_total.GetValue()
    
    # --- 6. Print Counts ---
    print("-" * 50)
    print(f"Total Input Events:    {total_events}")
    print("-" * 50)
    print(f"Dielectron Selected:   {count_el.GetValue()}")
    print(f"Dimuon Selected:       {count_mu.GetValue()}")
    print(f"Electron-Muon Selected:{count_elmu.GetValue()}")
    print(f"Muon + >=2 Jets:       {count_mu1jet.GetValue()}")
    print(f"Electron + >=2 Jets:   {count_el1jet.GetValue()}")
    print("-" * 50)

    # --- 7. Extract & Save ---
    results = {}

    def extract_histo(h_ptr):
        h = h_ptr.GetValue()
        n_bins = h.GetNbinsX()
        counts = np.array([h.GetBinContent(i) for i in range(1, n_bins + 1)])
        edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, n_bins + 2)])
        return counts, edges

    for key in h_pt:
        c, e = extract_histo(h_pt[key])
        results[f"{key}_pt_counts"], results[f"{key}_pt_edges"] = c, e
        
        c, e = extract_histo(h_eta[key])
        results[f"{key}_eta_counts"], results[f"{key}_eta_edges"] = c, e
        
        c, e = extract_histo(h_njet[key])
        results[f"{key}_njet_counts"], results[f"{key}_njet_edges"] = c, e

    # Save FCal
    c, e = extract_histo(h_fcal)
    results["fcal_counts"] = c
    results["fcal_edges"] = e

    print(f"Saving histograms to {OUTPUT_FILE}...")
    np.savez_compressed(OUTPUT_FILE, **results)
    print("Done.")

if __name__ == "__main__":
    analyze_data()