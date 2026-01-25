import ROOT
import numpy as np

# --- Configuration ---
FILE_MC = "output_mcv4.root"
OUTPUT_NPZ = "analysis_output.npz"
OUTPUT_TXT = "counts_mc.txt"

# Cuts
PT_EL_CUT = 15000
PT_MU_CUT = 15000
PT_JET_CUT = 30000
MET_CUT_VAL = 30000  # 30 GeV MET cut

# Enable multi-threading
ROOT.EnableImplicitMT()

def analyze_mc():
    print(f"Initializing RDataFrame on {FILE_MC}...")
    
    # 1. Load Trees and Friend them
    chain = ROOT.TChain("reco")
    chain.Add(FILE_MC)
    chain.AddFriend("particleLevel", FILE_MC) 
    
    df = ROOT.RDataFrame(chain)

    # 2. Define GenMET (C++ Helper)
    ROOT.gInterpreter.Declare("""
    template <typename T>
    float compute_gen_met(const ROOT::RVec<T>& pts, const ROOT::RVec<T>& phis) {
        if (pts.empty()) return 0.0f;
        
        double px = 0.0;
        double py = 0.0;
        
        for(size_t i=0; i < pts.size(); ++i) {
            double pt = (pts[i] > 500.0) ? (pts[i] / 1000.0) : pts[i];
            px += pt * std::cos(phis[i]);
            py += pt * std::sin(phis[i]);
        }
        return (float)std::sqrt(px*px + py*py);
    }
    
    template <typename T>
    float compute_gen_met_phi(const ROOT::RVec<T>& pts, const ROOT::RVec<T>& phis) {
        if (pts.empty()) return 0.0f;
        
        double px = 0.0;
        double py = 0.0;
        
        for(size_t i=0; i < pts.size(); ++i) {
            double pt = (pts[i] > 500.0) ? (pts[i] / 1000.0) : pts[i];
            px += pt * std::cos(phis[i]);
            py += pt * std::sin(phis[i]);
        }
        return (float)std::atan2(py, px);
    }
    """)

    # 3. Define Variables & Masks
    df = df.Define("gen_met", "compute_gen_met(PL_nu_pt, PL_nu_phi)") \
           .Define("gen_met_phi", "compute_gen_met_phi(PL_nu_pt, PL_nu_phi)") \
           .Define("met_reco", "met1_tot") \
           .Define("met_reco_phi", "met1_phi") \
           .Define("el_kin", f"el_pt_NOSYS > {PT_EL_CUT}") \
           .Define("mu_kin", f"mu_pt_NOSYS > {PT_MU_CUT}") \
           .Define("jet_mask", f"jet_pt_NOSYS > {PT_JET_CUT}") \
           .Define("n_good_jets", "Sum(jet_mask)") \
           .Define("el_loose",  "el_kin && el_select_loose_HI_NOSYS") \
           .Define("el_medium", "el_kin && el_select_medium_HI_NOSYS") \
           .Define("el_tight",  "el_kin && el_select_tight_HI_NOSYS") \
           .Define("mu_loose",  "mu_kin && mu_select_loose_NOSYS") \
           .Define("mu_medium", "mu_kin && mu_select_medium_NOSYS") \
           .Define("mu_tight",  "mu_kin && mu_select_tight_NOSYS")

    # 4. Helper to Get MET Distributions After Cuts
    def get_met_after_cuts(el_mask, mu_mask):
        # Dielectron channel
        df_el = df.Filter(f"Sum({el_mask}) == 2", "2 El") \
                  .Filter(f"el_charge[{el_mask}][0] * el_charge[{el_mask}][1] < 0", "OppQ") \
                  .Filter("n_good_jets > 0")
        
        # Dimuon channel (with MET cut)
        df_mu = df.Filter(f"Sum({mu_mask}) == 2", "2 Mu") \
                  .Filter(f"mu_charge[{mu_mask}][0] * mu_charge[{mu_mask}][1] < 0", "OppQ") \
                  .Filter("n_good_jets > 0")

        # Electron-Muon channel (with MET cut)
        df_elmu = df.Filter(f"Sum({el_mask}) == 1 && Sum({mu_mask}) == 1", "1El 1Mu") \
                    .Filter(f"el_charge[{el_mask}][0] * mu_charge[{mu_mask}][0] < 0", "OppQ") \
                    .Filter("n_good_jets > 0")
                   
        # Muon + >=2 Jets
        df_mu1j = df.Filter(f"Sum({mu_mask}) == 1", "1 Mu") \
                    .Filter("n_good_jets >= 2")

        # Electron + >=2 Jets
        df_el1j = df.Filter(f"Sum({el_mask}) == 1", "1 El") \
                    .Filter("n_good_jets >= 2")
        
        # Get counts
        counts = [
            df_el.Count(),
            df_mu.Count(),
            df_elmu.Count(),
            df_mu1j.Count(),
            df_el1j.Count()
        ]
        
        # Get MET distributions for each channel
        met_data = {
            'dielectron': df_el.AsNumpy(columns=["gen_met", "gen_met_phi", "met_reco", "met_reco_phi"]),
            'dimuon': df_mu.AsNumpy(columns=["gen_met", "gen_met_phi", "met_reco", "met_reco_phi"]),
            'elmu': df_elmu.AsNumpy(columns=["gen_met", "gen_met_phi", "met_reco", "met_reco_phi"]),
            'mu1j': df_mu1j.AsNumpy(columns=["gen_met", "gen_met_phi", "met_reco", "met_reco_phi"]),
            'el1j': df_el1j.AsNumpy(columns=["gen_met", "gen_met_phi", "met_reco", "met_reco_phi"])
        }
        
        return counts, met_data

    # 5. Execute Analysis
    print("Analyzing Loose selection...")
    cnts_loose, met_loose = get_met_after_cuts("el_loose", "mu_loose")
    
    # Trigger execution for counts
    vals_l = [c.GetValue() for c in cnts_loose]

    # 6. Diagnostic Output
    print("\n--- Diagnostic: Sample of Post-Cut MET Values ---")
    channels = ["Dielectron", "Dimuon (+MET)", "El-Mu (+MET)", "Muon+>=2Jets", "Electron+>=2Jets"]
    
    for i, ch_name in enumerate(["dielectron", "dimuon", "elmu", "mu1j", "el1j"]):
        data = met_loose[ch_name]
        n_events = len(data["gen_met"])
        print(f"\n{channels[i]}: {n_events} events")
        
        if n_events > 0:
            print(f"  First 3 events:")
            print(f"  {'GenMET':<12} {'GenMET_phi':<12} {'RecoMET':<12} {'RecoMET_phi':<12}")
            for j in range(min(3, n_events)):
                print(f"  {data['gen_met'][j]:<12.2f} {data['gen_met_phi'][j]:<12.3f} "
                      f"{data['met_reco'][j]:<12.2f} {data['met_reco_phi'][j]:<12.3f}")
            
            print(f"  Mean GenMET: {np.mean(data['gen_met']):.2f}")
            print(f"  Mean RecoMET: {np.mean(data['met_reco']):.2f}")

    # 7. Save Counts
    with open(OUTPUT_TXT, "w") as f:
        f.write(f"{'Channel':<20} {'Loose':<10}\n")
        f.write("-" * 35 + "\n")
        for i, ch in enumerate(channels):
            f.write(f"{ch:<20} {vals_l[i]:<10}\n")
    
    print(f"\nCounts saved to {OUTPUT_TXT}")

    # 8. Save NPZ with all channels and angular distributions
    save_dict = {
        'channels': channels,
        'counts_loose': vals_l
    }
    
    # Add MET data for each channel
    for ch_name, ch_label in zip(["dielectron", "dimuon", "elmu", "mu1j", "el1j"], 
                                   ["diel", "dimu", "elmu", "mu1j", "el1j"]):
        data = met_loose[ch_name]
        save_dict[f'{ch_label}_gen_met'] = data['gen_met']
        save_dict[f'{ch_label}_gen_met_phi'] = data['gen_met_phi']
        save_dict[f'{ch_label}_reco_met'] = data['met_reco']
        save_dict[f'{ch_label}_reco_met_phi'] = data['met_reco_phi']
    
    np.savez_compressed(OUTPUT_NPZ, **save_dict)
    print(f"MET distributions (magnitude and phi) saved to {OUTPUT_NPZ}")
    print("\nAnalysis complete!")

if __name__ == "__main__":
    analyze_mc()