import ROOT
import numpy as np

# --- Configuration ---
FILE_PATH = "output_data1.root"
TREE_NAME = "reco"
OUTPUT_TXT = "counts_data.txt"
OUTPUT_NPZ = "data_analysis_output.npz"
PT_EL_CUT = 15000
PT_MU_CUT = 15000
PT_JET_CUT = 30000
MET_CUT_VAL = 30000

ROOT.EnableImplicitMT()

def analyze_data():
    print(f"Initializing RDataFrame on {FILE_PATH}...")
    df = ROOT.RDataFrame(TREE_NAME, FILE_PATH)
    
    # 1. Define Variables & Masks
    df = df.Define("met_reco", "met1_tot") \
           .Define("met_reco_phi", "met1_phi") \
           .Define("el_kin", f"el_pt_NOSYS > {PT_EL_CUT}") \
           .Define("mu_kin", f"mu_pt_NOSYS > {PT_MU_CUT}") \
           .Define("jet_mask", f"jetR2_pt_NOSYS > {PT_JET_CUT}") \
           .Define("n_good_jets", "Sum(jet_mask)") \
           .Define("el_loose",  "el_kin && el_select_loose_HI_NOSYS") \
           .Define("el_medium", "el_kin && el_select_medium_HI_NOSYS") \
           .Define("el_tight",  "el_kin && el_select_tight_HI_NOSYS") \
           .Define("mu_loose",  "mu_kin && mu_select_loose_NOSYS") \
           .Define("mu_medium", "mu_kin && mu_select_medium_NOSYS") \
           .Define("mu_tight",  "mu_kin && mu_select_tight_NOSYS") \
           .Define("pass_met_cut", f"met_reco > {MET_CUT_VAL}")
    
    # 2. Define C++ helper functions for extracting lepton kinematics
    ROOT.gInterpreter.Declare("""
    template <typename T>
    float get_lead_pt(const ROOT::RVec<T>& pts, const ROOT::RVec<bool>& mask) {
        auto selected = pts[mask];
        return selected.size() > 0 ? selected[0] : 0.0f;
    }
    
    template <typename T>
    float get_lead_eta(const ROOT::RVec<T>& etas, const ROOT::RVec<bool>& mask) {
        auto selected = etas[mask];
        return selected.size() > 0 ? selected[0] : 0.0f;
    }
    
    template <typename T>
    float get_lead_phi(const ROOT::RVec<T>& phis, const ROOT::RVec<bool>& mask) {
        auto selected = phis[mask];
        return selected.size() > 0 ? selected[0] : 0.0f;
    }
    
    template <typename T>
    float get_sublead_pt(const ROOT::RVec<T>& pts, const ROOT::RVec<bool>& mask) {
        auto selected = pts[mask];
        return selected.size() > 1 ? selected[1] : 0.0f;
    }
    
    template <typename T>
    float get_sublead_eta(const ROOT::RVec<T>& etas, const ROOT::RVec<bool>& mask) {
        auto selected = etas[mask];
        return selected.size() > 1 ? selected[1] : 0.0f;
    }
    
    template <typename T>
    float get_sublead_phi(const ROOT::RVec<T>& phis, const ROOT::RVec<bool>& mask) {
        auto selected = phis[mask];
        return selected.size() > 1 ? selected[1] : 0.0f;
    }
    """)
    
    # 3. Helper to Extract Kinematics After Cuts
    def get_kinematics(el_mask, mu_mask):
        # Dielectron channel
        df_el = df.Filter(f"Sum({el_mask}) == 2", "2 El") \
                  .Filter(f"el_charge[{el_mask}][0] * el_charge[{el_mask}][1] < 0", "OppQ") \
                  .Filter("n_good_jets > 0") \
                  .Define("el_lead_pt", f"get_lead_pt(el_pt_NOSYS, {el_mask})") \
                  .Define("el_lead_eta", f"get_lead_eta(el_eta_NOSYS, {el_mask})") \
                  .Define("el_lead_phi", f"get_lead_phi(el_phi_NOSYS, {el_mask})") \
                  .Define("el_sublead_pt", f"get_sublead_pt(el_pt_NOSYS, {el_mask})") \
                  .Define("el_sublead_eta", f"get_sublead_eta(el_eta_NOSYS, {el_mask})") \
                  .Define("el_sublead_phi", f"get_sublead_phi(el_phi_NOSYS, {el_mask})")
        
        # Dimuon channel (with MET cut)
        df_mu = df.Filter(f"Sum({mu_mask}) == 2", "2 Mu") \
                  .Filter(f"mu_charge[{mu_mask}][0] * mu_charge[{mu_mask}][1] < 0", "OppQ") \
                  .Filter("pass_met_cut", "MET Cut") \
                  .Filter("n_good_jets > 0") \
                  .Define("mu_lead_pt", f"get_lead_pt(mu_pt_NOSYS, {mu_mask})") \
                  .Define("mu_lead_eta", f"get_lead_eta(mu_eta_NOSYS, {mu_mask})") \
                  .Define("mu_lead_phi", f"get_lead_phi(mu_phi_NOSYS, {mu_mask})") \
                  .Define("mu_sublead_pt", f"get_sublead_pt(mu_pt_NOSYS, {mu_mask})") \
                  .Define("mu_sublead_eta", f"get_sublead_eta(mu_eta_NOSYS, {mu_mask})") \
                  .Define("mu_sublead_phi", f"get_sublead_phi(mu_phi_NOSYS, {mu_mask})")
        
        # Electron-Muon channel (with MET cut)
        df_elmu = df.Filter(f"Sum({el_mask}) == 1 && Sum({mu_mask}) == 1", "1El 1Mu") \
                    .Filter(f"el_charge[{el_mask}][0] * mu_charge[{mu_mask}][0] < 0", "OppQ") \
                    .Filter("pass_met_cut", "MET Cut") \
                    .Filter("n_good_jets > 0") \
                    .Define("el_pt", f"get_lead_pt(el_pt_NOSYS, {el_mask})") \
                    .Define("el_eta", f"get_lead_eta(el_eta_NOSYS, {el_mask})") \
                    .Define("el_phi", f"get_lead_phi(el_phi_NOSYS, {el_mask})") \
                    .Define("mu_pt", f"get_lead_pt(mu_pt_NOSYS, {mu_mask})") \
                    .Define("mu_eta", f"get_lead_eta(mu_eta_NOSYS, {mu_mask})") \
                    .Define("mu_phi", f"get_lead_phi(mu_phi_NOSYS, {mu_mask})")
        
        # Muon + >=2 Jets
        df_mu1j = df.Filter(f"Sum({mu_mask}) == 1", "1 Mu") \
                    .Filter("n_good_jets >= 2") \
                    .Define("mu_pt", f"get_lead_pt(mu_pt_NOSYS, {mu_mask})") \
                    .Define("mu_eta", f"get_lead_eta(mu_eta_NOSYS, {mu_mask})") \
                    .Define("mu_phi", f"get_lead_phi(mu_phi_NOSYS, {mu_mask})")
        
        # Electron + >=2 Jets
        df_el1j = df.Filter(f"Sum({el_mask}) == 1", "1 El") \
                    .Filter("n_good_jets >= 2") \
                    .Define("el_pt", f"get_lead_pt(el_pt_NOSYS, {el_mask})") \
                    .Define("el_eta", f"get_lead_eta(el_eta_NOSYS, {el_mask})") \
                    .Define("el_phi", f"get_lead_phi(el_phi_NOSYS, {el_mask})")
        
        # Get counts
        counts = [
            df_el.Count(),
            df_mu.Count(),
            df_elmu.Count(),
            df_mu1j.Count(),
            df_el1j.Count()
        ]
        
        # Extract kinematics for each channel
        kinematics = {}
        
        # Dielectron
        kinematics['dielectron'] = df_el.AsNumpy(columns=[
            "met_reco", "met_reco_phi", 
            "el_lead_pt_sel", "el_lead_eta_sel", "el_lead_phi_sel",
            "el_sublead_pt_sel", "el_sublead_eta_sel", "el_sublead_phi_sel"
        ])
        
        # Dimuon
        kinematics['dimuon'] = df_mu.AsNumpy(columns=[
            "met_reco", "met_reco_phi",
            "mu_lead_pt_sel", "mu_lead_eta_sel", "mu_lead_phi_sel",
            "mu_sublead_pt_sel", "mu_sublead_eta_sel", "mu_sublead_phi_sel"
        ])
        
        # Electron-Muon
        kinematics['elmu'] = df_elmu.AsNumpy(columns=[
            "met_reco", "met_reco_phi",
            "el_pt_sel", "el_eta_sel", "el_phi_sel",
            "mu_pt_sel", "mu_eta_sel", "mu_phi_sel"
        ])
        
        # Muon + >=2 Jets
        kinematics['mu1j'] = df_mu1j.AsNumpy(columns=[
            "met_reco", "met_reco_phi",
            "mu_pt_sel", "mu_eta_sel", "mu_phi_sel"
        ])
        
        # Electron + >=2 Jets
        kinematics['el1j'] = df_el1j.AsNumpy(columns=[
            "met_reco", "met_reco_phi",
            "el_pt_sel", "el_eta_sel", "el_phi_sel"
        ])
        
        return counts, kinematics
    
    # 4. Execute for all Operating Points
    print("Analyzing Loose selection...")
    cnts_l, kin_l = get_kinematics("el_loose", "mu_loose")
    
    print("Analyzing Medium selection...")
    cnts_m, kin_m = get_kinematics("el_medium", "mu_medium")
    
    print("Analyzing Tight selection...")
    cnts_t, kin_t = get_kinematics("el_tight", "mu_tight")
    
    # Trigger execution
    vals_l = [c.GetValue() for c in cnts_l]
    vals_m = [c.GetValue() for c in cnts_m]
    vals_t = [c.GetValue() for c in cnts_t]
    
    # 5. Print Diagnostics
    print("\n--- Diagnostic: Sample Data (Loose Selection) ---")
    channels_info = ["Dielectron", "Dimuon (+MET)", "Electron-Muon (+MET)", 
                     "Muon+>=2Jets", "Electron+>=2Jets"]
    
    for i, (ch_key, ch_name) in enumerate([("dielectron", "Dielectron"), 
                                            ("dimuon", "Dimuon"),
                                            ("elmu", "Electron-Muon"),
                                            ("mu1j", "Muon+>=2Jets"),
                                            ("el1j", "Electron+>=2Jets")]):
        data = kin_l[ch_key]
        n_events = len(data["met_reco"])
        print(f"\n{ch_name}: {n_events} events")
        
        if n_events > 0:
            print(f"  MET mean: {np.mean(data['met_reco'])/1000:.2f} GeV")
            
            if 'el_lead_pt_sel' in data:
                print(f"  Lead electron pT mean: {np.mean(data['el_lead_pt_sel'])/1000:.2f} GeV")
                sublead = data['el_sublead_pt_sel'][data['el_sublead_pt_sel'] > 0]
                if len(sublead) > 0:
                    print(f"  Sublead electron pT mean: {np.mean(sublead)/1000:.2f} GeV")
            elif 'mu_lead_pt_sel' in data:
                print(f"  Lead muon pT mean: {np.mean(data['mu_lead_pt_sel'])/1000:.2f} GeV")
                sublead = data['mu_sublead_pt_sel'][data['mu_sublead_pt_sel'] > 0]
                if len(sublead) > 0:
                    print(f"  Sublead muon pT mean: {np.mean(sublead)/1000:.2f} GeV")
            elif ch_key == 'elmu':
                print(f"  Electron pT mean: {np.mean(data['el_pt_sel'])/1000:.2f} GeV")
                print(f"  Muon pT mean: {np.mean(data['mu_pt_sel'])/1000:.2f} GeV")
            elif 'mu_pt_sel' in data:
                print(f"  Muon pT mean: {np.mean(data['mu_pt_sel'])/1000:.2f} GeV")
            elif 'el_pt_sel' in data:
                print(f"  Electron pT mean: {np.mean(data['el_pt_sel'])/1000:.2f} GeV")
    
    # 6. Save Counts to Text
    with open(OUTPUT_TXT, "w") as f:
        f.write(f"{'Channel':<20} {'Loose':<10} {'Medium':<10} {'Tight':<10}\n")
        f.write("-" * 55 + "\n")
        for i, ch in enumerate(channels_info):
            f.write(f"{ch:<20} {vals_l[i]:<10} {vals_m[i]:<10} {vals_t[i]:<10}\n")
    
    print(f"\nCounts saved to {OUTPUT_TXT}")
    
    # 7. Save NPZ with all kinematics (Loose selection)
    save_dict = {
        'channels': channels_info,
        'counts_loose': vals_l,
        'counts_medium': vals_m,
        'counts_tight': vals_t
    }
    
    # Add kinematics data for loose selection
    for ch_key, ch_label in [("dielectron", "diel"), ("dimuon", "dimu"), 
                              ("elmu", "elmu"), ("mu1j", "mu1j"), ("el1j", "el1j")]:
        data = kin_l[ch_key]
        for key, value in data.items():
            save_dict[f'{ch_label}_{key}'] = value
    
    np.savez_compressed(OUTPUT_NPZ, **save_dict)
    print(f"Kinematics and MET distributions saved to {OUTPUT_NPZ}")
    print("\nData Analysis complete!")

if __name__ == "__main__":
    analyze_data()