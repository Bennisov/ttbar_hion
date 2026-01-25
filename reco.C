#define reco_cxx
#include "reco.h"
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLorentzVector.h>
#include <fstream>
#include <iostream>
#include <iomanip>

void reco::Loop()
{
    // 2D arrays for each channel: [electron_WP][muon_WP]
    // WP levels: 0=Loose, 1=Medium, 2=Tight
    int count_dielectron[3][3] = {{0}};
    int count_electronmuon[3][3] = {{0}};
    int count_dimuon[3][3] = {{0}};
    int count_mujet[3][3] = {{0}};
    int count_eljet[3][3] = {{0}};
    
    // Create histograms for MET distributions (Z-mass window: 66-116 GeV)
    TH1F *h_met1 = new TH1F("h_met1", "MET Distributions for Z-mass Window;MET [MeV];Events", 50, 0, 200000);
    TH1F *h_met2 = new TH1F("h_met2", "MET Distributions for Z-mass Window;MET [MeV];Events", 50, 0, 200000);
    TH1F *h_met3 = new TH1F("h_met3", "MET Distributions for Z-mass Window;MET [MeV];Events", 50, 0, 200000);
    TH1F *h_met4 = new TH1F("h_met4", "MET Distributions for Z-mass Window;MET [MeV];Events", 50, 0, 200000);
    TH1F *h_met5 = new TH1F("h_met5", "MET Distributions for Z-mass Window;MET [MeV];Events", 50, 0, 200000);
    TH1F *h_mt_e = new TH1F("h_mt_e", "Transverse Mass Electron;M_{T} [GeV];Events", 50, 0, 200);
    TH1F *h_mt_mu = new TH1F("h_mt_mu", "Transverse Mass Muon;M_{T} [GeV];Events", 50, 0, 200);
    TH1F *h_mt_emu = new TH1F("h_mt_emu", "Transverse Mass Electron-Muon;M_{T} [GeV];Events", 50, 0, 200);
    
    // Set histogram styles
    h_met1->SetLineColor(kRed);
    h_met1->SetLineWidth(2);
    h_met2->SetLineColor(kBlue);
    h_met2->SetLineWidth(2);
    h_met3->SetLineColor(kGreen+2);
    h_met3->SetLineWidth(2);
    h_met4->SetLineColor(kMagenta);
    h_met4->SetLineWidth(2);
    h_met5->SetLineColor(kOrange);
    h_met5->SetLineWidth(2);
    
    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        
        // Loop over electron WP levels
        for (int el_wp = 0; el_wp < 3; el_wp++) {
            // Loop over muon WP levels
            for (int mu_wp = 0; mu_wp < 3; mu_wp++) {
                
                // Select electrons based on WP and pt > 15 GeV
                std::vector<int> good_el_indices;
                for (size_t i = 0; i < el_pt_NOSYS->size(); i++) {
                    if (((*el_pt_NOSYS)[i] > 15000.0) &&
                    (((*el_pt2cone20)[i] / (*el_pt_NOSYS)[i]) < 0.16) &&
                    (((*el_topoetcone20)[i] / (*el_pt_NOSYS)[i]) < 0.22))
                    {
                        bool pass_wp = false;
                        if (el_wp == 0 && (*el_select_loose_HI_NOSYS)[i] != 0) pass_wp = true;
                        else if (el_wp == 1 && (*el_select_medium_HI_NOSYS)[i] != 0) pass_wp = true;
                        else if (el_wp == 2 && (*el_select_tight_HI_NOSYS)[i] != 0) pass_wp = true;
                        
                        if (pass_wp) {
                            good_el_indices.push_back(i);
                        }
                    }
                }
                
                // Select muons based on WP and pt > 15 GeV
                std::vector<int> good_mu_indices;
                for (size_t i = 0; i < mu_pt_NOSYS->size(); i++) {
                    if (((*mu_pt_NOSYS)[i] > 15000.0) &&
                    (((*mu_pt2cone20)[i] / (*mu_pt_NOSYS)[i]) < 0.32) &&
                    (((*mu_topoetcone20)[i] / (*mu_pt_NOSYS)[i]) < 0.50))
                    {
                        bool pass_wp = false;
                        if (mu_wp == 0 && (*mu_select_loose_NOSYS)[i] != 0) pass_wp = true;
                        else if (mu_wp == 1 && (*mu_select_medium_NOSYS)[i] != 0) pass_wp = true;
                        else if (mu_wp == 2 && (*mu_select_tight_NOSYS)[i] != 0) pass_wp = true;
                        
                        if (pass_wp) {
                            good_mu_indices.push_back(i);
                        }
                    }
                }
                
                int n_good_el = good_el_indices.size();
                int n_good_mu = good_mu_indices.size();
                
                // Count jets with pt > 30 GeV
                int jet_mult = 0;
                for (size_t i = 0; i < jetR4_pt_NOSYS->size(); i++) {
                    if ((*jetR4_pt_NOSYS)[i] > 30000.0) {
                        jet_mult++;
                    }
                }
                
                // Check dielectron channel
                bool pass_dielectron = false;
                if (n_good_el == 2 && jet_mult > 0) {
                    int idx_0 = good_el_indices[0];
                    int idx_1 = good_el_indices[1];
                    if ((*el_charge)[idx_0] * (*el_charge)[idx_1] < 0) {
                        pass_dielectron = true;
                        
                        // Calculate invariant mass using TLorentzVector
                        TLorentzVector lep1, lep2;
                        lep1.SetPtEtaPhiM((*el_pt_NOSYS)[idx_0], (*el_eta)[idx_0], 
                                         (*el_phi)[idx_0], 0.000511); // electron mass in GeV
                        lep2.SetPtEtaPhiM((*el_pt_NOSYS)[idx_1], (*el_eta)[idx_1], 
                                         (*el_phi)[idx_1], 0.000511);
                        
                        double m_inv = (lep1 + lep2).M() / 1000.0; // Convert MeV to GeV
                        h_mt_e->Fill(m_inv);
                        
                        // Fill MET histograms if in Z-mass window
                        if (m_inv >= 66.0 && m_inv <= 116.0) {
                            h_met1->Fill(met1_tot * 1000.0);
                            h_met2->Fill(met2_tot * 1000.0);
                            h_met3->Fill(met3_tot * 1000.0);
                            h_met4->Fill(met4_tot * 1000.0);
                            h_met5->Fill(met5_tot * 1000.0);
                        }
                    }
                }
                if (pass_dielectron) count_dielectron[el_wp][mu_wp]++;
                
                // Check electron-muon channel
                bool pass_electronmuon = false;
                if (n_good_el == 1 && n_good_mu == 1 && jet_mult > 0) {
                    int idx_el = good_el_indices[0];
                    int idx_mu = good_mu_indices[0];
                    if ((*el_charge)[idx_el] * (*mu_charge)[idx_mu] < 0) {
                        pass_electronmuon = true;
                        
                        // Calculate invariant mass using TLorentzVector
                        TLorentzVector lep_el, lep_mu;
                        lep_el.SetPtEtaPhiM((*el_pt_NOSYS)[idx_el], (*el_eta)[idx_el], 
                                           (*el_phi)[idx_el], 0.000511); // electron mass
                        lep_mu.SetPtEtaPhiM((*mu_pt_NOSYS)[idx_mu], (*mu_eta)[idx_mu], 
                                           (*mu_phi)[idx_mu], 0.10566); // muon mass in GeV
                        
                        double m_inv = (lep_el + lep_mu).M() / 1000.0; // Convert MeV to GeV
                        h_mt_emu->Fill(m_inv);
                        // Fill MET histograms if in Z-mass window
                        if (m_inv >= 66.0 && m_inv <= 116.0) {
                            h_met1->Fill(met1_tot * 1000.0);
                            h_met2->Fill(met2_tot * 1000.0);
                            h_met3->Fill(met3_tot * 1000.0);
                            h_met4->Fill(met4_tot * 1000.0);
                            h_met5->Fill(met5_tot * 1000.0);
                        }
                    }
                }
                if (pass_electronmuon) count_electronmuon[el_wp][mu_wp]++;
                
                // Check dimuon channel
                bool pass_dimuon = false;
                if (n_good_mu == 2 && jet_mult > 0) {
                    int idx_0 = good_mu_indices[0];
                    int idx_1 = good_mu_indices[1];
                    if ((*mu_charge)[idx_0] * (*mu_charge)[idx_1] < 0) {
                        pass_dimuon = true;
                        
                        // Calculate invariant mass using TLorentzVector
                        TLorentzVector lep1, lep2;
                        lep1.SetPtEtaPhiM((*mu_pt_NOSYS)[idx_0], (*mu_eta)[idx_0], 
                                         (*mu_phi)[idx_0], 0.10566); // muon mass in GeV
                        lep2.SetPtEtaPhiM((*mu_pt_NOSYS)[idx_1], (*mu_eta)[idx_1], 
                                         (*mu_phi)[idx_1], 0.10566);
                        
                        double m_inv = (lep1 + lep2).M() / 1000.0; // Convert MeV to GeV
                        h_mt_mu->Fill(m_inv);
                        // Fill MET histograms if in Z-mass window
                        if (m_inv >= 66.0 && m_inv <= 116.0) {
                            h_met1->Fill(met1_tot * 1000.0);
                            h_met2->Fill(met2_tot * 1000.0);
                            h_met3->Fill(met3_tot * 1000.0);
                            h_met4->Fill(met4_tot * 1000.0);
                            h_met5->Fill(met5_tot * 1000.0);
                        }
                    }
                }
                if (pass_dimuon) count_dimuon[el_wp][mu_wp]++;
                
                // Check muon-jet channel
                if (n_good_mu == 1 && jet_mult >= 2) {
                    count_mujet[el_wp][mu_wp]++;
                }
                
                // Check electron-jet channel
                if (n_good_el == 1 && jet_mult >= 2) {
                    count_eljet[el_wp][mu_wp]++;
                }
                
            } // end muon WP loop
        } // end electron WP loop
        
    } // end event loop
    
    // Create canvas and plot MET distributions
    TCanvas *c1 = new TCanvas("c1", "MET Distributions", 800, 600);
    gStyle->SetOptStat(0);
    
    // Find maximum to scale y-axis
    double max_val = h_met1->GetMaximum();
    max_val = std::max(max_val, h_met2->GetMaximum());
    max_val = std::max(max_val, h_met3->GetMaximum());
    max_val = std::max(max_val, h_met4->GetMaximum());
    max_val = std::max(max_val, h_met5->GetMaximum());
    
    h_met1->SetMaximum(max_val * 1.2);
    h_met1->Draw("HIST");
    h_met2->Draw("HIST SAME");
    h_met3->Draw("HIST SAME");
    h_met4->Draw("HIST SAME");
    h_met5->Draw("HIST SAME");
    
    // Add legend
    TLegend *leg = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg->AddEntry(h_met1, "MET1", "l");
    leg->AddEntry(h_met2, "MET2", "l");
    leg->AddEntry(h_met3, "MET3", "l");
    leg->AddEntry(h_met4, "MET4", "l");
    leg->AddEntry(h_met5, "MET5", "l");
    leg->Draw();
    
    c1->SaveAs("MET_distributions_Zmass.pdf");
    c1->SaveAs("MET_distributions_Zmass.png");
    
    std::cout << "MET distribution plots saved!" << std::endl;

    // Create canvas and plot MT distributions
    TCanvas *c2 = new TCanvas("c2", "MT Distributions", 800, 600);
    gStyle->SetOptStat(0);
    
    // Find maximum to scale y-axis
    max_val = h_mt_e->GetMaximum();
    max_val = std::max(max_val, h_mt_mu->GetMaximum());
    max_val = std::max(max_val, h_mt_emu->GetMaximum());
    
    h_mt_e->SetMaximum(max_val * 1.2);
    h_mt_e->Draw("HIST");
    h_mt_mu->Draw("HIST SAME");
    h_mt_emu->Draw("HIST SAME");
    
    // Add legend
    TLegend *leg2 = new TLegend(0.65, 0.65, 0.88, 0.88);
    leg2->AddEntry(h_mt_e, "MT_Electron", "l");
    leg2->AddEntry(h_mt_mu, "MT_Muon", "l");
    leg2->AddEntry(h_mt_emu, "MT_Electron-Muon", "l");
    leg2->Draw();
    
    c2->SaveAs("MT_distributions.pdf");
    c2->SaveAs("MT_distributions.png");
    
    std::cout << "MT distribution plots saved!" << std::endl;
   
    // Write results to output file
    std::ofstream outfile("channel_counts_2D.txt");
    
    outfile << "Channel Analysis Results with Operating Points" << std::endl;
    outfile << "==============================================" << std::endl;
    outfile << "Total events analyzed: " << nentries << std::endl;
    outfile << std::endl;
    outfile << "2D Array Format: [Electron WP][Muon WP]" << std::endl;
    outfile << "Electron WP: 0=HILoose, 1=HIMedium, 2=HITight" << std::endl;
    outfile << "Muon WP: 0=Loose, 1=Medium, 2=Tight" << std::endl;
    outfile << std::endl;
    outfile << "Selection Criteria:" << std::endl;
    outfile << "  - Dielectron: exactly 2 e, opposite charges, pt > 15 GeV, >= 1 jet (pt > 30 GeV)" << std::endl;
    outfile << "  - Electron-Muon: exactly 1 e + 1 mu, opposite charges, pt > 15 GeV, >= 1 jet (pt > 30 GeV)" << std::endl;
    outfile << "  - Dimuon: exactly 2 mu, opposite charges, pt > 15 GeV, >= 1 jet (pt > 30 GeV)" << std::endl;
    outfile << "  - Muon-Jet: exactly 1 mu, pt > 15 GeV, >= 2 jets (pt > 30 GeV)" << std::endl;
    outfile << "  - Electron-Jet: exactly 1 e, pt > 15 GeV, >= 2 jets (pt > 30 GeV)" << std::endl;
    outfile << std::endl;
    outfile << "MET Distributions:" << std::endl;
    outfile << "  - Plotted for dilepton events with invariant mass in Z-window (66-116 GeV)" << std::endl;
    outfile << "  - All three dilepton channels combined (ee, emu, mumu)" << std::endl;
    outfile << std::endl;
    
    // Helper lambda to print 2D array
    auto print_2d_array = [&outfile](const char* channel_name, int arr[3][3]) {
        outfile << "===== " << channel_name << " =====" << std::endl;
        outfile << std::setw(15) << "" << std::setw(12) << "Mu_Loose" 
                << std::setw(12) << "Mu_Medium" << std::setw(12) << "Mu_Tight" << std::endl;
        
        const char* el_labels[3] = {"El_HILoose", "El_HIMedium", "El_HITight"};
        for (int i = 0; i < 3; i++) {
            outfile << std::setw(15) << el_labels[i];
            for (int j = 0; j < 3; j++) {
                outfile << std::setw(12) << arr[i][j];
            }
            outfile << std::endl;
        }
        outfile << std::endl;
    };
    
    print_2d_array("Dielectron Channel", count_dielectron);
    print_2d_array("Electron-Muon Channel", count_electronmuon);
    print_2d_array("Dimuon Channel", count_dimuon);
    print_2d_array("Muon-Jet Channel", count_mujet);
    print_2d_array("Electron-Jet Channel", count_eljet);
    
    // Calculate totals for each WP across all channels
    int total_loose = 0, total_medium = 0, total_tight = 0;
    
    // Sum along diagonal (matching el and mu WPs)
    total_loose = count_dielectron[0][0] + count_electronmuon[0][0] + count_dimuon[0][0] + 
                  count_mujet[0][0] + count_eljet[0][0];
    total_medium = count_dielectron[1][1] + count_electronmuon[1][1] + count_dimuon[1][1] + 
                   count_mujet[1][1] + count_eljet[1][1];
    total_tight = count_dielectron[2][2] + count_electronmuon[2][2] + count_dimuon[2][2] + 
                  count_mujet[2][2] + count_eljet[2][2];
    
    outfile << "===== Summary Statistics (matching el and mu WPs) =====" << std::endl;
    outfile << "Total events across all channels:" << std::endl;
    outfile << "  Loose:  " << total_loose << std::endl;
    outfile << "  Medium: " << total_medium << std::endl;
    outfile << "  Tight:  " << total_tight << std::endl;
    outfile << std::endl;
    
    auto pct = [](int part, int whole) -> double {
        return (whole > 0) ? (part * 100.0 / whole) : 0.0;
    };
    
    outfile << "Channel percentages of totals (per WP):" << std::endl;
    outfile << std::fixed << std::setprecision(2);
    outfile << "  Dielectron:    loose=" << pct(count_dielectron[0][0], total_loose) << "%, "
            << "medium=" << pct(count_dielectron[1][1], total_medium) << "%, "
            << "tight=" << pct(count_dielectron[2][2], total_tight) << "%" << std::endl;
    outfile << "  Electron-Muon: loose=" << pct(count_electronmuon[0][0], total_loose) << "%, "
            << "medium=" << pct(count_electronmuon[1][1], total_medium) << "%, "
            << "tight=" << pct(count_electronmuon[2][2], total_tight) << "%" << std::endl;
    outfile << "  Dimuon:        loose=" << pct(count_dimuon[0][0], total_loose) << "%, "
            << "medium=" << pct(count_dimuon[1][1], total_medium) << "%, "
            << "tight=" << pct(count_dimuon[2][2], total_tight) << "%" << std::endl;
    outfile << "  Muon-Jet:      loose=" << pct(count_mujet[0][0], total_loose) << "%, "
            << "medium=" << pct(count_mujet[1][1], total_medium) << "%, "
            << "tight=" << pct(count_mujet[2][2], total_tight) << "%" << std::endl;
    outfile << "  Electron-Jet:  loose=" << pct(count_eljet[0][0], total_loose) << "%, "
            << "medium=" << pct(count_eljet[1][1], total_medium) << "%, "
            << "tight=" << pct(count_eljet[2][2], total_tight) << "%" << std::endl;
    
    // Add MET histogram statistics to file
    outfile << std::endl;
    outfile << "===== MET Statistics (Z-mass window: 66-116 GeV) =====" << std::endl;
    outfile << "Total events in Z-mass window across all dilepton channels:" << std::endl;
    outfile << "  MET1: " << h_met1->GetEntries() << " events, mean = " 
            << h_met1->GetMean() << " MeV" << std::endl;
    outfile << "  MET2: " << h_met2->GetEntries() << " events, mean = " 
            << h_met2->GetMean() << " MeV" << std::endl;
    outfile << "  MET3: " << h_met3->GetEntries() << " events, mean = " 
            << h_met3->GetMean() << " MeV" << std::endl;
    outfile << "  MET4: " << h_met4->GetEntries() << " events, mean = " 
            << h_met4->GetMean() << " MeV" << std::endl;
    outfile << "  MET5: " << h_met5->GetEntries() << " events, mean = " 
            << h_met5->GetMean() << " MeV" << std::endl;
    
    outfile.close();
    
    // Also print summary to screen
    std::cout << "Channel Analysis with Operating Points completed!" << std::endl;
    std::cout << "Results written to channel_counts_2D.txt" << std::endl;
    std::cout << std::endl;
    std::cout << "Quick Summary (matching el/mu WPs):" << std::endl;
    std::cout << "  Loose  - Total: " << total_loose << std::endl;
    std::cout << "    Dielectron: " << count_dielectron[0][0] << std::endl;
    std::cout << "    Electron-Muon: " << count_electronmuon[0][0] << std::endl;
    std::cout << "    Dimuon: " << count_dimuon[0][0] << std::endl;
    std::cout << "    Muon-Jet: " << count_mujet[0][0] << std::endl;
    std::cout << "    Electron-Jet: " << count_eljet[0][0] << std::endl;
    std::cout << std::endl;
    std::cout << "  Medium - Total: " << total_medium << std::endl;
    std::cout << "    Dielectron: " << count_dielectron[1][1] << std::endl;
    std::cout << "    Electron-Muon: " << count_electronmuon[1][1] << std::endl;
    std::cout << "    Dimuon: " << count_dimuon[1][1] << std::endl;
    std::cout << "    Muon-Jet: " << count_mujet[1][1] << std::endl;
    std::cout << "    Electron-Jet: " << count_eljet[1][1] << std::endl;
    std::cout << std::endl;
    std::cout << "  Tight  - Total: " << total_tight << std::endl;
    std::cout << "    Dielectron: " << count_dielectron[2][2] << std::endl;
    std::cout << "    Electron-Muon: " << count_electronmuon[2][2] << std::endl;
    std::cout << "    Dimuon: " << count_dimuon[2][2] << std::endl;
    std::cout << "    Muon-Jet: " << count_mujet[2][2] << std::endl;
    std::cout << "    Electron-Jet: " << count_eljet[2][2] << std::endl;
    
    std::cout << std::endl;
    std::cout << "Events in Z-mass window (66-116 GeV): " << h_met1->GetEntries() << std::endl;
}