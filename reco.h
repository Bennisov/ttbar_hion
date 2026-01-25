//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan 13 17:38:57 2026 by ROOT version 6.30/04
// from TTree reco/xAOD->NTuple tree
// found on file: output_data1.root
//////////////////////////////////////////////////////////

#ifndef reco_h
#define reco_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
#include "vector"
#include "vector"
#include "vector"

class reco {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   vector<unsigned char> *bootstrapWeights;
   vector<float>   *el_charge;
   vector<float>   *el_eta;
   vector<float>   *el_phi;
   vector<float>   *el_pt2cone20;
   vector<float>   *el_pt2cone30;
   vector<float>   *el_pt2cone40;
   vector<float>   *el_ptcone20;
   vector<float>   *el_ptcone30;
   vector<float>   *el_ptvarcone20;
   vector<float>   *el_ptvarcone30;
   vector<float>   *el_topoetcone20;
   vector<float>   *el_topoetcone30;
   vector<float>   *el_topoetcone40;
   ULong64_t       eventNumber;
   Float_t         fcal_sum_et;
   vector<float>   *jetR2_eta;
   vector<float>   *jetR2_phi;
   vector<float>   *jetR4_eta;
   vector<float>   *jetR4_phi;
   UInt_t          mcChannelNumber;
   Float_t         met1_phi;
   Float_t         met1_tot;
   Float_t         met2_phi;
   Float_t         met2_tot;
   Float_t         met3_phi;
   Float_t         met3_tot;
   Float_t         met4_phi;
   Float_t         met4_tot;
   Float_t         met5_phi;
   Float_t         met5_tot;
   vector<float>   *mu_charge;
   vector<float>   *mu_eta;
   vector<float>   *mu_phi;
   vector<float>   *mu_pt2cone20;
   vector<float>   *mu_pt2cone30;
   vector<float>   *mu_pt2cone40;
   vector<float>   *mu_ptcone20;
   vector<float>   *mu_ptcone30;
   vector<float>   *mu_ptcone40;
   vector<float>   *mu_ptvarcone20;
   vector<float>   *mu_ptvarcone30;
   vector<float>   *mu_ptvarcone40;
   vector<float>   *mu_topoetcone20;
   vector<float>   *mu_topoetcone30;
   vector<float>   *mu_topoetcone40;
   UInt_t          runNumber;
   vector<unsigned char> *track_hiloose_PixelHits;
   vector<unsigned char> *track_hiloose_SCTHits;
   vector<int>     *track_hiloose_cutflow;
   vector<float>   *track_hiloose_d0;
   vector<float>   *track_hiloose_eta;
   Int_t           track_hiloose_n;
   vector<float>   *track_hiloose_phi;
   vector<float>   *track_hiloose_pt;
   vector<float>   *track_hiloose_qOverP;
   vector<float>   *track_hiloose_theta;
   vector<float>   *track_hiloose_z0;
   vector<unsigned char> *track_hitight_PixelHits;
   vector<unsigned char> *track_hitight_SCTHits;
   vector<int>     *track_hitight_cutflow;
   vector<float>   *track_hitight_d0;
   vector<float>   *track_hitight_eta;
   Int_t           track_hitight_n;
   vector<float>   *track_hitight_phi;
   vector<float>   *track_hitight_pt;
   vector<float>   *track_hitight_qOverP;
   vector<float>   *track_hitight_theta;
   vector<float>   *track_hitight_z0;
   Bool_t          trigPassed_HLT_e120_lhloose;
   Bool_t          trigPassed_HLT_e140_lhloose_L1EM22VHI;
   Bool_t          trigPassed_HLT_e140_lhloose_nod0;
   Bool_t          trigPassed_HLT_e15_lhloose_ion_L1EM12;
   Bool_t          trigPassed_HLT_e15_loose_nogsf_ion_L1eEM15;
   Bool_t          trigPassed_HLT_e24_lhmedium_L1EM20VH;
   Bool_t          trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI;
   Bool_t          trigPassed_HLT_e26_lhtight_nod0_ivarloose;
   Bool_t          trigPassed_HLT_e60_lhmedium;
   Bool_t          trigPassed_HLT_e60_lhmedium_L1EM22VHI;
   Bool_t          trigPassed_HLT_e60_lhmedium_nod0;
   Bool_t          trigPassed_HLT_mu20_iloose_L1MU15;
   Bool_t          trigPassed_HLT_mu24_ivarmedium_L1MU14FCH;
   Bool_t          trigPassed_HLT_mu26_ivarmedium;
   Bool_t          trigPassed_HLT_mu40;
   Bool_t          trigPassed_HLT_mu4_L1MU3V;
   Bool_t          trigPassed_HLT_mu50;
   Bool_t          trigPassed_HLT_mu50_L1MU14FCH;
   Bool_t          trigPassed_HLT_mu8;
   vector<float>   *el_e_NOSYS;
   vector<float>   *el_pt_NOSYS;
   vector<char>    *el_select_loose_HI_NOSYS;
   vector<char>    *el_select_medium_HI_NOSYS;
   vector<char>    *el_select_outputSelect_NOSYS;
   vector<char>    *el_select_tight_HI_NOSYS;
   vector<float>   *jetR2_e_NOSYS;
   vector<float>   *jetR2_pt_NOSYS;
   vector<char>    *jetR2_select_outputSelect_NOSYS;
   vector<float>   *jetR4_e_NOSYS;
   vector<float>   *jetR4_pt_NOSYS;
   vector<char>    *jetR4_select_outputSelect_NOSYS;
   vector<float>   *mu_e_NOSYS;
   vector<float>   *mu_pt_NOSYS;
   vector<char>    *mu_select_loose_NOSYS;
   vector<char>    *mu_select_medium_NOSYS;
   vector<char>    *mu_select_outputSelect_NOSYS;
   vector<char>    *mu_select_tight_NOSYS;
   Char_t          pass_SUBcommon_NOSYS;
   Char_t          pass_ee_NOSYS;
   Char_t          pass_ejet_NOSYS;
   Char_t          pass_emu_NOSYS;
   Char_t          pass_mujet_NOSYS;
   Char_t          pass_mumu_NOSYS;

   // List of branches
   TBranch        *b_bootstrapWeights;   //!
   TBranch        *b_el_charge;   //!
   TBranch        *b_el_eta;   //!
   TBranch        *b_el_phi;   //!
   TBranch        *b_el_pt2cone20;   //!
   TBranch        *b_el_pt2cone30;   //!
   TBranch        *b_el_pt2cone40;   //!
   TBranch        *b_el_ptcone20;   //!
   TBranch        *b_el_ptcone30;   //!
   TBranch        *b_el_ptvarcone20;   //!
   TBranch        *b_el_ptvarcone30;   //!
   TBranch        *b_el_topoetcone20;   //!
   TBranch        *b_el_topoetcone30;   //!
   TBranch        *b_el_topoetcone40;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_fcal_sum_et;   //!
   TBranch        *b_jetR2_eta;   //!
   TBranch        *b_jetR2_phi;   //!
   TBranch        *b_jetR4_eta;   //!
   TBranch        *b_jetR4_phi;   //!
   TBranch        *b_mcChannelNumber;   //!
   TBranch        *b_met1_phi;   //!
   TBranch        *b_met1_tot;   //!
   TBranch        *b_met2_phi;   //!
   TBranch        *b_met2_tot;   //!
   TBranch        *b_met3_phi;   //!
   TBranch        *b_met3_tot;   //!
   TBranch        *b_met4_phi;   //!
   TBranch        *b_met4_tot;   //!
   TBranch        *b_met5_phi;   //!
   TBranch        *b_met5_tot;   //!
   TBranch        *b_mu_charge;   //!
   TBranch        *b_mu_eta;   //!
   TBranch        *b_mu_phi;   //!
   TBranch        *b_mu_pt2cone20;   //!
   TBranch        *b_mu_pt2cone30;   //!
   TBranch        *b_mu_pt2cone40;   //!
   TBranch        *b_mu_ptcone20;   //!
   TBranch        *b_mu_ptcone30;   //!
   TBranch        *b_mu_ptcone40;   //!
   TBranch        *b_mu_ptvarcone20;   //!
   TBranch        *b_mu_ptvarcone30;   //!
   TBranch        *b_mu_ptvarcone40;   //!
   TBranch        *b_mu_topoetcone20;   //!
   TBranch        *b_mu_topoetcone30;   //!
   TBranch        *b_mu_topoetcone40;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_track_hiloose_PixelHits;   //!
   TBranch        *b_track_hiloose_SCTHits;   //!
   TBranch        *b_track_hiloose_cutflow;   //!
   TBranch        *b_track_hiloose_d0;   //!
   TBranch        *b_track_hiloose_eta;   //!
   TBranch        *b_track_hiloose_n;   //!
   TBranch        *b_track_hiloose_phi;   //!
   TBranch        *b_track_hiloose_pt;   //!
   TBranch        *b_track_hiloose_qOverP;   //!
   TBranch        *b_track_hiloose_theta;   //!
   TBranch        *b_track_hiloose_z0;   //!
   TBranch        *b_track_hitight_PixelHits;   //!
   TBranch        *b_track_hitight_SCTHits;   //!
   TBranch        *b_track_hitight_cutflow;   //!
   TBranch        *b_track_hitight_d0;   //!
   TBranch        *b_track_hitight_eta;   //!
   TBranch        *b_track_hitight_n;   //!
   TBranch        *b_track_hitight_phi;   //!
   TBranch        *b_track_hitight_pt;   //!
   TBranch        *b_track_hitight_qOverP;   //!
   TBranch        *b_track_hitight_theta;   //!
   TBranch        *b_track_hitight_z0;   //!
   TBranch        *b_trigPassed_HLT_e120_lhloose;   //!
   TBranch        *b_trigPassed_HLT_e140_lhloose_L1EM22VHI;   //!
   TBranch        *b_trigPassed_HLT_e140_lhloose_nod0;   //!
   TBranch        *b_trigPassed_HLT_e15_lhloose_ion_L1EM12;   //!
   TBranch        *b_trigPassed_HLT_e15_loose_nogsf_ion_L1eEM15;   //!
   TBranch        *b_trigPassed_HLT_e24_lhmedium_L1EM20VH;   //!
   TBranch        *b_trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI;   //!
   TBranch        *b_trigPassed_HLT_e26_lhtight_nod0_ivarloose;   //!
   TBranch        *b_trigPassed_HLT_e60_lhmedium;   //!
   TBranch        *b_trigPassed_HLT_e60_lhmedium_L1EM22VHI;   //!
   TBranch        *b_trigPassed_HLT_e60_lhmedium_nod0;   //!
   TBranch        *b_trigPassed_HLT_mu20_iloose_L1MU15;   //!
   TBranch        *b_trigPassed_HLT_mu24_ivarmedium_L1MU14FCH;   //!
   TBranch        *b_trigPassed_HLT_mu26_ivarmedium;   //!
   TBranch        *b_trigPassed_HLT_mu40;   //!
   TBranch        *b_trigPassed_HLT_mu4_L1MU3V;   //!
   TBranch        *b_trigPassed_HLT_mu50;   //!
   TBranch        *b_trigPassed_HLT_mu50_L1MU14FCH;   //!
   TBranch        *b_trigPassed_HLT_mu8;   //!
   TBranch        *b_el_e_NOSYS;   //!
   TBranch        *b_el_pt_NOSYS;   //!
   TBranch        *b_el_select_loose_HI_NOSYS;   //!
   TBranch        *b_el_select_medium_HI_NOSYS;   //!
   TBranch        *b_el_select_outputSelect_NOSYS;   //!
   TBranch        *b_el_select_tight_HI_NOSYS;   //!
   TBranch        *b_jetR2_e_NOSYS;   //!
   TBranch        *b_jetR2_pt_NOSYS;   //!
   TBranch        *b_jetR2_select_outputSelect_NOSYS;   //!
   TBranch        *b_jetR4_e_NOSYS;   //!
   TBranch        *b_jetR4_pt_NOSYS;   //!
   TBranch        *b_jetR4_select_outputSelect_NOSYS;   //!
   TBranch        *b_mu_e_NOSYS;   //!
   TBranch        *b_mu_pt_NOSYS;   //!
   TBranch        *b_mu_select_loose_NOSYS;   //!
   TBranch        *b_mu_select_medium_NOSYS;   //!
   TBranch        *b_mu_select_outputSelect_NOSYS;   //!
   TBranch        *b_mu_select_tight_NOSYS;   //!
   TBranch        *b_pass_SUBcommon_NOSYS;   //!
   TBranch        *b_pass_ee_NOSYS;   //!
   TBranch        *b_pass_ejet_NOSYS;   //!
   TBranch        *b_pass_emu_NOSYS;   //!
   TBranch        *b_pass_mujet_NOSYS;   //!
   TBranch        *b_pass_mumu_NOSYS;   //!

   reco(TTree *tree=0);
   virtual ~reco();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef reco_cxx
reco::reco(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output_data1.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("output_data1.root");
      }
      f->GetObject("reco",tree);

   }
   Init(tree);
}

reco::~reco()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t reco::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t reco::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void reco::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   bootstrapWeights = 0;
   el_charge = 0;
   el_eta = 0;
   el_phi = 0;
   el_pt2cone20 = 0;
   el_pt2cone30 = 0;
   el_pt2cone40 = 0;
   el_ptcone20 = 0;
   el_ptcone30 = 0;
   el_ptvarcone20 = 0;
   el_ptvarcone30 = 0;
   el_topoetcone20 = 0;
   el_topoetcone30 = 0;
   el_topoetcone40 = 0;
   jetR2_eta = 0;
   jetR2_phi = 0;
   jetR4_eta = 0;
   jetR4_phi = 0;
   mu_charge = 0;
   mu_eta = 0;
   mu_phi = 0;
   mu_pt2cone20 = 0;
   mu_pt2cone30 = 0;
   mu_pt2cone40 = 0;
   mu_ptcone20 = 0;
   mu_ptcone30 = 0;
   mu_ptcone40 = 0;
   mu_ptvarcone20 = 0;
   mu_ptvarcone30 = 0;
   mu_ptvarcone40 = 0;
   mu_topoetcone20 = 0;
   mu_topoetcone30 = 0;
   mu_topoetcone40 = 0;
   track_hiloose_PixelHits = 0;
   track_hiloose_SCTHits = 0;
   track_hiloose_cutflow = 0;
   track_hiloose_d0 = 0;
   track_hiloose_eta = 0;
   track_hiloose_phi = 0;
   track_hiloose_pt = 0;
   track_hiloose_qOverP = 0;
   track_hiloose_theta = 0;
   track_hiloose_z0 = 0;
   track_hitight_PixelHits = 0;
   track_hitight_SCTHits = 0;
   track_hitight_cutflow = 0;
   track_hitight_d0 = 0;
   track_hitight_eta = 0;
   track_hitight_phi = 0;
   track_hitight_pt = 0;
   track_hitight_qOverP = 0;
   track_hitight_theta = 0;
   track_hitight_z0 = 0;
   el_e_NOSYS = 0;
   el_pt_NOSYS = 0;
   el_select_loose_HI_NOSYS = 0;
   el_select_medium_HI_NOSYS = 0;
   el_select_outputSelect_NOSYS = 0;
   el_select_tight_HI_NOSYS = 0;
   jetR2_e_NOSYS = 0;
   jetR2_pt_NOSYS = 0;
   jetR2_select_outputSelect_NOSYS = 0;
   jetR4_e_NOSYS = 0;
   jetR4_pt_NOSYS = 0;
   jetR4_select_outputSelect_NOSYS = 0;
   mu_e_NOSYS = 0;
   mu_pt_NOSYS = 0;
   mu_select_loose_NOSYS = 0;
   mu_select_medium_NOSYS = 0;
   mu_select_outputSelect_NOSYS = 0;
   mu_select_tight_NOSYS = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("bootstrapWeights", &bootstrapWeights, &b_bootstrapWeights);
   fChain->SetBranchAddress("el_charge", &el_charge, &b_el_charge);
   fChain->SetBranchAddress("el_eta", &el_eta, &b_el_eta);
   fChain->SetBranchAddress("el_phi", &el_phi, &b_el_phi);
   fChain->SetBranchAddress("el_pt2cone20", &el_pt2cone20, &b_el_pt2cone20);
   fChain->SetBranchAddress("el_pt2cone30", &el_pt2cone30, &b_el_pt2cone30);
   fChain->SetBranchAddress("el_pt2cone40", &el_pt2cone40, &b_el_pt2cone40);
   fChain->SetBranchAddress("el_ptcone20", &el_ptcone20, &b_el_ptcone20);
   fChain->SetBranchAddress("el_ptcone30", &el_ptcone30, &b_el_ptcone30);
   fChain->SetBranchAddress("el_ptvarcone20", &el_ptvarcone20, &b_el_ptvarcone20);
   fChain->SetBranchAddress("el_ptvarcone30", &el_ptvarcone30, &b_el_ptvarcone30);
   fChain->SetBranchAddress("el_topoetcone20", &el_topoetcone20, &b_el_topoetcone20);
   fChain->SetBranchAddress("el_topoetcone30", &el_topoetcone30, &b_el_topoetcone30);
   fChain->SetBranchAddress("el_topoetcone40", &el_topoetcone40, &b_el_topoetcone40);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("fcal_sum_et", &fcal_sum_et, &b_fcal_sum_et);
   fChain->SetBranchAddress("jetR2_eta", &jetR2_eta, &b_jetR2_eta);
   fChain->SetBranchAddress("jetR2_phi", &jetR2_phi, &b_jetR2_phi);
   fChain->SetBranchAddress("jetR4_eta", &jetR4_eta, &b_jetR4_eta);
   fChain->SetBranchAddress("jetR4_phi", &jetR4_phi, &b_jetR4_phi);
   fChain->SetBranchAddress("mcChannelNumber", &mcChannelNumber, &b_mcChannelNumber);
   fChain->SetBranchAddress("met1_phi", &met1_phi, &b_met1_phi);
   fChain->SetBranchAddress("met1_tot", &met1_tot, &b_met1_tot);
   fChain->SetBranchAddress("met2_phi", &met2_phi, &b_met2_phi);
   fChain->SetBranchAddress("met2_tot", &met2_tot, &b_met2_tot);
   fChain->SetBranchAddress("met3_phi", &met3_phi, &b_met3_phi);
   fChain->SetBranchAddress("met3_tot", &met3_tot, &b_met3_tot);
   fChain->SetBranchAddress("met4_phi", &met4_phi, &b_met4_phi);
   fChain->SetBranchAddress("met4_tot", &met4_tot, &b_met4_tot);
   fChain->SetBranchAddress("met5_phi", &met5_phi, &b_met5_phi);
   fChain->SetBranchAddress("met5_tot", &met5_tot, &b_met5_tot);
   fChain->SetBranchAddress("mu_charge", &mu_charge, &b_mu_charge);
   fChain->SetBranchAddress("mu_eta", &mu_eta, &b_mu_eta);
   fChain->SetBranchAddress("mu_phi", &mu_phi, &b_mu_phi);
   fChain->SetBranchAddress("mu_pt2cone20", &mu_pt2cone20, &b_mu_pt2cone20);
   fChain->SetBranchAddress("mu_pt2cone30", &mu_pt2cone30, &b_mu_pt2cone30);
   fChain->SetBranchAddress("mu_pt2cone40", &mu_pt2cone40, &b_mu_pt2cone40);
   fChain->SetBranchAddress("mu_ptcone20", &mu_ptcone20, &b_mu_ptcone20);
   fChain->SetBranchAddress("mu_ptcone30", &mu_ptcone30, &b_mu_ptcone30);
   fChain->SetBranchAddress("mu_ptcone40", &mu_ptcone40, &b_mu_ptcone40);
   fChain->SetBranchAddress("mu_ptvarcone20", &mu_ptvarcone20, &b_mu_ptvarcone20);
   fChain->SetBranchAddress("mu_ptvarcone30", &mu_ptvarcone30, &b_mu_ptvarcone30);
   fChain->SetBranchAddress("mu_ptvarcone40", &mu_ptvarcone40, &b_mu_ptvarcone40);
   fChain->SetBranchAddress("mu_topoetcone20", &mu_topoetcone20, &b_mu_topoetcone20);
   fChain->SetBranchAddress("mu_topoetcone30", &mu_topoetcone30, &b_mu_topoetcone30);
   fChain->SetBranchAddress("mu_topoetcone40", &mu_topoetcone40, &b_mu_topoetcone40);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("track_hiloose_PixelHits", &track_hiloose_PixelHits, &b_track_hiloose_PixelHits);
   fChain->SetBranchAddress("track_hiloose_SCTHits", &track_hiloose_SCTHits, &b_track_hiloose_SCTHits);
   fChain->SetBranchAddress("track_hiloose_cutflow", &track_hiloose_cutflow, &b_track_hiloose_cutflow);
   fChain->SetBranchAddress("track_hiloose_d0", &track_hiloose_d0, &b_track_hiloose_d0);
   fChain->SetBranchAddress("track_hiloose_eta", &track_hiloose_eta, &b_track_hiloose_eta);
   fChain->SetBranchAddress("track_hiloose_n", &track_hiloose_n, &b_track_hiloose_n);
   fChain->SetBranchAddress("track_hiloose_phi", &track_hiloose_phi, &b_track_hiloose_phi);
   fChain->SetBranchAddress("track_hiloose_pt", &track_hiloose_pt, &b_track_hiloose_pt);
   fChain->SetBranchAddress("track_hiloose_qOverP", &track_hiloose_qOverP, &b_track_hiloose_qOverP);
   fChain->SetBranchAddress("track_hiloose_theta", &track_hiloose_theta, &b_track_hiloose_theta);
   fChain->SetBranchAddress("track_hiloose_z0", &track_hiloose_z0, &b_track_hiloose_z0);
   fChain->SetBranchAddress("track_hitight_PixelHits", &track_hitight_PixelHits, &b_track_hitight_PixelHits);
   fChain->SetBranchAddress("track_hitight_SCTHits", &track_hitight_SCTHits, &b_track_hitight_SCTHits);
   fChain->SetBranchAddress("track_hitight_cutflow", &track_hitight_cutflow, &b_track_hitight_cutflow);
   fChain->SetBranchAddress("track_hitight_d0", &track_hitight_d0, &b_track_hitight_d0);
   fChain->SetBranchAddress("track_hitight_eta", &track_hitight_eta, &b_track_hitight_eta);
   fChain->SetBranchAddress("track_hitight_n", &track_hitight_n, &b_track_hitight_n);
   fChain->SetBranchAddress("track_hitight_phi", &track_hitight_phi, &b_track_hitight_phi);
   fChain->SetBranchAddress("track_hitight_pt", &track_hitight_pt, &b_track_hitight_pt);
   fChain->SetBranchAddress("track_hitight_qOverP", &track_hitight_qOverP, &b_track_hitight_qOverP);
   fChain->SetBranchAddress("track_hitight_theta", &track_hitight_theta, &b_track_hitight_theta);
   fChain->SetBranchAddress("track_hitight_z0", &track_hitight_z0, &b_track_hitight_z0);
   fChain->SetBranchAddress("trigPassed_HLT_e120_lhloose", &trigPassed_HLT_e120_lhloose, &b_trigPassed_HLT_e120_lhloose);
   fChain->SetBranchAddress("trigPassed_HLT_e140_lhloose_L1EM22VHI", &trigPassed_HLT_e140_lhloose_L1EM22VHI, &b_trigPassed_HLT_e140_lhloose_L1EM22VHI);
   fChain->SetBranchAddress("trigPassed_HLT_e140_lhloose_nod0", &trigPassed_HLT_e140_lhloose_nod0, &b_trigPassed_HLT_e140_lhloose_nod0);
   fChain->SetBranchAddress("trigPassed_HLT_e15_lhloose_ion_L1EM12", &trigPassed_HLT_e15_lhloose_ion_L1EM12, &b_trigPassed_HLT_e15_lhloose_ion_L1EM12);
   fChain->SetBranchAddress("trigPassed_HLT_e15_loose_nogsf_ion_L1eEM15", &trigPassed_HLT_e15_loose_nogsf_ion_L1eEM15, &b_trigPassed_HLT_e15_loose_nogsf_ion_L1eEM15);
   fChain->SetBranchAddress("trigPassed_HLT_e24_lhmedium_L1EM20VH", &trigPassed_HLT_e24_lhmedium_L1EM20VH, &b_trigPassed_HLT_e24_lhmedium_L1EM20VH);
   fChain->SetBranchAddress("trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI", &trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI, &b_trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI);
   fChain->SetBranchAddress("trigPassed_HLT_e26_lhtight_nod0_ivarloose", &trigPassed_HLT_e26_lhtight_nod0_ivarloose, &b_trigPassed_HLT_e26_lhtight_nod0_ivarloose);
   fChain->SetBranchAddress("trigPassed_HLT_e60_lhmedium", &trigPassed_HLT_e60_lhmedium, &b_trigPassed_HLT_e60_lhmedium);
   fChain->SetBranchAddress("trigPassed_HLT_e60_lhmedium_L1EM22VHI", &trigPassed_HLT_e60_lhmedium_L1EM22VHI, &b_trigPassed_HLT_e60_lhmedium_L1EM22VHI);
   fChain->SetBranchAddress("trigPassed_HLT_e60_lhmedium_nod0", &trigPassed_HLT_e60_lhmedium_nod0, &b_trigPassed_HLT_e60_lhmedium_nod0);
   fChain->SetBranchAddress("trigPassed_HLT_mu20_iloose_L1MU15", &trigPassed_HLT_mu20_iloose_L1MU15, &b_trigPassed_HLT_mu20_iloose_L1MU15);
   fChain->SetBranchAddress("trigPassed_HLT_mu24_ivarmedium_L1MU14FCH", &trigPassed_HLT_mu24_ivarmedium_L1MU14FCH, &b_trigPassed_HLT_mu24_ivarmedium_L1MU14FCH);
   fChain->SetBranchAddress("trigPassed_HLT_mu26_ivarmedium", &trigPassed_HLT_mu26_ivarmedium, &b_trigPassed_HLT_mu26_ivarmedium);
   fChain->SetBranchAddress("trigPassed_HLT_mu40", &trigPassed_HLT_mu40, &b_trigPassed_HLT_mu40);
   fChain->SetBranchAddress("trigPassed_HLT_mu4_L1MU3V", &trigPassed_HLT_mu4_L1MU3V, &b_trigPassed_HLT_mu4_L1MU3V);
   fChain->SetBranchAddress("trigPassed_HLT_mu50", &trigPassed_HLT_mu50, &b_trigPassed_HLT_mu50);
   fChain->SetBranchAddress("trigPassed_HLT_mu50_L1MU14FCH", &trigPassed_HLT_mu50_L1MU14FCH, &b_trigPassed_HLT_mu50_L1MU14FCH);
   fChain->SetBranchAddress("trigPassed_HLT_mu8", &trigPassed_HLT_mu8, &b_trigPassed_HLT_mu8);
   fChain->SetBranchAddress("el_e_NOSYS", &el_e_NOSYS, &b_el_e_NOSYS);
   fChain->SetBranchAddress("el_pt_NOSYS", &el_pt_NOSYS, &b_el_pt_NOSYS);
   fChain->SetBranchAddress("el_select_loose_HI_NOSYS", &el_select_loose_HI_NOSYS, &b_el_select_loose_HI_NOSYS);
   fChain->SetBranchAddress("el_select_medium_HI_NOSYS", &el_select_medium_HI_NOSYS, &b_el_select_medium_HI_NOSYS);
   fChain->SetBranchAddress("el_select_outputSelect_NOSYS", &el_select_outputSelect_NOSYS, &b_el_select_outputSelect_NOSYS);
   fChain->SetBranchAddress("el_select_tight_HI_NOSYS", &el_select_tight_HI_NOSYS, &b_el_select_tight_HI_NOSYS);
   fChain->SetBranchAddress("jetR2_e_NOSYS", &jetR2_e_NOSYS, &b_jetR2_e_NOSYS);
   fChain->SetBranchAddress("jetR2_pt_NOSYS", &jetR2_pt_NOSYS, &b_jetR2_pt_NOSYS);
   fChain->SetBranchAddress("jetR2_select_outputSelect_NOSYS", &jetR2_select_outputSelect_NOSYS, &b_jetR2_select_outputSelect_NOSYS);
   fChain->SetBranchAddress("jetR4_e_NOSYS", &jetR4_e_NOSYS, &b_jetR4_e_NOSYS);
   fChain->SetBranchAddress("jetR4_pt_NOSYS", &jetR4_pt_NOSYS, &b_jetR4_pt_NOSYS);
   fChain->SetBranchAddress("jetR4_select_outputSelect_NOSYS", &jetR4_select_outputSelect_NOSYS, &b_jetR4_select_outputSelect_NOSYS);
   fChain->SetBranchAddress("mu_e_NOSYS", &mu_e_NOSYS, &b_mu_e_NOSYS);
   fChain->SetBranchAddress("mu_pt_NOSYS", &mu_pt_NOSYS, &b_mu_pt_NOSYS);
   fChain->SetBranchAddress("mu_select_loose_NOSYS", &mu_select_loose_NOSYS, &b_mu_select_loose_NOSYS);
   fChain->SetBranchAddress("mu_select_medium_NOSYS", &mu_select_medium_NOSYS, &b_mu_select_medium_NOSYS);
   fChain->SetBranchAddress("mu_select_outputSelect_NOSYS", &mu_select_outputSelect_NOSYS, &b_mu_select_outputSelect_NOSYS);
   fChain->SetBranchAddress("mu_select_tight_NOSYS", &mu_select_tight_NOSYS, &b_mu_select_tight_NOSYS);
   fChain->SetBranchAddress("pass_SUBcommon_NOSYS", &pass_SUBcommon_NOSYS, &b_pass_SUBcommon_NOSYS);
   fChain->SetBranchAddress("pass_ee_NOSYS", &pass_ee_NOSYS, &b_pass_ee_NOSYS);
   fChain->SetBranchAddress("pass_ejet_NOSYS", &pass_ejet_NOSYS, &b_pass_ejet_NOSYS);
   fChain->SetBranchAddress("pass_emu_NOSYS", &pass_emu_NOSYS, &b_pass_emu_NOSYS);
   fChain->SetBranchAddress("pass_mujet_NOSYS", &pass_mujet_NOSYS, &b_pass_mujet_NOSYS);
   fChain->SetBranchAddress("pass_mumu_NOSYS", &pass_mumu_NOSYS, &b_pass_mumu_NOSYS);
   Notify();
}

Bool_t reco::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void reco::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t reco::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef reco_cxx
