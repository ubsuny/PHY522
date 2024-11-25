// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

// Header file to access Pythia 8 program elements.
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"

// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

int main(int argc, char ** argv) {

  if ( argc < 2 ) {
    std::cout << "usage: " << argv[0] << " root_file" << std::endl;
    return 0;
  }


  // Define the AK8 jet finder.
  double R = 0.8, ptmin = 20.0;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  bool verbose = false;

  char * infile = argv[1];

  double z_cut = 0.10;
  double beta  = 0.0;
  fastjet::contrib::SoftDrop sd(beta, z_cut);

  // HOTVR configuration
  

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open(infile);
  TTree *T = (TTree *)file->Get("T");
  const Int_t kMaxJet = 10;                       // Stores leading 10 jets
  const Int_t kMaxGen = 5000;                     // and 1000 of the generator particles
  const Int_t kMaxConstituent = 5000;             // and 1000 of the jet constituents
  const Int_t kMaxNsjBeta = 4;                    // Various tau beta values
  ULong64_t eventNum = 0;                         // need to store an event number for uproot access

  Int_t nGen=0;
  Int_t   gen_orig[kMaxGen]; // original index for debugging
  Float_t gen_pt[kMaxGen];
  Float_t gen_eta[kMaxGen];
  Float_t gen_phi[kMaxGen];
  Float_t gen_m[kMaxGen];
  Int_t   gen_id[kMaxGen];
  Int_t   gen_flags[kMaxGen];
  Int_t   gen_status[kMaxGen];
  Int_t   gen_mother1[kMaxGen];
  Int_t   gen_mother2[kMaxGen];
  Int_t   gen_daughter1[kMaxGen];
  Int_t   gen_daughter2[kMaxGen];
  Int_t   gen_col[kMaxGen];
  Float_t   gen_vxx[kMaxGen];
  Float_t   gen_vyy[kMaxGen];
  Float_t   gen_vzz[kMaxGen];
  Float_t   gen_tau[kMaxGen];


  Int_t nConstituent=0;
  Int_t   constituent_orig[kMaxConstituent]; // original index for debugging
  Float_t constituent_pt[kMaxConstituent];
  Float_t constituent_eta[kMaxConstituent];
  Float_t constituent_phi[kMaxConstituent];
  Float_t constituent_m[kMaxConstituent];
  Int_t   constituent_jetndx[kMaxConstituent];
  Int_t   constituent_subjetndx[kMaxConstituent];
  Int_t   constituent_id[kMaxConstituent];
  Int_t   constituent_flags[kMaxConstituent];
  Int_t   constituent_status[kMaxConstituent];
  Int_t   constituent_mother1[kMaxConstituent];
  Int_t   constituent_mother2[kMaxConstituent];
  Int_t   constituent_daughter1[kMaxConstituent];
  Int_t   constituent_daughter2[kMaxConstituent];
  Int_t   constituent_col[kMaxConstituent];
  Float_t   constituent_vxx[kMaxConstituent];
  Float_t   constituent_vyy[kMaxConstituent];
  Float_t   constituent_vzz[kMaxConstituent];
  Float_t   constituent_tau[kMaxConstituent];

  T->SetBranchAddress("eventNum",    &eventNum);
  T->SetBranchAddress("nGen",    &nGen);
  T->SetBranchAddress("gen_pt",        &gen_pt);
  T->SetBranchAddress("gen_eta",       &gen_eta);
  T->SetBranchAddress("gen_phi",       &gen_phi);
  T->SetBranchAddress("gen_m",         &gen_m);
  T->SetBranchAddress("gen_flags",     &gen_flags);
  T->SetBranchAddress("gen_id",        &gen_id);
  T->SetBranchAddress("gen_status",    &gen_status);
  T->SetBranchAddress("gen_mother1",   &gen_mother1);   
  T->SetBranchAddress("gen_mother2",   &gen_mother2);  
  T->SetBranchAddress("gen_daughter1", &gen_daughter1); 
  T->SetBranchAddress("gen_daughter2", &gen_daughter2);
  T->SetBranchAddress("gen_col",       &gen_col);
  T->SetBranchAddress("gen_vxx",       &gen_vxx);      
  T->SetBranchAddress("gen_vyy",       &gen_vyy);      
  T->SetBranchAddress("gen_vzz",       &gen_vzz);      
  T->SetBranchAddress("gen_tau",       &gen_tau);      


  T->SetBranchAddress("nConstituent",    &nConstituent);
  T->SetBranchAddress("constituent_pt",        &constituent_pt);
  T->SetBranchAddress("constituent_eta",       &constituent_eta);
  T->SetBranchAddress("constituent_phi",       &constituent_phi);
  T->SetBranchAddress("constituent_m",         &constituent_m);
  T->SetBranchAddress("constituent_flags",     &constituent_flags);     
  T->SetBranchAddress("constituent_id",        &constituent_id);        
  T->SetBranchAddress("constituent_jetndx",    &constituent_jetndx);    
  T->SetBranchAddress("constituent_subjetndx", &constituent_subjetndx); 
  T->SetBranchAddress("constituent_status",    &constituent_status);    
  T->SetBranchAddress("constituent_mother1",   &constituent_mother1);   
  T->SetBranchAddress("constituent_mother2",   &constituent_mother2);   
  T->SetBranchAddress("constituent_daughter1", &constituent_daughter1); 
  T->SetBranchAddress("constituent_daughter2", &constituent_daughter2); 
  T->SetBranchAddress("constituent_col",       &constituent_col);       
  T->SetBranchAddress("constituent_vxx",       &constituent_vxx);       
  T->SetBranchAddress("constituent_vyy",       &constituent_vyy);       
  T->SetBranchAddress("constituent_vzz",       &constituent_vzz);       
  T->SetBranchAddress("constituent_tau",       &constituent_tau);       

  
  
 // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < T->GetEntries(); ++iEvent) {
    T->GetEntry(iEvent);

    // Create the list of input 4-vectors
    std::vector<fastjet::PseudoJet> inputs(nConstituent);
    for ( auto i = 0; i < nConstituent; ++i ){
      fastjet::PseudoJet ijet;
      ijet.reset_PtYPhiM( constituent_pt[i], constituent_eta[i], constituent_phi[i], constituent_m[i] );
      ijet.set_user_index( i );
      inputs.emplace_back( ijet );
    }
    // Define the cluster sequence
    fastjet::ClusterSequence cs(inputs, jet_def);
    // Cluster the jets
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(ptmin));
    std::cout << "nConstituent = " << nConstituent << std::endl;
    std::cout << "jets.size() = " << jets.size() << std::endl;
  }

  file->Close();
  // DO NOT delete T. 

  // Done.
  return 0;
}
