// main92.cc is a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program.
// Modified by Rene Brun and Axel Naumann to put the Pythia::event
// into a TTree.

// Header file to access Pythia 8 program elements.
#include "Pythia8/Pythia.h"
#include "Pythia8/Basics.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh" // In external code, this should be fastjet/contrib/Nsubjettiness.hh
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"


// ROOT, for saving Pythia events as trees in a file.
#include "TTree.h"
#include "TFile.h"

using namespace Pythia8;

class CompareIndex {
public:
  CompareIndex( fastjet::PseudoJet j ) : i_(j.user_index()){
  }
  CompareIndex( unsigned int j = 0 ) : i_(j){
  }

  bool operator() ( fastjet::PseudoJet const & j ) const { return j.user_index() == i_; }
  bool operator() ( unsigned int j ) const { return j == i_; }

protected : 
  unsigned int i_; 
};

int main(int argc, char ** argv) {

  if ( argc < 4 ) {
    std::cout << "usage: " << argv[0] << " config_file root_file n_events <optional: seed (-1 = default, 0=use time, or input your own)> <optional: ptcut>" << std::endl;
    return 0;
  }


  // Define the AK8 jet finder.
  double R = 0.8, ptmin = 50.0, lepfrac = 0.9;
  fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);

  bool verbose = false;

  char * configfile = argv[1];
  char * outfile = argv[2];
  unsigned int nEvents = atol(argv[3]);
  long seed = -1; 
  if ( argc > 5 ) {
    seed = atol(argv[4]);
  }
  if ( argc > 6 ) {
    ptmin = atof( argv[5]);
  }


  double z_cut = 0.10;
  double beta  = 0.0;
  fastjet::contrib::SoftDrop sd(beta, z_cut);

  // Create Pythia instance. Read config from a text file. 
  Pythia pythia;
  char buff[1000];
  sprintf(buff, "Random:seed = %d", seed);
  pythia.readString("Random:setSeed = on");
  pythia.readString(buff);
  std::ifstream config( configfile );
  while (!config.eof() ) {
    std::string line;
    std::getline( config, line );
    if ( line[0] != '!' && line != "" && line != "\n" ){
      pythia.readString(line);
    }
  }
  pythia.init();

  // Set up the ROOT TFile and TTree.
  TFile *file = TFile::Open(outfile,"recreate");
  Event *event = &pythia.event;
  const Int_t kMaxJet = 10;                       // Stores leading 10 jets
  const Int_t kMaxGen = 5000;                     // and 1000 of the generator particles
  const Int_t kMaxConstituent = 5000;             // and 1000 of the jet constituents
  const Int_t kMaxNsjBeta = 4;                    // Various tau beta values
  ULong64_t eventNum = 0;                         // need to store an event number for uproot access
  Int_t nJet=0;
  Float_t jet_pt[kMaxJet];
  Float_t jet_eta[kMaxJet];
  Float_t jet_phi[kMaxJet];
  Float_t jet_m[kMaxJet];
  Float_t jet_msd[kMaxJet];
  Float_t jet_tau1   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau2   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau3   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau4   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau5   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau6   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau7   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau8   [kMaxJet][kMaxNsjBeta];
  Float_t jet_tau1_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau2_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau3_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau4_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau5_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau6_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau7_sd[kMaxJet][kMaxNsjBeta];
  Float_t jet_tau8_sd[kMaxJet][kMaxNsjBeta];


  
  Int_t   jet_nc[kMaxJet];
  Int_t   jet_ic[kMaxJet][50];
  Int_t   jet_nsubjet[kMaxJet];
  Float_t jet_subjet0_pt[kMaxJet];
  Float_t jet_subjet0_eta[kMaxJet];
  Float_t jet_subjet0_phi[kMaxJet];
  Float_t jet_subjet0_m[kMaxJet];
  Float_t jet_subjet1_pt[kMaxJet];
  Float_t jet_subjet1_eta[kMaxJet];
  Float_t jet_subjet1_phi[kMaxJet];
  Float_t jet_subjet1_m[kMaxJet];
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

  TTree * T = new TTree("T","ev1 Tree");         // Allocate the tree, but DO NOT DELETE IT since ROOT takes ownership magically.
  T->Branch("eventNum",    &eventNum,  "eventNum/l");
  T->Branch("nJet",    &nJet,  "nJet/I");
  T->Branch("jet_pt",  &jet_pt,  "jet_pt[nJet]/F");
  T->Branch("jet_eta", &jet_eta, "jet_eta[nJet]/F");
  T->Branch("jet_phi", &jet_phi, "jet_phi[nJet]/F");
  T->Branch("jet_m",   &jet_m,   "jet_m[nJet]/F");
  T->Branch("jet_msd",   &jet_msd,   "jet_msd[nJet]/F");
  T->Branch("jet_tau1",   &jet_tau1,   "jet_tau1[nJet][4]/F");
  T->Branch("jet_tau2",   &jet_tau2,   "jet_tau2[nJet][4]/F");
  T->Branch("jet_tau3",   &jet_tau3,   "jet_tau3[nJet][4]/F");
  T->Branch("jet_tau4",   &jet_tau4,   "jet_tau4[nJet][4]/F");
  T->Branch("jet_tau5",   &jet_tau5,   "jet_tau5[nJet][4]/F");
  T->Branch("jet_tau6",   &jet_tau6,   "jet_tau6[nJet][4]/F");
  T->Branch("jet_tau7",   &jet_tau7,   "jet_tau7[nJet][4]/F");
  T->Branch("jet_tau8",   &jet_tau8,   "jet_tau8[nJet][4]/F");
  T->Branch("jet_tau1_sd",   &jet_tau1_sd,   "jet_tau1_sd[nJet][4]/F");
  T->Branch("jet_tau2_sd",   &jet_tau2_sd,   "jet_tau2_sd[nJet][4]/F");
  T->Branch("jet_tau3_sd",   &jet_tau3_sd,   "jet_tau3_sd[nJet][4]/F");
  T->Branch("jet_tau4_sd",   &jet_tau4_sd,   "jet_tau4_sd[nJet][4]/F");
  T->Branch("jet_tau5_sd",   &jet_tau5_sd,   "jet_tau5_sd[nJet][4]/F");
  T->Branch("jet_tau6_sd",   &jet_tau6_sd,   "jet_tau6_sd[nJet][4]/F");
  T->Branch("jet_tau7_sd",   &jet_tau7_sd,   "jet_tau7_sd[nJet][4]/F");
  T->Branch("jet_tau8_sd",   &jet_tau8_sd,   "jet_tau8_sd[nJet][4]/F");
  T->Branch("jet_nc",  &jet_nc,  "jet_nc[nJet]/I");
  T->Branch("jet_ic",  &jet_ic,  "jet_ic[nJet][50]/I");
  T->Branch("jet_nsubjet",  &jet_nsubjet,  "jet_nsubjet[nJet]/I");
  T->Branch("jet_subjet0_pt",  &jet_subjet0_pt,  "jet_subjet0_pt[nJet]/F");
  T->Branch("jet_subjet0_eta", &jet_subjet0_eta, "jet_subjet0_eta[nJet]/F");
  T->Branch("jet_subjet0_phi", &jet_subjet0_phi, "jet_subjet0_phi[nJet]/F");
  T->Branch("jet_subjet0_m",   &jet_subjet0_m,   "jet_subjet0_m[nJet]/F");
  T->Branch("jet_subjet1_pt",  &jet_subjet1_pt,  "jet_subjet1_pt[nJet]/F");
  T->Branch("jet_subjet1_eta", &jet_subjet1_eta, "jet_subjet1_eta[nJet]/F");
  T->Branch("jet_subjet1_phi", &jet_subjet1_phi, "jet_subjet1_phi[nJet]/F");
  T->Branch("jet_subjet1_m",   &jet_subjet1_m,   "jet_subjet1_m[nJet]/F");    
  T->Branch("nGen",    &nGen,  "nGen/I");
  T->Branch("gen_pt",        &gen_pt,  "gen_pt[nGen]/F");
  T->Branch("gen_eta",       &gen_eta, "gen_eta[nGen]/F");
  T->Branch("gen_phi",       &gen_phi, "gen_phi[nGen]/F");
  T->Branch("gen_m",         &gen_m,   "gen_m[nGen]/F");
  T->Branch("gen_flags",     &gen_flags,     "gen_flags[nGen]/I");
  T->Branch("gen_id",        &gen_id,        "gen_id[nGen]/I"       );
  T->Branch("gen_status",    &gen_status,    "gen_status[nGen]/I"   );
  T->Branch("gen_mother1",   &gen_mother1,   "gen_mother1[nGen]/I"  );
  T->Branch("gen_mother2",   &gen_mother2,   "gen_mother2[nGen]/I"  );
  T->Branch("gen_daughter1", &gen_daughter1, "gen_daughter1[nGen]/I");
  T->Branch("gen_daughter2", &gen_daughter2, "gen_daughter2[nGen]/I");
  T->Branch("gen_col",       &gen_col,       "gen_col[nGen]/I"      );
  T->Branch("gen_vxx",       &gen_vxx,       "gen_vxx[nGen]/F"      );
  T->Branch("gen_vyy",       &gen_vyy,       "gen_vyy[nGen]/F"      );
  T->Branch("gen_vzz",       &gen_vzz,       "gen_vzz[nGen]/F"      );
  T->Branch("gen_tau",       &gen_tau,       "gen_tau[nGen]/F"      );


  T->Branch("nConstituent",    &nConstituent,  "nConstituent/I");
  T->Branch("constituent_pt",        &constituent_pt,  "constituent_pt[nConstituent]/F");
  T->Branch("constituent_eta",       &constituent_eta, "constituent_eta[nConstituent]/F");
  T->Branch("constituent_phi",       &constituent_phi, "constituent_phi[nConstituent]/F");
  T->Branch("constituent_m",         &constituent_m,   "constituent_m[nConstituent]/F");
  T->Branch("constituent_flags",     &constituent_flags,     "constituent_flags[nConstituent]/I");
  T->Branch("constituent_id",        &constituent_id,        "constituent_id[nConstituent]/I"       );
  T->Branch("constituent_jetndx",    &constituent_jetndx,    "constituent_jetndx[nConstituent]/I"   );
  T->Branch("constituent_subjetndx", &constituent_subjetndx, "constituent_subjetndx[nConstituent]/I");
  T->Branch("constituent_status",    &constituent_status,    "constituent_status[nConstituent]/I"   );
  T->Branch("constituent_mother1",   &constituent_mother1,   "constituent_mother1[nConstituent]/I"  );
  T->Branch("constituent_mother2",   &constituent_mother2,   "constituent_mother2[nConstituent]/I"  );
  T->Branch("constituent_daughter1", &constituent_daughter1, "constituent_daughter1[nConstituent]/I");
  T->Branch("constituent_daughter2", &constituent_daughter2, "constituent_daughter2[nConstituent]/I");
  T->Branch("constituent_col",       &constituent_col,       "constituent_col[nConstituent]/I"      );
  T->Branch("constituent_vxx",       &constituent_vxx,       "constituent_vxx[nConstituent]/F"      );
  T->Branch("constituent_vyy",       &constituent_vyy,       "constituent_vyy[nConstituent]/F"      );
  T->Branch("constituent_vzz",       &constituent_vzz,       "constituent_vzz[nConstituent]/F"      );
  T->Branch("constituent_tau",       &constituent_tau,       "constituent_tau[nConstituent]/F"      );

  
  
 // Begin event loop. Generate event; skip if generation aborted.
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    eventNum = iEvent; 
    nConstituent = nGen = nJet = 0;
    if (!pythia.next()) continue;
    if ( verbose ) 
      std::cout << "Generating event " << iEvent << std::endl;
    for ( auto x : jet_pt ) x=0.0;
    for ( auto x : jet_eta ) x=0.0;
    for ( auto x : jet_phi ) x=0.0;
    for ( auto x : jet_m ) x=0.0;
    for ( auto x : jet_msd ) x=0.0;


    for ( auto x : jet_tau1    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau2    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau3    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau4    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau5    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau6    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau7    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau8    ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau1_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau2_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau3_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau4_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau5_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau6_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau7_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_tau8_sd ) for ( unsigned int j = 0; j < kMaxNsjBeta; ++j )  x[j]=0.0;
    for ( auto x : jet_nc ) x=0;
    for ( auto x : jet_nsubjet ) x=0;
    for ( auto i = 0; i < kMaxJet; ++i )
      for ( auto j = 0; j < 50; ++j )
	jet_ic[i][j] = 0;
    for ( auto x : jet_subjet0_pt ) x=0.0;
    for ( auto x : jet_subjet0_eta ) x=0.0;
    for ( auto x : jet_subjet0_phi ) x=0.0;
    for ( auto x : jet_subjet0_m ) x=0.0;
    for ( auto x : jet_subjet1_pt ) x=0.0;
    for ( auto x : jet_subjet1_eta ) x=0.0;
    for ( auto x : jet_subjet1_phi ) x=0.0;
    for ( auto x : jet_subjet1_m ) x=0.0;    
    for ( auto x : gen_pt ) x=0.0;
    for ( auto x : gen_eta ) x=0.0;
    for ( auto x : gen_phi ) x=0.0;
    for ( auto x : gen_m ) x=0.0;
    for ( auto x : gen_flags ) x = 0;
    for ( auto x : gen_id ) x = 0;
    for ( auto x : gen_status ) x = 0;
    for ( auto x : gen_mother1 ) x = 0;
    for ( auto x : gen_mother2 ) x = 0;
    for ( auto x : gen_daughter1 ) x = 0;
    for ( auto x : gen_daughter2 ) x = 0;
    for ( auto x : gen_col ) x = 0;
    for ( auto x : gen_vxx ) x = 0.;
    for ( auto x : gen_vyy ) x = 0.;
    for ( auto x : gen_vzz ) x = 0.;
    for ( auto x : gen_tau ) x = 0.;
    for ( auto x : constituent_pt ) x=0.0;
    for ( auto x : constituent_eta ) x=0.0;
    for ( auto x : constituent_phi ) x=0.0;
    for ( auto x : constituent_m ) x=0.0;
    for ( auto x : constituent_flags ) x = 0;
    for ( auto x : constituent_jetndx ) x = -1;
    for ( auto x : constituent_subjetndx ) x = -1;
    for ( auto x : constituent_id ) x = 0;
    for ( auto x : constituent_status ) x = 0;
    for ( auto x : constituent_mother1 ) x = 0;
    for ( auto x : constituent_mother2 ) x = 0;
    for ( auto x : constituent_daughter1 ) x = 0;
    for ( auto x : constituent_daughter2 ) x = 0;
    for ( auto x : constituent_col ) x = 0;
    for ( auto x : constituent_vxx ) x = 0.;
    for ( auto x : constituent_vyy ) x = 0.;
    for ( auto x : constituent_vzz ) x = 0.;
    for ( auto x : constituent_tau ) x = 0.;

    // Dump the PYTHIA8 content.     
    // Create AK8 jets with pt > 170 GeV
    std::vector<fastjet::PseudoJet> fj_particles;
    std::map<int,int> constituentmap;
    for (int i = 0; i < event->size(); ++i){
      auto const & p = pythia.event[i];
      if ( p.isFinalPartonLevel() || p.isResonance()) {
	if ( verbose ) {
	  char buff[1000];
	  sprintf( buff, "  ndx=%6d, id=%6d, status=%6d, p4=(%6.4f,%6.2f,%6.2f,%6.4f)", i, p.id(), p.status(), p.pT(), p.eta(), p.phi(), p.m() );
	  std::cout << buff << std::endl; 
	}
	gen_pt[nGen] = p.pT();
	gen_eta[nGen] = p.eta();
	gen_phi[nGen] = p.phi();
	gen_m[nGen] = p.m();
	gen_orig[nGen] = i;
	gen_id[nGen] =         p.id();
	gen_flags[nGen] =      p.isHadron() << 3 | p.isFinal() << 2 | p.isFinalPartonLevel() << 1 | p.isVisible() << 0; 
	gen_status[nGen] =     p.status();    
	gen_mother1[nGen] =    p.mother1();   
	gen_mother2[nGen] =    p.mother2();   
	gen_daughter1[nGen] =  p.daughter1(); 
	gen_daughter2[nGen] =  p.daughter2(); 
	gen_col[nGen] =        p.col();       
	gen_vxx[nGen] =        p.xProd();       
	gen_vyy[nGen] =        p.yProd();       
	gen_vzz[nGen] =        p.zProd();       
	gen_tau[nGen] =        p.tau();             
	++nGen;
	if ( nGen >= kMaxGen ){
	  std::cout << "too many particles in event " << iEvent << ", storing first " << kMaxGen << std::endl;
	  break;
	}
      } else if ( p.isFinal() ) {

	constituent_pt[nConstituent] = p.pT();
	constituent_eta[nConstituent] = p.eta();
	constituent_phi[nConstituent] = p.phi();
	constituent_m[nConstituent] = p.m();
	constituent_orig[nConstituent] = i;
	constituent_jetndx[nConstituent] =     -1; // set later
	constituent_subjetndx[nConstituent] =  -1; // set later
	constituent_id[nConstituent] =         p.id();
	constituent_flags[nConstituent] =      p.isHadron() << 3 | p.isFinal() << 2 | p.isFinalPartonLevel() << 1 | p.isVisible() << 0; 
	constituent_status[nConstituent] =     p.status();    
	constituent_mother1[nConstituent] =    p.mother1();   
	constituent_mother2[nConstituent] =    p.mother2();   
	constituent_daughter1[nConstituent] =  p.daughter1(); 
	constituent_daughter2[nConstituent] =  p.daughter2(); 
	constituent_col[nConstituent] =        p.col();       
	constituent_vxx[nConstituent] =        p.xProd();       
	constituent_vyy[nConstituent] =        p.yProd();       
	constituent_vzz[nConstituent] =        p.zProd();       
	constituent_tau[nConstituent] =        p.tau();             
	if ( p.isFinal() ) {
	  fj_particles.emplace_back( p.px(), p.py(), p.pz(), p.e()  );
	  fj_particles.back().set_user_index( i );
	  constituentmap[i] = nConstituent;
	} 
	++nConstituent;
	if ( nConstituent >= kMaxConstituent ){
	  std::cout << "too many jet constituents in event " << iEvent << ", storing first " << kMaxConstituent << std::endl;
	  break;
	}
      }
    }
    if ( verbose) std::cout << "About to cluster" << std::endl;
    fastjet::ClusterSequence cs(fj_particles, jet_def);
    std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(ptmin));

    if ( verbose ) std::cout << "About to loop over jets" << std::endl;
    nJet = 0;
    if ( jets.size() > 0 ) {
      auto ibegin = jets.begin();
      auto iend = jets.end();
      for ( auto ijet=ibegin;ijet!=iend;++ijet ) {
	if ( verbose ) std::cout << "processing jet " << ijet - ibegin << std::endl;

	auto sd_jet =  sd(*ijet);
	auto constituents = ijet->constituents();

	// Get the fraction of the jet originating from leptons.
	// This is to remove jets that are comprised entirely of isolated leptons
	// such as Z->ll.
	// We ignore jets with >90% of their energy from leptons. 
	unsigned nlepton = 0;
	auto lepp4 = fastjet::PseudoJet();
	for ( auto icon = constituents.begin(); icon != constituents.end(); ++icon ) {
	  auto const & py8part = pythia.event[ icon->user_index() ];
	  if ( std::abs( py8part.id() ) > 10 && std::abs(py8part.id()) < 16){
	    if ( verbose ){
	      char buff[1000];
	      sprintf( buff, "  lepton  :  id=%6d  p4=(%6.4f,%6.2f,%6.2f,%6.4f)", py8part.id(), icon->pt(), icon->eta(), icon->phi(), icon->m() );
	      std::cout << buff << std::endl; 	      
	    }
	    lepp4 += *icon;
	  }
	}
	if ( lepp4.e() / ijet->e() > lepfrac){
	  if ( verbose ){
	    char buff[1000];
	    sprintf( buff, "  skip jet:  ndx=%6d, nc=%6d  p4=(%6.4f,%6.2f,%6.2f,%6.4f)", ijet-ibegin, constituents.size(), ijet->pt(), ijet->eta(), ijet->phi(), ijet->m() );
	    std::cout << buff << std::endl; 
	  }
	  continue;
	}

	if ( verbose ) {
	  std::cout << "getting constituents for jet :" << std::endl;
	  char buff[1000];
	  sprintf( buff, "  add  jet:  ndx=%6d, nc=%6d  p4=(%6.4f,%6.2f,%6.2f,%6.4f)", ijet-ibegin, constituents.size(), ijet->pt(), ijet->eta(), ijet->phi(), ijet->m() );
	  std::cout << buff << std::endl;
	}
	if ( nJet < kMaxJet ) { 
	  
	  jet_pt[nJet]=ijet->perp();
	  jet_eta[nJet]=ijet->eta();
	  jet_phi[nJet]=ijet->phi();
	  jet_m[nJet]=ijet->m();	  
	  jet_msd[nJet] = sd_jet.m();

	  if ( nJet < 20 ) { //N-jettiness is hard-coded to only allow up to 20 jets


	    // Define Nsubjettiness functions for beta = 1.0 using one-pass WTA KT axes

	    
	    for ( unsigned int nsj_index = 0; nsj_index < kMaxNsjBeta; ++nsj_index) {
	      double beta_nsj = 0.5 + 0.5*nsj_index;
	      const Int_t max_nsj = 8;

	      fastjet::contrib::Nsubjettiness nSub1_beta1(1, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub2_beta1(2, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub3_beta1(3, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub4_beta1(4, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub5_beta1(5, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub6_beta1(6, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub7_beta1(7, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));
	      fastjet::contrib::Nsubjettiness nSub8_beta1(8, fastjet::contrib::OnePass_WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta_nsj));

	      jet_tau1[nJet][nsj_index] = nSub1_beta1(*ijet);
	      jet_tau2[nJet][nsj_index] = nSub2_beta1(*ijet);
	      jet_tau3[nJet][nsj_index] = nSub3_beta1(*ijet);
	      jet_tau4[nJet][nsj_index] = nSub4_beta1(*ijet);
	      jet_tau5[nJet][nsj_index] = nSub5_beta1(*ijet);
	      jet_tau6[nJet][nsj_index] = nSub6_beta1(*ijet);
	      jet_tau7[nJet][nsj_index] = nSub7_beta1(*ijet);
	      jet_tau8[nJet][nsj_index] = nSub8_beta1(*ijet);

	      jet_tau1_sd[nJet][nsj_index] = nSub1_beta1(sd_jet);
	      jet_tau2_sd[nJet][nsj_index] = nSub2_beta1(sd_jet);
	      jet_tau3_sd[nJet][nsj_index] = nSub3_beta1(sd_jet);
	      jet_tau4_sd[nJet][nsj_index] = nSub4_beta1(sd_jet);
	      jet_tau5_sd[nJet][nsj_index] = nSub5_beta1(sd_jet);
	      jet_tau6_sd[nJet][nsj_index] = nSub6_beta1(sd_jet);
	      jet_tau7_sd[nJet][nsj_index] = nSub7_beta1(sd_jet);
	      jet_tau8_sd[nJet][nsj_index] = nSub8_beta1(sd_jet);
	    }

	  }
	  
	  jet_nc[nJet] = constituents.size();
	  auto subjets = sd_jet.pieces();
	  
	  jet_nsubjet[nJet] = subjets.size(); 
	  std::vector<fastjet::PseudoJet> sj0_pieces, sj1_pieces; 

	  if ( subjets.size() >= 1 ) {
	    jet_subjet0_pt[nJet]  = subjets[0].perp();
	    jet_subjet0_eta[nJet] = subjets[0].eta();
	    jet_subjet0_phi[nJet] = subjets[0].phi();
	    jet_subjet0_m[nJet]   = subjets[0].m();	    
	    auto ipieces = subjets[0].constituents();
	    sj0_pieces.insert( sj0_pieces.begin(), ipieces.begin(), ipieces.end() );
	  }
	  if ( subjets.size() >= 2 ) {
	    jet_subjet1_pt[nJet]  = subjets[1].perp();
	    jet_subjet1_eta[nJet] = subjets[1].eta();
	    jet_subjet1_phi[nJet] = subjets[1].phi();
	    jet_subjet1_m[nJet]   = subjets[1].m();
	    auto ipieces = subjets[1].constituents();
	    sj1_pieces.insert( sj1_pieces.begin(), ipieces.begin(), ipieces.end() );
	  } else{
	    jet_subjet1_pt[nJet]  = 0;
	    jet_subjet1_eta[nJet] = 0;
	    jet_subjet1_phi[nJet] = 0;
	    jet_subjet1_m[nJet]   = 0;
	  }
	  if ( constituents.size() > 0 ) {	    
	    auto jbegin = constituents.begin();
	    auto jend = constituents.end();
	    Int_t nParticle = 0;
	    for ( auto iparticle=jbegin; iparticle != jend;++iparticle ){
	      
	      auto index = iparticle->user_index();
	      jet_ic[nJet][nParticle] = constituentmap[index];
	      constituent_jetndx[constituentmap[index]] = nJet;
	      
	      auto sj0_find = std::find_if( sj0_pieces.begin(), sj0_pieces.end(), CompareIndex(*iparticle) );
	      auto sj1_find = std::find_if( sj1_pieces.begin(), sj1_pieces.end(), CompareIndex(*iparticle) );
	      if ( sj0_find != sj0_pieces.end() ){
		constituent_subjetndx[constituentmap[index]]=0;
	      } else if ( sj1_find != sj1_pieces.end() ) {
		constituent_subjetndx[constituentmap[index]]=1;
	      } else {
		constituent_subjetndx[constituentmap[index]]=-1;
	      }
	      if ( verbose ) std::cout << index << " ";
	      ++nParticle;
	    }
	    if ( verbose) std::cout << endl;
	  }
	  ++nJet;
	}
      }
      if ( verbose ) 
	std::cout << "About to write" << std::endl;
      if ( nJet > 0 && jet_pt[0] > ptmin ) {
	// Fill the pythia event into the TTree.
	T->Fill();
      }
      if ( verbose ) 
	std::cout << "Done writing." << std::endl;
    } // end check if jets.size() > 0
  // End event loop.
  }

  // Statistics on event generation.
  pythia.stat();

  //  Write tree.
  T->Write();
  file->Close();
  // DO NOT delete T. 

  // Done.
  return 0;
}
