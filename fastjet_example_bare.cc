// taken from main92.cc from pythia8 as a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to access Pythia 8 program elements.
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include <fstream>
#include <cmath>


int main(int argc, char ** argv) {

  std::vector<double> e = {1, 40, 5, 5, 20, 10, 30, 1};
  std::vector<double> phi = {-2.6, -2.4, -2.0, -0.5, 0.0, 0.5, 2.4, 2.7}; 

  std::vector<fastjet::PseudoJet> fj_particles; 
  for ( unsigned int i = 0; i < e.size(); ++i ) {
    fj_particles.push_back( fastjet::PseudoJet( e[i] * cos(phi[i]), e[i] * sin(phi[i]), 0, e[i]) );
  }

  bool verbose = true;
  std::ofstream output("fastjet_example_bare_output.txt");
  
  // Define the AK8 jet finder.
  double R = 0.8, ptmin = 0.0;
  fastjet::JetDefinition jet_def_antikt(fastjet::antikt_algorithm, R);

  // Define some groomed jets: soft drop beta=0 (mmdt), beta=1, beta=2
  // give the soft drop groomer a short name
  // Use a symmetry cut z > z_cut R^beta
  // By default, there is no mass-drop requirement
  double z_cut = 0.10;
  double beta0  = 0.0;
  fastjet::contrib::SoftDrop sdb0(beta0, z_cut);
  double beta1  = 1.0;
  fastjet::contrib::SoftDrop sdb1(beta1, z_cut);
  double beta2  = 2.0;
  fastjet::contrib::SoftDrop sdb2(beta2, z_cut);
  
  fastjet::ClusterSequence cs_antikt(fj_particles, jet_def_antikt);
  std::vector<fastjet::PseudoJet> jets_antikt = fastjet::sorted_by_pt(cs_antikt.inclusive_jets(ptmin));

  // Loop over all different jet clusterings
  for ( auto jets : {jets_antikt} ) {

    output << "Processing jet collection" << std::endl;
    if ( jets.size() > 0 ) {
      auto ibegin = jets.begin();
      auto iend = jets.end();      
      for ( auto ijet=ibegin;ijet!=iend;++ijet ) {
	auto constituents = ijet->constituents();
      

	fastjet::PseudoJet jet = *ijet;
	fastjet::PseudoJet jetsdb0 = sdb0( jet );
	fastjet::PseudoJet jetsdb1 = sdb1( jet );
	fastjet::PseudoJet jetsdb2 = sdb2( jet );
      
	char buff[1000];
	sprintf(buff, " %8.4f %6.2f %6.2f %8.4f %8.4f %8.4f %8.4f", jet.perp(), jet.eta(), jet.phi(), jet.m(), jetsdb0.m(), jetsdb1.m(), jetsdb2.m() );
	output << buff << std::endl;
	if ( verbose ) {

	  sprintf(buff, " ungroomed: %8.4f %6.2f %6.2f %8.4f", jet.perp(), jet.eta(), jet.phi(), jet.m() );
	  std::cout << buff << std::endl;
	  sprintf(buff, " sd beta=0: %8.4f %6.2f %6.2f %8.4f", jetsdb0.perp(), jetsdb0.eta(), jetsdb0.phi(), jetsdb0.m() );
	  std::cout << buff << std::endl;
	  sprintf(buff, " sd beta=1: %8.4f %6.2f %6.2f %8.4f", jetsdb1.perp(), jetsdb1.eta(), jetsdb1.phi(), jetsdb1.m() );
	  std::cout << buff << std::endl;
	  sprintf(buff, " sd beta=2: %8.4f %6.2f %6.2f %8.4f", jetsdb2.perp(), jetsdb2.eta(), jetsdb2.phi(), jetsdb2.m() );
	  std::cout << buff << std::endl;
	}
      }
    }
  }  
  // End event loop.
  
  // Done.
  return 0;
}
