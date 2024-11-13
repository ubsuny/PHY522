// taken from main92.cc from pythia8 as a part of the PYTHIA event generator.
// Copyright (C) 2018 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL v2 or later, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// Header file to access Pythia 8 program elements.
#include "fastjet/ClusterSequence.hh"
#include "fastjet/contrib/SoftDrop.hh"
#include "fastjet/contrib/Nsubjettiness.hh"
#include <fstream>
#include <cmath>

int main(int argc, char ** argv) {

  
  // Define the AK8 jet finder.
  double R = 1.0, ptmin = 0.0;
  fastjet::JetDefinition jet_def_antikt(fastjet::antikt_algorithm, R);
  fastjet::JetDefinition jet_def_kt(fastjet::kt_algorithm, R);
  fastjet::JetDefinition jet_def_ca(fastjet::cambridge_algorithm, R);

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

  // Now make the n-subjettiness ratios
  double beta=1; 
  fastjet::contrib::Nsubjettiness calc_tau1(1, fastjet::contrib::WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
  fastjet::contrib::Nsubjettiness calc_tau2(2, fastjet::contrib::WTA_KT_Axes(), fastjet::contrib::UnnormalizedMeasure(beta));
  
  bool verbose = true;
  std::ofstream output("fastjet_example_bare_output.txt");
  
  std::vector<double> e = {1, 40, 5, 5, 20, 10, 30, 1};
  std::vector<double> phi = {-2.6, -2.4, -2.0, -0.5, 0.0, 0.5, 2.4, 2.7}; 

  std::vector<fastjet::PseudoJet> fj_particles;
  std::cout << "Inputs: " << std::endl;
  output << "Inputs: " << std::endl;
  for ( unsigned int i = 0; i < e.size(); ++i ) {
    fastjet::PseudoJet jet(  fastjet::PseudoJet( e[i] * cos(phi[i]), e[i] * sin(phi[i]), 0, e[i]) );
    jet.set_user_index(i);
    fj_particles.push_back( jet );
    char buff[1000];
    sprintf(buff, " i=%3d pt=%8.4f eta=%6.2f phi=%6.2f m=%8.4f", i, jet.perp(), jet.eta(), jet.phi(), jet.m() );
    output << buff << std::endl;
    std::cout << buff << std::endl;
  }

  std::cout << "Anti-kt: " << std::endl;
  fastjet::ClusterSequence cs_antikt(fj_particles, jet_def_antikt, true);
  std::vector<fastjet::PseudoJet> jets_antikt = fastjet::sorted_by_pt(cs_antikt.inclusive_jets(ptmin));

  std::cout << "kt: " << std::endl;
  fastjet::ClusterSequence cs_kt(fj_particles, jet_def_kt, true);
  std::vector<fastjet::PseudoJet> jets_kt = fastjet::sorted_by_pt(cs_kt.inclusive_jets(ptmin));

  std::cout << "C/A: " << std::endl;
  fastjet::ClusterSequence cs_ca(fj_particles, jet_def_ca, true);
  std::vector<fastjet::PseudoJet> jets_ca = fastjet::sorted_by_pt(cs_ca.inclusive_jets(ptmin));

  // Loop over all different jet clusterings
  for ( auto colls :
	  {std::make_tuple("antikt", jets_antikt),
	   std::make_tuple("kt", jets_kt),
	   std::make_tuple("ca", jets_ca)
	  }
	) {

    auto name = std::get<0>(colls);
    auto jets = std::get<1>(colls); 

    output << name << std::endl;
    std::cout << name << std::endl;
    if ( jets.size() > 0 ) {
      auto ibegin = jets.begin();
      auto iend = jets.end();      
      for ( auto ijet=ibegin;ijet!=iend;++ijet ) {
	auto constituents = ijet->constituents();
      

	fastjet::PseudoJet jet = *ijet;
	fastjet::PseudoJet jetsdb0 = sdb0( jet );
	double tau1 = calc_tau1( jet );
	double tau2 = calc_tau2( jet );
	double tau2_sdb0 = calc_tau2( jetsdb0 );

	std::vector<fastjet::PseudoJet> subjets = jetsdb0.pieces();
	std::vector<fastjet::PseudoJet> groomed_constituents;
	for ( auto ig : subjets ){
	  auto consts = ig.constituents();
	  groomed_constituents.insert( groomed_constituents.end(), consts.begin(), consts.end() );
	}
      
	char buff[1000];
	sprintf(buff, " i=%2d pt=%4.2f eta=%4.2f phi=%4.2f m=%5.2f m(sd0)=%5.2f tau1=%4.2f tau2=%4.2f tau2_sdb0=%4.2f, constituents: ", ijet-ibegin, jet.perp(), jet.eta(), jet.phi(), jet.m(), jetsdb0.m(), tau1, tau2, tau2_sdb0 );
	output << buff;
	std::cout << buff;
	for ( auto ic : constituents ){
	  sprintf(buff, " %2d", ic.user_index() );
	  output << buff;
	  std::cout << buff;
	}
	output << ", groomed : ";
	std::cout << ", groomed : ";
	for ( auto ic : groomed_constituents ){
	  sprintf(buff, " %2d", ic.user_index() );
	  output << buff;
	  std::cout << buff;
	}	
	output << std::endl;
	std::cout << std::endl;
      }
    }
  }  
  // End event loop.
  
  // Done.
  return 0;
}
