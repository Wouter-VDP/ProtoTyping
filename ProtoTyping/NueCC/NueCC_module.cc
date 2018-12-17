////////////////////////////////////////////////////////////////////////
// Class:       NueCC
// Plugin Type: analyzer (art v2_11_03)
// File:        NueCC_module.cc
//
// Generated at Sun Oct 28 14:38:04 2018 by Wouter Van de pontseele using cetskelgen
// from cetlib version v3_03_01.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "TTree.h"

class NueCC;


class NueCC : public art::EDAnalyzer {
public:
  explicit NueCC(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  NueCC(NueCC const &) = delete;
  NueCC(NueCC &&) = delete;
  NueCC & operator = (NueCC const &) = delete;
  NueCC & operator = (NueCC &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:
  std::string m_pfp_producer;
};

DEFINE_ART_MODULE(NueCC)

NueCC::NueCC(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{}

void NueCC::analyze(art::Event const & e)
{
  std::cout << "[NueCC::analyze]: Hello World" << std::endl;
  // Implementation of required member function here.
}
