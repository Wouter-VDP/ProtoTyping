////////////////////////////////////////////////////////////////////////
// Class:       CRTdata
// Plugin Type: analyzer (art v2_11_03)
// File:        CRTdata_module.cc
// Mon Nov 12 2018 by Wouter Van De Pontseele
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"

#include "ubobj/CRT/CRTHit.hh"
#include "ubobj/RawData/DAQHeaderTimeUBooNE.h"
#include <TTree.h>

class CRTdata;

class CRTdata : public art::EDAnalyzer
{
public:
  explicit CRTdata(fhicl::ParameterSet const &p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  CRTdata(CRTdata const &) = delete;
  CRTdata(CRTdata &&) = delete;
  CRTdata &operator=(CRTdata const &) = delete;
  CRTdata &operator=(CRTdata &&) = delete;

  // Required functions.
  void analyze(art::Event const &e) override;

private:
  const std::string m_CRTHitProducer = "CRTHitProducer";
  const std::string m_DAQHeaderProducer = "DAQHeaderProducer";
  const float m_DTOffset = 68600.; // in microseconds

  // TTree Declaration.
  TTree *fCRTTree;

  // Event info.
  uint fRun, fSubrun, fEvent;
  uint fNhits_in_event;

  // CRT hit info.
  uint fPlane;
  float fTime;
  float fPE;
  float fX;
  float fX_err;
  float fY;
  float fY_err;
  float fZ;
  float fZ_err;

  uint32_t ts0_s;
  int8_t ts0_s_corr;

  uint32_t ts0_ns;
  int32_t ts0_ns_corr;
  int32_t ts1_ns;
};

CRTdata::CRTdata(fhicl::ParameterSet const &p)
    : EDAnalyzer(p) // ,
                    // More initializers here.
{
  art::ServiceHandle<art::TFileService> tfs;

  //// Tree for every CRT hit
  fCRTTree = tfs->make<TTree>("CRT", "CRT Tree");

  fCRTTree->Branch("event", &fEvent, "event/i");
  fCRTTree->Branch("run", &fRun, "run/i");
  fCRTTree->Branch("subrun", &fSubrun, "subrun/i");
  fCRTTree->Branch("nhits", &fNhits_in_event, "nhits/i");

  fCRTTree->Branch("plane", &fPlane, "plane/i");
  fCRTTree->Branch("time", &fTime, "time/F");
  fCRTTree->Branch("pe", &fPE, "pe/F");
  fCRTTree->Branch("x", &fX, "x/F");
  fCRTTree->Branch("x_err", &fX_err, "x_err/F");
  fCRTTree->Branch("y", &fY, "y/F");
  fCRTTree->Branch("y_err", &fY_err, "y_err/F");
  fCRTTree->Branch("z", &fZ, "z/F");
  fCRTTree->Branch("z_err", &fZ_err, "z_err/F");
}

void CRTdata::analyze(art::Event const &e)
{
  // Fill the event info.
  fEvent = e.event();
  fSubrun = e.subRun();
  fRun = e.run();

  // Declare an object for the GPS timestamp of the event so that you can offset the cosmic t0 times.
  art::Handle<raw::DAQHeaderTimeUBooNE> rawHandle_DAQHeader;
  e.getByLabel(m_DAQHeaderProducer, rawHandle_DAQHeader);

  raw::DAQHeaderTimeUBooNE const &my_DAQHeader(*rawHandle_DAQHeader);
  art::Timestamp evtTimeGPS = my_DAQHeader.gps_time();
  float evt_timeGPS_nsec = evtTimeGPS.timeLow();

  // Load CRT hits.
  art::Handle<std::vector<crt::CRTHit>> crthit_h;
  e.getByLabel(m_CRTHitProducer, crthit_h);
  if (!crthit_h.isValid())
  {
    std::cerr << "... could not locate CRT Hits!" << std::endl;
    throw std::exception();
  }

  // Set the variable for the number of CRT hits in the event.
  fNhits_in_event = crthit_h->size();
  std::cout << "Number of CRT hits in event:" << fNhits_in_event << std::endl;

  for (size_t j = 0; j < fNhits_in_event; j++)
  {
    fPlane = crthit_h->at(j).plane;
    fPE = crthit_h->at(j).peshit;
    fX = crthit_h->at(j).x_pos;
    fY = crthit_h->at(j).y_pos;
    fZ = crthit_h->at(j).z_pos;
    fX_err = crthit_h->at(j).x_err;
    fY_err = crthit_h->at(j).y_err;
    fZ_err = crthit_h->at(j).z_err;

    // Time of the CRT Hit wrt the event timestamp
    fTime = ((crthit_h->at(j).ts0_ns - evt_timeGPS_nsec + m_DTOffset) / 1000.);

    uint32_t ts0_s;
    int8_t ts0_s_corr;

    uint32_t ts0_ns;
    int32_t ts0_ns_corr;
    int32_t ts1_ns;

    ts0_s = crthit_h->at(j).ts0_ns;
    ts0_s_corr = crthit_h->at(j).ts0_s_corr;
    ts0_ns = crthit_h->at(j).ts0_ns;
    ts1_ns = crthit_h->at(j).ts1_ns;
    ts0_ns_corr = crthit_h->at(j).ts0_ns_corr;

    std::cout << "fTime: " << fTime << "\t ts0_s: " << ts0_s << "\t ts0_ns: " << ts0_ns << "\t ts1_ns: " << ts1_ns;
    std::cout << "\t ts0_s_corr: " << ts0_s_corr << "\t ts0_ns_corr: " << ts0_ns_corr << std::endl;

    fCRTTree->Fill();
  }
}

DEFINE_ART_MODULE(CRTdata)
