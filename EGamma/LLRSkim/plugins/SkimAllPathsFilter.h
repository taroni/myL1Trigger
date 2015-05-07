#ifndef SkimAllPathsFilter_h
#define SkimAllPathsFilter_h

#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"

#include <iostream>
#include <memory>

//#include "TTree.h"
//#include "TLorentzVector.h"
//#include "TClonesArray.h"
//#include "TParticle.h"

//
// class declaration
//

class SkimAllPathsFilter : public edm::EDFilter
{
 public:
  
  explicit SkimAllPathsFilter(const edm::ParameterSet&);
  ~SkimAllPathsFilter();
  
  virtual bool filter(edm::Event&, const edm::EventSetup& );
      
 private:

  edm::InputTag electronCollection_;
  edm::InputTag muonCollection_;

  int n_init;
  int n_lepton1_skimmed;
  int n_lepton2_skimmed;
  int n_ele2TP_skimmed;
  int n_ele2TP_skimmed_nadir;
  int n_lepton3_skimmed;
  int n_ml_skimmed;
  int n_ep_skimmed;
  int n_old_skimmed;

  // 
  std::string _mode;  // Skimming mode
  std::string _eleID; // Choose type of eleID

};

#endif
