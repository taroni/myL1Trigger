#include <iostream>
#include <sstream>
#include <istream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <functional>


#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"  
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"  
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalCleaningAlgo.h"  

#include "DataFormats/EcalDigi/interface/EcalTrigPrimCompactColl.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "DataFormats/EcalDigi/interface/EcalTriggerPrimitiveSample.h"
 
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetupFwd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutSetup.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"

#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h"

#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"


#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTriggerRecord.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GtTechnicalTrigger.h"

#include "Geometry/Records/interface/IdealGeometryRecord.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFLayer.h"



#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"


#include "JetMETCorrections/MCJet/plugins/PFDataTreeProducer.h" 
#include "JetMETCorrections/MCJet/plugins/JetUtilMC.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"




#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"




#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"


#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/EcalElectronicsId.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "RecoCaloTools/Navigation/interface/CaloNavigator.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"

#include "Math/VectorUtil.h"
#include "TVector3.h"

using namespace edm;
using namespace reco;
using namespace std;


PFDataTreeProducer::PFDataTreeProducer(const edm::ParameterSet& cfg)
  : bunchstartbx_(cfg.getParameter<std::vector<int> >("bunchstartbx"))
{
  jets_          = cfg.getParameter<std::string> ("jets");
  histogramFile_ = cfg.getParameter<std::string> ("histogramFile");
  tracks_        = cfg.getParameter<std::string> ("tracks");
  vertex_coll_   = cfg.getParameter<std::string> ("vertex");
  //  jetcorr_       = cfg.getParameter<std::string> ("JetCorrectionService");
  ebhitcoll_     = cfg.getParameter<std::string> ("EBRecHitCollection");
  eehitcoll_     = cfg.getParameter<std::string> ("EERecHitCollection");
  //  pfhitcoll1_    = cfg.getParameter<std::string> ("PFRecHitCollection1");
  //  pfhitcoll2_    = cfg.getParameter<std::string> ("PFRecHitCollection2");
  isMC_          = cfg.getParameter<bool>("IsMC");
  onlineTPTag_   = cfg.getParameter<edm::InputTag>("OnlineTPs"); //modif-alex

  edm::ParameterSet cleaningPs = 
    cfg.getParameter<edm::ParameterSet>("cleaningConfig");
  cleaningAlgo_ = new EcalCleaningAlgo(cleaningPs);


}
//////////////////////////////////////////////////////////////////////////////////////////
void PFDataTreeProducer::beginJob() 
{


  file_tt = new TFile("/afs/cern.ch/cms/CAF/CMSCOMM/COMM_ECAL/azabi/CMSSW_5315_TPG/src/JetMETCorrections/MCJet/plugins/ietaphiTT.root","READ");

  ietamapTT=(TH2F*)file_tt->Get("ietaTT");
  iphimapTT=(TH2F*)file_tt->Get("iphiTT");


  file_          = new TFile(histogramFile_.c_str(),"RECREATE");

  histo_event_   = new TH1D("histo_event","",350,0.,3500.);

  eb_rechitenergy_   = new TH1D("eb_rechitenergy","",1100,-2.,20.);
  ee_rechitenergy_   = new TH1D("ee_rechitenergy","",1100,-2.,20.);

 ee_rechitenergy_notypeb_   = new TH1D("ee_rechitenergy_notypeb","",1100,-2.,20.);


 Emin_=1000.0;
 side_=5;


  eb_rechitenergy_02    = new TH1D("eb_rechitenergy_02","",1100,-2.,20.);
  eb_rechitenergy_04    = new TH1D("eb_rechitenergy_04","",1100,-2.,20.);
  eb_rechitenergy_06    = new TH1D("eb_rechitenergy_06","",1100,-2.,20.);
  eb_rechitenergy_08    = new TH1D("eb_rechitenergy_08","",1100,-2.,20.);
  eb_rechitenergy_10    = new TH1D("eb_rechitenergy_10","",1100,-2.,20.);
  eb_rechitenergy_12    = new TH1D("eb_rechitenergy_12","",1100,-2.,20.);
  eb_rechitenergy_14    = new TH1D("eb_rechitenergy_14","",1100,-2.,20.);
  eb_rechitenergy_148   = new TH1D("eb_rechitenergy_148","",1100,-2.,20.);
 
  ee_rechitenergy_16    = new TH1D("ee_rechitenergy_16","",1100,-2.,20.);
  ee_rechitenergy_18    = new TH1D("ee_rechitenergy_18","",1100,-2.,20.);
  ee_rechitenergy_20    = new TH1D("ee_rechitenergy_20","",1100,-2.,20.);
  ee_rechitenergy_22    = new TH1D("ee_rechitenergy_22","",1100,-2.,20.);
  ee_rechitenergy_24    = new TH1D("ee_rechitenergy_24","",1100,-2.,20.);
  ee_rechitenergy_26    = new TH1D("ee_rechitenergy_26","",1100,-2.,20.);
  ee_rechitenergy_28    = new TH1D("ee_rechitenergy_28","",1100,-2.,20.);
  ee_rechitenergy_30    = new TH1D("ee_rechitenergy_30","",1100,-2.,20.);
 

  eb_rechitet_02    = new TH1D("eb_rechitet_02","",1100,-2.,20.);
  eb_rechitet_04    = new TH1D("eb_rechitet_04","",1100,-2.,20.);
  eb_rechitet_06    = new TH1D("eb_rechitet_06","",1100,-2.,20.);
  eb_rechitet_08    = new TH1D("eb_rechitet_08","",1100,-2.,20.);
  eb_rechitet_10    = new TH1D("eb_rechitet_10","",1100,-2.,20.);
  eb_rechitet_12    = new TH1D("eb_rechitet_12","",1100,-2.,20.);
  eb_rechitet_14    = new TH1D("eb_rechitet_14","",1100,-2.,20.);
  eb_rechitet_148   = new TH1D("eb_rechitet_148","",1100,-2.,20.);
 
  ee_rechitet_16    = new TH1D("ee_rechitet_16","",1100,-2.,20.);
  ee_rechitet_18    = new TH1D("ee_rechitet_18","",1100,-2.,20.);
  ee_rechitet_20    = new TH1D("ee_rechitet_20","",1100,-2.,20.);
  ee_rechitet_22    = new TH1D("ee_rechitet_22","",1100,-2.,20.);
  ee_rechitet_24    = new TH1D("ee_rechitet_24","",1100,-2.,20.);
  ee_rechitet_26    = new TH1D("ee_rechitet_26","",1100,-2.,20.);
  ee_rechitet_28    = new TH1D("ee_rechitet_28","",1100,-2.,20.);
  ee_rechitet_30    = new TH1D("ee_rechitet_30","",1100,-2.,20.);
 


  eb_rechitetvspu_05 = new TH2D("eb_rechitetvspu_05","",25,0,50,500,-0.5,20);
  eb_rechitetvspu_10 = new TH2D("eb_rechitetvspu_10","",25,0,50,500,-0.5,20);
  eb_rechitetvspu_15 = new TH2D("eb_rechitetvspu_15","",25,0,50,500,-0.5,20);
 

  ee_rechitetvspu_20 = new TH2D("ee_rechitetvspu_20","",25,0,50,500,-0.5,20);
  ee_rechitetvspu_25 = new TH2D("ee_rechitetvspu_25","",25,0,50,500,-0.5,20);
  ee_rechitetvspu_30 = new TH2D("ee_rechitetvspu_30","",25,0,50,500,-0.5,20);
    

 




  eb_rechitet_   = new TH1D("eb_rechitet","",1100,-2.,20.);
  ee_rechitet_   = new TH1D("ee_rechitet","",1100,-2.,20.);

  eb_rechiten_vs_eta = new TH2D("ebrechiten_vs_eta","",30,-1.5,1.5,350,-2,5);
  eb_rechitet_vs_eta = new TH2D("ebrechitet_vs_eta","",30,-1.5,1.5,350,-2,5);


  eep_rechiten_vs_eta = new TH2D("eep_rechiten_vs_eta","",18,1.4,3.2,350,-2,5);
  eep_rechiten_vs_phi = new TH2D("eep_rechiten_vs_phi","",18,-3.1416,3.1416,350,-2,5);
  eem_rechiten_vs_eta = new TH2D("eem_rechiten_vs_eta","",18,1.4,3.2,350,-2,5);
  eem_rechiten_vs_phi = new TH2D("eem_rechiten_vs_phi","",18,-3.1416,3.1416,350,-2,5);


  eep_rechitet_vs_eta = new TH2D("eep_rechitet_vs_eta","",18,1.4,3.2,350,-2,5);
  eep_rechitet_vs_phi = new TH2D("eep_rechitet_vs_phi","",18,-3.1416,3.1416,350,-2,5);
  eem_rechitet_vs_eta = new TH2D("eem_rechitet_vs_eta","",18,1.4,3.2,350,-2,5);
  eem_rechitet_vs_phi = new TH2D("eem_rechitet_vs_phi","",18,-3.1416,3.1416,350,-2,5);

 
 
  ebocc = new TH2F("ebocc","",360,0,360,170,-85,85);
  eboccgt1 = new TH2F("eboccgt1","",360,0,360,170,-85,85);
  eboccgt1et = new TH2F("eboccgt1et","",360,0,360,170,-85,85);
  eboccet = new TH2F("eboccet","",360,0,360,170,-85,85);
  eboccetgt1et = new TH2F("eboccetgt1et","",360,0,360,170,-85,85);
  eboccen = new TH2F("eboccen","",360,0,360,170,-85,85);
  eboccengt1 = new TH2F("eboccengt1","",360,0,360,170,-85,85);


  eeocc = new TH2F("eeocc","",200,0,200,100,0,100);
  eeoccgt1 = new TH2F("eeoccgt1","",200,0,200,100,0,100);
  eeoccgt1et = new TH2F("eeoccgt1et","",200,0,200,100,0,100);
  eeoccet = new TH2F("eeoccet","",200,0,200,100,0,100);
  eeoccetgt1et = new TH2F("eeoccetgt1et","",200,0,200,100,0,100);
  eeoccen = new TH2F("eeoccen","",200,0,200,100,0,100);
  eeoccengt1 = new TH2F("eeoccengt1","",200,0,200,100,0,100);



  scocc_eb_gt50 = new TH2F("scocc_eb_gt50","",100,-3.14,3.14,100,-3,3);
  scocc_ee_gt50 = new TH2F("scocc_ee_gt50","",100,-3.14,3.14,100,-3,3);




  eb_timing_0   = new TH1D("eb_timing_0","",300,-50.,100.);
  eb_timing_200   = new TH1D("eb_timing_200","",300,-50.,100.);
  eb_timing_400   = new TH1D("eb_timing_400","",300,-50.,100.);
  eb_timing_600   = new TH1D("eb_timing_600","",300,-50.,100.);
  eb_timing_800   = new TH1D("eb_timing_800","",300,-50.,100.);
  eb_timing_1000   = new TH1D("eb_timing_1000","",300,-50.,100.);
  eb_timing_2000   = new TH1D("eb_timing_2000","",300,-50.,100.);
  eb_timing_3000   = new TH1D("eb_timing_3000","",300,-50.,100.);
  eb_timing_5000   = new TH1D("eb_timing_5000","",300,-50.,100.);


  eb_r4_0   = new TH1D("eb_r4_0","",200,0.2,1.2);
  eb_r4_200   = new TH1D("eb_r4_200","",200,0.2,1.2);
  eb_r4_400   = new TH1D("eb_r4_400","",200,0.2,1.2);
  eb_r4_600   = new TH1D("eb_r4_600","",200,0.2,1.2);
  eb_r4_800   = new TH1D("eb_r4_800","",200,0.2,1.2);
  eb_r4_1000   = new TH1D("eb_r4_1000","",200,0.2,1.2);
  eb_r4_2000   = new TH1D("eb_r4_2000","",200,0.2,1.2);
  eb_r4_3000   = new TH1D("eb_r4_3000","",200,0.2,1.2);
  eb_r4_5000   = new TH1D("eb_r4_5000","",200,0.2,1.2);

  eb_timing_r4_0   = new TH1D("eb_timing_r4_0","",300,-50.,100.);
  eb_timing_r4_200   = new TH1D("eb_timing_r4_200","",300,-50.,100.);
  eb_timing_r4_400   = new TH1D("eb_timing_r4_400","",300,-50.,100.);
  eb_timing_r4_600   = new TH1D("eb_timing_r4_600","",300,-50.,100.);
  eb_timing_r4_800   = new TH1D("eb_timing_r4_800","",300,-50.,100.);
  eb_timing_r4_1000   = new TH1D("eb_timing_r4_1000","",300,-50.,100.);
  eb_timing_r4_2000   = new TH1D("eb_timing_r4_2000","",300,-50.,100.);
  eb_timing_r4_3000   = new TH1D("eb_timing_r4_3000","",300,-50.,100.);
  eb_timing_r4_5000   = new TH1D("eb_timing_r4_5000","",300,-50.,100.);


  eb_timing_vs_r4_0   = new TH2D("eb_timing_vs_r4_0","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_200   = new TH2D("eb_timing_vs_r4_200","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_400   = new TH2D("eb_timing_vs_r4_400","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_600   = new TH2D("eb_timing_vs_r4_600","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_800   = new TH2D("eb_timing_vs_r4_800","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_1000   = new TH2D("eb_timing_vs_r4_1000","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_2000   = new TH2D("eb_timing_vs_r4_2000","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_3000   = new TH2D("eb_timing_vs_r4_3000","",100,0.2,1.2,150,-50.,100.);
  eb_timing_vs_r4_5000   = new TH2D("eb_timing_vs_r4_5000","",100,0.2,1.2,150,-50.,100.);


  numtp_vs_ieta=new TH1D("numtp_vs_ieta","",56,-28,28);
  numtp_vs_ieta_samp2=new TH1D("numtp_vs_ieta_samp2","",56,-28,28);
  numtp_geq1_vs_ieta=new TH1D("numtp_geq1_vs_ieta","",56,-28,28);
  numtp_etweighted_vs_ieta=new TH1D("numtp_etweighted_vs_ieta","",56,-28,28);


  rechiteta_vs_bxtrain_01 = new TH2D("rechiteta_vs_bxtrain_01","",40,-2,38,60,-3,3);
  rechiteta_vs_bxtrain_05 = new TH2D("rechiteta_vs_bxtrain_05","",40,-2,38,60,-3,3);
  sceta_vs_bxtrain = new TH2D("sceta_vs_bxtrain","",40,-2,38,60,-3,3);


  ebtime_vs_bxtrain_01 = new TH2D("ebtime_vx_bxtrain_01","",40,-2,38,200,-100,100);
  ebtime_vs_bxtrain_05 = new TH2D("ebtime_vx_bxtrain_05","",40,-2,38,200,-100,100);

  eetime_vs_bxtrain_01 = new TH2D("eetime_vx_bxtrain_01","",40,-2,38,200,-100,100);
  eetime_vs_bxtrain_05 = new TH2D("eetime_vx_bxtrain_05","",40,-2,38,200,-100,100);



  matched_timevset=new TH2D("matched_timevset","",120,-60,60,100,-0.5,9.5);
  unmatched_timevset=new TH2D("unmatched_timevset","",120,-60,60,100,-0.5,9.5);

  matched_timevsieta=new TH2D("matched_timevsieta","",120,-60,60,11,17.5,28.5);
  unmatched_timevsieta=new TH2D("unmatched_timevsieta","",120,-60,60,11,17.5,28.5);

  matched_timevsieta_etweighted=new TH2D("matched_timevsieta_etweighted","",120,-60,60,11,17.5,28.5);
  unmatched_timevsieta_etweighted=new TH2D("unmatched_timevsieta_etweighted","",120,-60,60,11,17.5,28.5);

  matched_timevsetvsieta=new TH3D("matched_timevsetvsieta","",60,-60,60,100,-0.5,9.5,11,17.5,28.5);
  unmatched_timevsetvsieta=new TH3D("unmatched_timevsetvsieta","",60,-60,60,100,-0.5,9.5,11,17.5,28.5);

  ietamatched=new TH1D("ietamatched","",57,-28.5,28.5);
  ietatotal=new TH1D("ietatotal","",57,-28.5,28.5);




  rechiteta_all=new TH1F("rechiteta_all","",100,-3,3);
  rechiteta_gt1et=new TH1F("rechiteta_gt1et","",100,-3,3);
  rechiteta_etweight=new TH1F("rechiteta_etweight","",100,-3,3);
  rechiteta_etweight_gt1et=new TH1F("rechiteta_etweight_gt1et","",100,-3,3);


  rechiteta_gt1et_pu10=new TH1F("rechiteta_gt1et_pu10","",100,-3,3);
  rechiteta_gt1et_pu20=new TH1F("rechiteta_gt1et_pu20","",100,-3,3);
  rechiteta_gt1et_pu30=new TH1F("rechiteta_gt1et_pu30","",100,-3,3);





  calotowereta_all=new TH1F("calotowereta_all","",100,-3,3);
  calotowereta_gt1et=new TH1F("calotowereta_gt1et","",100,-3,3);
  calotowereta_etweight=new TH1F("calotowereta_etweight","",100,-3,3);
  calotowereta_etweight_gt1et=new TH1F("calotowereta_etweight_gt1et","",100,-3,3);


  calotowereta_gt1et_pu10=new TH1F("calotowereta_gt1et_pu10","",100,-3,3);
  calotowereta_gt1et_pu20=new TH1F("calotowereta_gt1et_pu20","",100,-3,3);
  calotowereta_gt1et_pu30=new TH1F("calotowereta_gt1et_pu30","",100,-3,3);



  sceta_all=new TH1F("sceta_all","",100,-3,3);
  sceta_severity0=new TH1F("sceta_severity0","",100,-3,3);
  sceta_koot0=new TH1F("sceta_koot0","",100,-3,3);


  sceta_all_gt2=new TH1F("sceta_all_gt2","",100,-3,3);
  sceta_severity0_gt2=new TH1F("sceta_severity0_gt2","",100,-3,3);
  sceta_koot0_gt2=new TH1F("sceta_koot0_gt2","",100,-3,3);


  sceta_all_gt5=new TH1F("sceta_all_gt5","",100,-3,3);
  sceta_severity0_gt5=new TH1F("sceta_severity0_gt5","",100,-3,3);
  sceta_koot0_gt5=new TH1F("sceta_koot0_gt5","",100,-3,3);

  sceta_all_gt10=new TH1F("sceta_all_gt10","",100,-3,3);
  sceta_severity0_gt10=new TH1F("sceta_severity0_gt10","",100,-3,3);
  sceta_koot0_gt10=new TH1F("sceta_koot0_gt10","",100,-3,3);
  


  sceta_all_pueq01=new TH1F("sceta_all_pueq01","",100,-3,3);
  sceta_severity0_pueq01=new TH1F("sceta_severity0_pueq01","",100,-3,3);

  sceta_all_pueq02=new TH1F("sceta_all_pueq02","",100,-3,3);
  sceta_severity0_pueq02=new TH1F("sceta_severity0_pueq02","",100,-3,3);

  sceta_all_pueq03=new TH1F("sceta_all_pueq03","",100,-3,3);
  sceta_severity0_pueq03=new TH1F("sceta_severity0_pueq03","",100,-3,3);

  sceta_all_pueq04=new TH1F("sceta_all_pueq04","",100,-3,3);
  sceta_severity0_pueq04=new TH1F("sceta_severity0_pueq04","",100,-3,3);

  sceta_all_pueq05=new TH1F("sceta_all_pueq05","",100,-3,3);
  sceta_severity0_pueq05=new TH1F("sceta_severity0_pueq05","",100,-3,3);

  sceta_all_pueq06=new TH1F("sceta_all_pueq06","",100,-3,3);
  sceta_severity0_pueq06=new TH1F("sceta_severity0_pueq06","",100,-3,3);

  sceta_all_pueq07=new TH1F("sceta_all_pueq07","",100,-3,3);
  sceta_severity0_pueq07=new TH1F("sceta_severity0_pueq07","",100,-3,3);

  sceta_all_pueq08=new TH1F("sceta_all_pueq08","",100,-3,3);
  sceta_severity0_pueq08=new TH1F("sceta_severity0_pueq08","",100,-3,3);

  sceta_all_pueq09=new TH1F("sceta_all_pueq09","",100,-3,3);
  sceta_severity0_pueq09=new TH1F("sceta_severity0_pueq09","",100,-3,3);



  scet_eb_all=new TH1F("scet_eb_all","",200,0,100);
  scet_eb_severity0=new TH1F("scet_eb_severity0","",200,0,100);
  scet_eb_koot0=new TH1F("scet_eb_koot0","",200,0,100);

  scet_ee_all=new TH1F("scet_ee_all","",200,0,100);
  scet_ee_severity0=new TH1F("scet_ee_severity0","",200,0,100);
  scet_ee_koot0=new TH1F("scet_ee_koot0","",200,0,100);



  scet_eb_all_eta15=new TH1F("scet_eb_all_eta15","",200,0,100);
  scet_eb_all_eta20=new TH1F("scet_eb_all_eta20","",200,0,100);
  scet_eb_all_eta25=new TH1F("scet_eb_all_eta25","",200,0,100);


  scet_eb_all_eta15_pu10=new TH1F("scet_eb_all_eta15_pu10","",200,0,100);
  scet_eb_all_eta20_pu10=new TH1F("scet_eb_all_eta20_pu10","",200,0,100);
  scet_eb_all_eta25_pu10=new TH1F("scet_eb_all_eta25_pu10","",200,0,100);

  scet_eb_all_eta15_pu20=new TH1F("scet_eb_all_eta15_pu20","",200,0,100);
  scet_eb_all_eta20_pu20=new TH1F("scet_eb_all_eta20_pu20","",200,0,100);
  scet_eb_all_eta25_pu20=new TH1F("scet_eb_all_eta25_pu20","",200,0,100);

  scet_eb_all_eta15_pu30=new TH1F("scet_eb_all_eta15_pu30","",200,0,100);
  scet_eb_all_eta20_pu30=new TH1F("scet_eb_all_eta20_pu30","",200,0,100);
  scet_eb_all_eta25_pu30=new TH1F("scet_eb_all_eta25_pu30","",200,0,100);


  scet_eb_all_eta15_pueq10=new TH1F("scet_eb_all_eta15_pueq10","",200,0,100);
  scet_eb_all_eta20_pueq10=new TH1F("scet_eb_all_eta20_pueq10","",200,0,100);
  scet_eb_all_eta25_pueq10=new TH1F("scet_eb_all_eta25_pueq10","",200,0,100);

  scet_eb_all_eta15_pueq20=new TH1F("scet_eb_all_eta15_pueq20","",200,0,100);
  scet_eb_all_eta20_pueq20=new TH1F("scet_eb_all_eta20_pueq20","",200,0,100);
  scet_eb_all_eta25_pueq20=new TH1F("scet_eb_all_eta25_pueq20","",200,0,100);


  
  scet_ee_all_eta15=new TH1F("scet_ee_all_eta15","",200,0,100);
  scet_ee_all_eta20=new TH1F("scet_ee_all_eta20","",200,0,100);
  scet_ee_all_eta25=new TH1F("scet_ee_all_eta25","",200,0,100);


  scet_ee_all_eta15_pu10=new TH1F("scet_ee_all_eta15_pu10","",200,0,100);
  scet_ee_all_eta20_pu10=new TH1F("scet_ee_all_eta20_pu10","",200,0,100);
  scet_ee_all_eta25_pu10=new TH1F("scet_ee_all_eta25_pu10","",200,0,100);

  scet_ee_all_eta15_pu20=new TH1F("scet_ee_all_eta15_pu20","",200,0,100);
  scet_ee_all_eta20_pu20=new TH1F("scet_ee_all_eta20_pu20","",200,0,100);
  scet_ee_all_eta25_pu20=new TH1F("scet_ee_all_eta25_pu20","",200,0,100);

  scet_ee_all_eta15_pu30=new TH1F("scet_ee_all_eta15_pu30","",200,0,100);
  scet_ee_all_eta20_pu30=new TH1F("scet_ee_all_eta20_pu30","",200,0,100);
  scet_ee_all_eta25_pu30=new TH1F("scet_ee_all_eta25_pu30","",200,0,100);


  scet_ee_all_eta15_pueq10=new TH1F("scet_ee_all_eta15_pueq10","",200,0,100);
  scet_ee_all_eta20_pueq10=new TH1F("scet_ee_all_eta20_pueq10","",200,0,100);
  scet_ee_all_eta25_pueq10=new TH1F("scet_ee_all_eta25_pueq10","",200,0,100);

  scet_ee_all_eta15_pueq20=new TH1F("scet_ee_all_eta15_pueq20","",200,0,100);
  scet_ee_all_eta20_pueq20=new TH1F("scet_ee_all_eta20_pueq20","",200,0,100);
  scet_ee_all_eta25_pueq20=new TH1F("scet_ee_all_eta25_pueq20","",200,0,100);



  scsumet_eb_all_eta15=new TH1F("scsumet_eb_all_eta15","",200,0,200);
  scsumet_eb_all_eta20=new TH1F("scsumet_eb_all_eta20","",200,0,200);
  scsumet_eb_all_eta25=new TH1F("scsumet_eb_all_eta25","",200,0,200);


  scsumet_eb_all_eta15_pu10=new TH1F("scsumet_eb_all_eta15_pu10","",200,0,200);
  scsumet_eb_all_eta20_pu10=new TH1F("scsumet_eb_all_eta20_pu10","",200,0,200);
  scsumet_eb_all_eta25_pu10=new TH1F("scsumet_eb_all_eta25_pu10","",200,0,200);

  scsumet_eb_all_eta15_pu20=new TH1F("scsumet_eb_all_eta15_pu20","",200,0,200);
  scsumet_eb_all_eta20_pu20=new TH1F("scsumet_eb_all_eta20_pu20","",200,0,200);
  scsumet_eb_all_eta25_pu20=new TH1F("scsumet_eb_all_eta25_pu20","",200,0,200);

  scsumet_eb_all_eta15_pu30=new TH1F("scsumet_eb_all_eta15_pu30","",200,0,200);
  scsumet_eb_all_eta20_pu30=new TH1F("scsumet_eb_all_eta20_pu30","",200,0,200);
  scsumet_eb_all_eta25_pu30=new TH1F("scsumet_eb_all_eta25_pu30","",200,0,200);


  scsumet_eb_all_eta15_pueq10=new TH1F("scsumet_eb_all_eta15_pueq10","",200,0,200);
  scsumet_eb_all_eta20_pueq10=new TH1F("scsumet_eb_all_eta20_pueq10","",200,0,200);
  scsumet_eb_all_eta25_pueq10=new TH1F("scsumet_eb_all_eta25_pueq10","",200,0,200);

  scsumet_eb_all_eta15_pueq20=new TH1F("scsumet_eb_all_eta15_pueq20","",200,0,200);
  scsumet_eb_all_eta20_pueq20=new TH1F("scsumet_eb_all_eta20_pueq20","",200,0,200);
  scsumet_eb_all_eta25_pueq20=new TH1F("scsumet_eb_all_eta25_pueq20","",200,0,200);


  
  scsumet_ee_all_eta15=new TH1F("scsumet_ee_all_eta15","",200,0,200);
  scsumet_ee_all_eta20=new TH1F("scsumet_ee_all_eta20","",200,0,200);
  scsumet_ee_all_eta25=new TH1F("scsumet_ee_all_eta25","",200,0,200);


  scsumet_ee_all_eta15_pu10=new TH1F("scsumet_ee_all_eta15_pu10","",200,0,200);
  scsumet_ee_all_eta20_pu10=new TH1F("scsumet_ee_all_eta20_pu10","",200,0,200);
  scsumet_ee_all_eta25_pu10=new TH1F("scsumet_ee_all_eta25_pu10","",200,0,200);

  scsumet_ee_all_eta15_pu20=new TH1F("scsumet_ee_all_eta15_pu20","",200,0,200);
  scsumet_ee_all_eta20_pu20=new TH1F("scsumet_ee_all_eta20_pu20","",200,0,200);
  scsumet_ee_all_eta25_pu20=new TH1F("scsumet_ee_all_eta25_pu20","",200,0,200);

  scsumet_ee_all_eta15_pu30=new TH1F("scsumet_ee_all_eta15_pu30","",200,0,200);
  scsumet_ee_all_eta20_pu30=new TH1F("scsumet_ee_all_eta20_pu30","",200,0,200);
  scsumet_ee_all_eta25_pu30=new TH1F("scsumet_ee_all_eta25_pu30","",200,0,200);


  scsumet_ee_all_eta15_pueq10=new TH1F("scsumet_ee_all_eta15_pueq10","",200,0,200);
  scsumet_ee_all_eta20_pueq10=new TH1F("scsumet_ee_all_eta20_pueq10","",200,0,200);
  scsumet_ee_all_eta25_pueq10=new TH1F("scsumet_ee_all_eta25_pueq10","",200,0,200);

  scsumet_ee_all_eta15_pueq20=new TH1F("scsumet_ee_all_eta15_pueq20","",200,0,200);
  scsumet_ee_all_eta20_pueq20=new TH1F("scsumet_ee_all_eta20_pueq20","",200,0,200);
  scsumet_ee_all_eta25_pueq20=new TH1F("scsumet_ee_all_eta25_pueq20","",200,0,200);




  scet_eb_all_eta15_pueq01=new TH1F("scet_eb_all_eta15_pueq01","",200,0,100);
  scet_eb_all_eta20_pueq01=new TH1F("scet_eb_all_eta20_pueq01","",200,0,100);
  scet_eb_all_eta25_pueq01=new TH1F("scet_eb_all_eta25_pueq01","",200,0,100);

  scet_ee_all_eta15_pueq01=new TH1F("scet_ee_all_eta15_pueq01","",200,0,100);
  scet_ee_all_eta20_pueq01=new TH1F("scet_ee_all_eta20_pueq01","",200,0,100);
  scet_ee_all_eta25_pueq01=new TH1F("scet_ee_all_eta25_pueq01","",200,0,100);

  scsumet_eb_all_eta15_pueq01=new TH1F("scsumet_eb_all_eta15_pueq01","",200,0,200);
  scsumet_eb_all_eta20_pueq01=new TH1F("scsumet_eb_all_eta20_pueq01","",200,0,200);
  scsumet_eb_all_eta25_pueq01=new TH1F("scsumet_eb_all_eta25_pueq01","",200,0,200);

  scsumet_ee_all_eta15_pueq01=new TH1F("scsumet_ee_all_eta15_pueq01","",200,0,200);
  scsumet_ee_all_eta20_pueq01=new TH1F("scsumet_ee_all_eta20_pueq01","",200,0,200);
  scsumet_ee_all_eta25_pueq01=new TH1F("scsumet_ee_all_eta25_pueq01","",200,0,200);


 scet_eb_all_eta15_pueq02=new TH1F("scet_eb_all_eta15_pueq02","",200,0,100);
  scet_eb_all_eta20_pueq02=new TH1F("scet_eb_all_eta20_pueq02","",200,0,100);
  scet_eb_all_eta25_pueq02=new TH1F("scet_eb_all_eta25_pueq02","",200,0,100);

  scet_ee_all_eta15_pueq02=new TH1F("scet_ee_all_eta15_pueq02","",200,0,100);
  scet_ee_all_eta20_pueq02=new TH1F("scet_ee_all_eta20_pueq02","",200,0,100);
  scet_ee_all_eta25_pueq02=new TH1F("scet_ee_all_eta25_pueq02","",200,0,100);

  scsumet_eb_all_eta15_pueq02=new TH1F("scsumet_eb_all_eta15_pueq02","",200,0,200);
  scsumet_eb_all_eta20_pueq02=new TH1F("scsumet_eb_all_eta20_pueq02","",200,0,200);
  scsumet_eb_all_eta25_pueq02=new TH1F("scsumet_eb_all_eta25_pueq02","",200,0,200);

  scsumet_ee_all_eta15_pueq02=new TH1F("scsumet_ee_all_eta15_pueq02","",200,0,200);
  scsumet_ee_all_eta20_pueq02=new TH1F("scsumet_ee_all_eta20_pueq02","",200,0,200);
  scsumet_ee_all_eta25_pueq02=new TH1F("scsumet_ee_all_eta25_pueq02","",200,0,200);


 scet_eb_all_eta15_pueq03=new TH1F("scet_eb_all_eta15_pueq03","",200,0,100);
  scet_eb_all_eta20_pueq03=new TH1F("scet_eb_all_eta20_pueq03","",200,0,100);
  scet_eb_all_eta25_pueq03=new TH1F("scet_eb_all_eta25_pueq03","",200,0,100);

  scet_ee_all_eta15_pueq03=new TH1F("scet_ee_all_eta15_pueq03","",200,0,100);
  scet_ee_all_eta20_pueq03=new TH1F("scet_ee_all_eta20_pueq03","",200,0,100);
  scet_ee_all_eta25_pueq03=new TH1F("scet_ee_all_eta25_pueq03","",200,0,100);

  scsumet_eb_all_eta15_pueq03=new TH1F("scsumet_eb_all_eta15_pueq03","",200,0,200);
  scsumet_eb_all_eta20_pueq03=new TH1F("scsumet_eb_all_eta20_pueq03","",200,0,200);
  scsumet_eb_all_eta25_pueq03=new TH1F("scsumet_eb_all_eta25_pueq03","",200,0,200);

  scsumet_ee_all_eta15_pueq03=new TH1F("scsumet_ee_all_eta15_pueq03","",200,0,200);
  scsumet_ee_all_eta20_pueq03=new TH1F("scsumet_ee_all_eta20_pueq03","",200,0,200);
  scsumet_ee_all_eta25_pueq03=new TH1F("scsumet_ee_all_eta25_pueq03","",200,0,200);


 scet_eb_all_eta15_pueq04=new TH1F("scet_eb_all_eta15_pueq04","",200,0,100);
  scet_eb_all_eta20_pueq04=new TH1F("scet_eb_all_eta20_pueq04","",200,0,100);
  scet_eb_all_eta25_pueq04=new TH1F("scet_eb_all_eta25_pueq04","",200,0,100);

  scet_ee_all_eta15_pueq04=new TH1F("scet_ee_all_eta15_pueq04","",200,0,100);
  scet_ee_all_eta20_pueq04=new TH1F("scet_ee_all_eta20_pueq04","",200,0,100);
  scet_ee_all_eta25_pueq04=new TH1F("scet_ee_all_eta25_pueq04","",200,0,100);

  scsumet_eb_all_eta15_pueq04=new TH1F("scsumet_eb_all_eta15_pueq04","",200,0,200);
  scsumet_eb_all_eta20_pueq04=new TH1F("scsumet_eb_all_eta20_pueq04","",200,0,200);
  scsumet_eb_all_eta25_pueq04=new TH1F("scsumet_eb_all_eta25_pueq04","",200,0,200);

  scsumet_ee_all_eta15_pueq04=new TH1F("scsumet_ee_all_eta15_pueq04","",200,0,200);
  scsumet_ee_all_eta20_pueq04=new TH1F("scsumet_ee_all_eta20_pueq04","",200,0,200);
  scsumet_ee_all_eta25_pueq04=new TH1F("scsumet_ee_all_eta25_pueq04","",200,0,200);


 scet_eb_all_eta15_pueq05=new TH1F("scet_eb_all_eta15_pueq05","",200,0,100);
  scet_eb_all_eta20_pueq05=new TH1F("scet_eb_all_eta20_pueq05","",200,0,100);
  scet_eb_all_eta25_pueq05=new TH1F("scet_eb_all_eta25_pueq05","",200,0,100);

  scet_ee_all_eta15_pueq05=new TH1F("scet_ee_all_eta15_pueq05","",200,0,100);
  scet_ee_all_eta20_pueq05=new TH1F("scet_ee_all_eta20_pueq05","",200,0,100);
  scet_ee_all_eta25_pueq05=new TH1F("scet_ee_all_eta25_pueq05","",200,0,100);

  scsumet_eb_all_eta15_pueq05=new TH1F("scsumet_eb_all_eta15_pueq05","",200,0,200);
  scsumet_eb_all_eta20_pueq05=new TH1F("scsumet_eb_all_eta20_pueq05","",200,0,200);
  scsumet_eb_all_eta25_pueq05=new TH1F("scsumet_eb_all_eta25_pueq05","",200,0,200);

  scsumet_ee_all_eta15_pueq05=new TH1F("scsumet_ee_all_eta15_pueq05","",200,0,200);
  scsumet_ee_all_eta20_pueq05=new TH1F("scsumet_ee_all_eta20_pueq05","",200,0,200);
  scsumet_ee_all_eta25_pueq05=new TH1F("scsumet_ee_all_eta25_pueq05","",200,0,200);


 scet_eb_all_eta15_pueq06=new TH1F("scet_eb_all_eta15_pueq06","",200,0,100);
  scet_eb_all_eta20_pueq06=new TH1F("scet_eb_all_eta20_pueq06","",200,0,100);
  scet_eb_all_eta25_pueq06=new TH1F("scet_eb_all_eta25_pueq06","",200,0,100);

  scet_ee_all_eta15_pueq06=new TH1F("scet_ee_all_eta15_pueq06","",200,0,100);
  scet_ee_all_eta20_pueq06=new TH1F("scet_ee_all_eta20_pueq06","",200,0,100);
  scet_ee_all_eta25_pueq06=new TH1F("scet_ee_all_eta25_pueq06","",200,0,100);

  scsumet_eb_all_eta15_pueq06=new TH1F("scsumet_eb_all_eta15_pueq06","",200,0,200);
  scsumet_eb_all_eta20_pueq06=new TH1F("scsumet_eb_all_eta20_pueq06","",200,0,200);
  scsumet_eb_all_eta25_pueq06=new TH1F("scsumet_eb_all_eta25_pueq06","",200,0,200);

  scsumet_ee_all_eta15_pueq06=new TH1F("scsumet_ee_all_eta15_pueq06","",200,0,200);
  scsumet_ee_all_eta20_pueq06=new TH1F("scsumet_ee_all_eta20_pueq06","",200,0,200);
  scsumet_ee_all_eta25_pueq06=new TH1F("scsumet_ee_all_eta25_pueq06","",200,0,200);


 scet_eb_all_eta15_pueq07=new TH1F("scet_eb_all_eta15_pueq07","",200,0,100);
  scet_eb_all_eta20_pueq07=new TH1F("scet_eb_all_eta20_pueq07","",200,0,100);
  scet_eb_all_eta25_pueq07=new TH1F("scet_eb_all_eta25_pueq07","",200,0,100);

  scet_ee_all_eta15_pueq07=new TH1F("scet_ee_all_eta15_pueq07","",200,0,100);
  scet_ee_all_eta20_pueq07=new TH1F("scet_ee_all_eta20_pueq07","",200,0,100);
  scet_ee_all_eta25_pueq07=new TH1F("scet_ee_all_eta25_pueq07","",200,0,100);

  scsumet_eb_all_eta15_pueq07=new TH1F("scsumet_eb_all_eta15_pueq07","",200,0,200);
  scsumet_eb_all_eta20_pueq07=new TH1F("scsumet_eb_all_eta20_pueq07","",200,0,200);
  scsumet_eb_all_eta25_pueq07=new TH1F("scsumet_eb_all_eta25_pueq07","",200,0,200);

  scsumet_ee_all_eta15_pueq07=new TH1F("scsumet_ee_all_eta15_pueq07","",200,0,200);
  scsumet_ee_all_eta20_pueq07=new TH1F("scsumet_ee_all_eta20_pueq07","",200,0,200);
  scsumet_ee_all_eta25_pueq07=new TH1F("scsumet_ee_all_eta25_pueq07","",200,0,200);


 scet_eb_all_eta15_pueq08=new TH1F("scet_eb_all_eta15_pueq08","",200,0,100);
  scet_eb_all_eta20_pueq08=new TH1F("scet_eb_all_eta20_pueq08","",200,0,100);
  scet_eb_all_eta25_pueq08=new TH1F("scet_eb_all_eta25_pueq08","",200,0,100);

  scet_ee_all_eta15_pueq08=new TH1F("scet_ee_all_eta15_pueq08","",200,0,100);
  scet_ee_all_eta20_pueq08=new TH1F("scet_ee_all_eta20_pueq08","",200,0,100);
  scet_ee_all_eta25_pueq08=new TH1F("scet_ee_all_eta25_pueq08","",200,0,100);

  scsumet_eb_all_eta15_pueq08=new TH1F("scsumet_eb_all_eta15_pueq08","",200,0,200);
  scsumet_eb_all_eta20_pueq08=new TH1F("scsumet_eb_all_eta20_pueq08","",200,0,200);
  scsumet_eb_all_eta25_pueq08=new TH1F("scsumet_eb_all_eta25_pueq08","",200,0,200);

  scsumet_ee_all_eta15_pueq08=new TH1F("scsumet_ee_all_eta15_pueq08","",200,0,200);
  scsumet_ee_all_eta20_pueq08=new TH1F("scsumet_ee_all_eta20_pueq08","",200,0,200);
  scsumet_ee_all_eta25_pueq08=new TH1F("scsumet_ee_all_eta25_pueq08","",200,0,200);


  scet_eb_all_eta15_pueq09=new TH1F("scet_eb_all_eta15_pueq09","",200,0,100);
  scet_eb_all_eta20_pueq09=new TH1F("scet_eb_all_eta20_pueq09","",200,0,100);
  scet_eb_all_eta25_pueq09=new TH1F("scet_eb_all_eta25_pueq09","",200,0,100);

  scet_ee_all_eta15_pueq09=new TH1F("scet_ee_all_eta15_pueq09","",200,0,100);
  scet_ee_all_eta20_pueq09=new TH1F("scet_ee_all_eta20_pueq09","",200,0,100);
  scet_ee_all_eta25_pueq09=new TH1F("scet_ee_all_eta25_pueq09","",200,0,100);

  scsumet_eb_all_eta15_pueq09=new TH1F("scsumet_eb_all_eta15_pueq09","",200,0,200);
  scsumet_eb_all_eta20_pueq09=new TH1F("scsumet_eb_all_eta20_pueq09","",200,0,200);
  scsumet_eb_all_eta25_pueq09=new TH1F("scsumet_eb_all_eta25_pueq09","",200,0,200);

  scsumet_ee_all_eta15_pueq09=new TH1F("scsumet_ee_all_eta15_pueq09","",200,0,200);
  scsumet_ee_all_eta20_pueq09=new TH1F("scsumet_ee_all_eta20_pueq09","",200,0,200);
  scsumet_ee_all_eta25_pueq09=new TH1F("scsumet_ee_all_eta25_pueq09","",200,0,200);



  scsumet_eb_all=new TH1F("scsumet_eb_all","",200,0,200);
  scsumet_eb_severity0=new TH1F("scsumet_eb_severity0","",200,0,200);
  scsumet_eb_koot0=new TH1F("scsumet_eb_koot0","",200,0,200);

  scsumet_ee_all=new TH1F("scsumet_ee_all","",200,0,200);
  scsumet_ee_severity0=new TH1F("scsumet_ee_severity0","",200,0,200);
  scsumet_ee_koot0=new TH1F("scsumet_ee_koot0","",200,0,200);


  scsumet_eb_all_gt2=new TH1F("scsumet_eb_all_gt2","",200,0,200);
  scsumet_eb_severity0_gt2=new TH1F("scsumet_eb_severity0_gt2","",200,0,200);
  scsumet_eb_koot0_gt2=new TH1F("scsumet_eb_koot0_gt2","",200,0,200);

  scsumet_ee_all_gt2=new TH1F("scsumet_ee_all_gt2","",200,0,200);
  scsumet_ee_severity0_gt2=new TH1F("scsumet_ee_severity0_gt2","",200,0,200);
  scsumet_ee_koot0_gt2=new TH1F("scsumet_ee_koot0_gt2","",200,0,200);



  scsumet_eb_all_gt5=new TH1F("scsumet_eb_all_gt5","",200,0,200);
  scsumet_eb_severity0_gt5=new TH1F("scsumet_eb_severity0_gt5","",200,0,200);
  scsumet_eb_koot0_gt5=new TH1F("scsumet_eb_koot0_gt5","",200,0,200);

  scsumet_ee_all_gt5=new TH1F("scsumet_ee_all_gt5","",200,0,200);
  scsumet_ee_severity0_gt5=new TH1F("scsumet_ee_severity0_gt5","",200,0,200);
  scsumet_ee_koot0_gt5=new TH1F("scsumet_ee_koot0_gt5","",200,0,200);



  scsumet_eb_all_gt10=new TH1F("scsumet_eb_all_gt10","",200,0,200);
  scsumet_eb_severity0_gt10=new TH1F("scsumet_eb_severity0_gt10","",200,0,200);
  scsumet_eb_koot0_gt10=new TH1F("scsumet_eb_koot0_gt10","",200,0,200);

  scsumet_ee_all_gt10=new TH1F("scsumet_ee_all_gt10","",200,0,200);
  scsumet_ee_severity0_gt10=new TH1F("scsumet_ee_severity0_gt10","",200,0,200);
  scsumet_ee_koot0_gt10=new TH1F("scsumet_ee_koot0_gt10","",200,0,200);



  calotower_map=new TH2F("calotower_map","",72,0.5,72.5,57,-28.5,28.5);
  tp_map=new TH2F("tp_map","",72,0.5,72.5,57,-28.5,28.5);
  tp_map_samp2=new TH2F("tp_map","",72,0.5,72.5,57,-28.5,28.5);
  tp_map2=new TH2F("tp_map2","",72,0.5,72.5,57,-28.5,28.5);
  rechit_map=new TH2F("rechit_map","",72,0.5,72.5,57,-28.5,28.5);
  rechit_map_time=new TH2F("rechit_map_time","",72,0.5,72.5,57,-28.5,28.5);
  rechit_map_ixiy=new TH2F("rechit_map_ixiy","",200,0,200,100,0,100);
  rechit_map_time_ixiy=new TH2F("rechit_map_time_ixiy","",200,0,200,100,0,100);


  calotower_map_cumul=new TH2F("calotower_map_cumul","",72,0.5,72.5,57,-28.5,28.5);
  tpmap_cumul=new TH2F("tpmap_cumul","",72,0.5,72.5,57,-28.5,28.5);
  tpmap_cumul2=new TH2F("tpmap_cumul2","",72,0.5,72.5,57,-28.5,28.5);
  rechitmap_cumul=new TH2F("rechitmap_cumul","",72,0.5,72.5,57,-28.5,28.5);

  tower_map_cumul_matched=new TH2F("tower_map_cumul_matched","",72,0.5,72.5,57,-28.5,28.5);
  tower_map_cumul_matched2=new TH2F("tower_map_cumul_matched2","",72,0.5,72.5,57,-28.5,28.5);
  tower_map_cumul_matched_rh=new TH2F("tower_map_cumul_matched_rh","",72,0.5,72.5,57,-28.5,28.5);
  tower_map_cumul_matched_rh2=new TH2F("tower_map_cumul_matched_rh2","",72,0.5,72.5,57,-28.5,28.5);

  tower_map_cumul_unmatched1=new TH2F("tower_map_cumul_unmatched1","",72,0.5,72.5,57,-28.5,28.5);

  tower_map_cumul_unmatched2=new TH2F("tower_map_cumul_unmatched2","",72,0.5,72.5,57,-28.5,28.5);



  ct_et_all=new TH1D("ct_et_all","",41,-0.25,20.25);
  ct_et_matched=new TH1D("ct_et_matched","",41,-0.25,20.25);
  ct_et_matched2=new TH1D("ct_et_matched2","",41,-0.25,20.25);
  ct_et_unmatched=new TH1D("ct_et_unmatched","",41,-0.25,20.25);

  rh_et_matched_rh=new TH1D("rh_et_matched_rh","",41,-0.25,20.25);
  rh_et_matched_rh2=new TH1D("rh_et_matched_rh2","",41,-0.25,20.25);


  tp_et_all=new TH1D("tp_et_all","",41,-0.25,20.25);
  tp_et_matched=new TH1D("tp_et_matched","",41,-0.25,20.25);
  tp_et_matched2=new TH1D("tp_et_matched2","",41,-0.25,20.25);
 
  tp_et_matched_rh=new TH1D("tp_et_matched_rh","",41,-0.25,20.25);
  tp_et_matched_rh2=new TH1D("tp_et_matched_rh2","",41,-0.25,20.25);
 
  tp_et_unmatched=new TH1D("tp_et_unmatched","",41,-0.25,20.25);

  ctminustp_et_matched=new TH1D("ctminustp_et_matched","",31,-3.1,3.1);
  ctminustp_et_matched2=new TH1D("ctminustp_et_matched2","",31,-3.1,3.1);
 
  rhminustp_et_matched_rh=new TH1D("rhminustp_et_matched_rh","",31,-3.1,3.1);
  rhminustp_et_matched_rh2=new TH1D("rhminustp_et_matched_rh2","",31,-3.1,3.1);

  ct_et_all_eb=new TH1D("ct_et_all_eb","",41,-0.25,20.25);
  ct_et_matched_eb=new TH1D("ct_et_matched_eb","",41,-0.25,20.25);
  ct_et_matched_eb2=new TH1D("ct_et_matched_eb2","",41,-0.25,20.25);


  rh_et_matched_eb_rh=new TH1D("rh_et_matched_eb_rh","",41,-0.25,20.25);
  rh_et_matched_eb_rh2=new TH1D("rh_et_matched_eb_rh2","",41,-0.25,20.25);

  ct_et_unmatched_eb=new TH1D("ct_et_unmatched_eb","",41,-0.25,20.25);


  tp_et_all_eb=new TH1D("tp_et_all_eb","",41,-0.25,20.25);
  tp_et_matched_eb=new TH1D("tp_et_matched_eb","",41,-0.25,20.25);
  tp_et_matched_eb2=new TH1D("tp_et_matched_eb2","",41,-0.25,20.25);
  tp_et_matched_eb_rh=new TH1D("tp_et_matched_eb_rh","",41,-0.25,20.25);
  tp_et_matched_eb_rh2=new TH1D("tp_et_matched_eb_rh2","",41,-0.25,20.25);
  tp_et_unmatched_eb=new TH1D("tp_et_unmatched_eb","",41,-0.25,20.25);


  ctminustp_et_matched_eb=new TH1D("ctminustp_et_matched_eb","",31,-3.1,3.1);
  ctminustp_et_matched_eb2=new TH1D("ctminustp_et_matched_eb2","",31,-3.1,3.1);

  rhminustp_et_matched_eb_rh=new TH1D("rhminustp_et_matched_eb_rh","",31,-3.1,3.1);
  rhminustp_et_matched_eb_rh2=new TH1D("rhminustp_et_matched_eb_rh2","",31,-3.1,3.1);


  ct_et_all_ee=new TH1D("ct_et_all_ee","",41,-0.25,20.25);
  ct_et_matched_ee=new TH1D("ct_et_matched_ee","",41,-0.25,20.25);
  ct_et_matched_ee2=new TH1D("ct_et_matched_ee2","",41,-0.25,20.25);
  ct_et_unmatched_ee=new TH1D("ct_et_unmatched_ee","",41,-0.25,20.25);



  rh_et_matched_ee_rh=new TH1D("rh_et_matched_ee_rh","",41,-0.25,20.25);
  rh_et_matched_ee24_rh=new TH1D("rh_et_matched_ee24_rh","",41,-0.25,20.25);
  rh_et_matched_ee_rh2=new TH1D("rh_et_matched_ee_rh2","",41,-0.25,20.25);

 
  rh_et_all_ee=new TH1D("rh_et_all_ee","",41,-0.25,20.25);
  rh_et_all_ee24=new TH1D("rh_et_all_ee24","",41,-0.25,20.25);

  tp_et_all_ee=new TH1D("tp_et_all_ee","",41,-0.25,20.25);
  tp_et_all_ee24=new TH1D("tp_et_all_ee24","",41,-0.25,20.25);
  tp_et_matched_ee=new TH1D("tp_et_matched_ee","",41,-0.25,20.25);
  tp_et_matched_ee2=new TH1D("tp_et_matched_ee2","",41,-0.25,20.25);
  tp_et_matched_ee_rh=new TH1D("tp_et_matched_ee_rh","",41,-0.25,20.25);
  tp_et_matched_ee24_rh=new TH1D("tp_et_matched_ee24_rh","",41,-0.25,20.25);
  tp_et_matched_ee_rh2=new TH1D("tp_et_matched_ee_rh2","",41,-0.25,20.25);
  tp_et_unmatched_ee=new TH1D("tp_et_unmatched_ee","",41,-0.25,20.25);

  ctminustp_et_matched_ee=new TH1D("ctminustp_et_matched_ee","",31,-3.1,3.1);
  ctminustp_et_matched_ee2=new TH1D("ctminustp_et_matched_ee2","",31,-3.1,3.1);

  rhminustp_et_matched_ee_rh=new TH1D("rhminustp_et_matched_ee_rh","",31,-3.1,3.1);
  rhminustp_et_matched_ee_rh2=new TH1D("rhminustp_et_matched_ee_rh2","",31,-3.1,3.1);


 
 


  

  ct_et_vs_ieta_all=new TH2D("ct_et_vs_ieta_all","",57,-28.5,28.5,41,-0.25,20.25);
  ct_et_vs_ieta_matched=new TH2D("ct_et_vs_ieta_matched","",57,-28.5,28.5,41,-0.25,20.25);
  ct_et_vs_ieta_matched2=new TH2D("ct_et_vs_ieta_matched2","",57,-28.5,28.5,41,-0.25,20.25);
  ct_et_vs_ieta_unmatched=new TH2D("ct_et_vs_ieta_unmatched","",57,-28.5,28.5,41,-0.25,20.25);

  rh_et_vs_ieta_matched_rh=new TH2D("rh_et_vs_ieta_matched_rh","",57,-28.5,28.5,41,-0.25,20.25);
  rh_et_vs_ieta_matched_rh2=new TH2D("rh_et_vs_ieta_matched_rh2","",57,-28.5,28.5,41,-0.25,20.25);
 

  tp_et_vs_ieta_all=new TH2D("tp_et_vs_ieta_all","",57,-28.5,28.5,41,-0.25,20.25);
  tp_et_vs_ieta_matched=new TH2D("tp_et_vs_ieta_matched","",57,-28.5,28.5,41,-0.25,20.25);
  tp_et_vs_ieta_matched2=new TH2D("tp_et_vs_ieta_matched2","",57,-28.5,28.5,41,-0.25,20.25);
  tp_et_vs_ieta_matched_rh=new TH2D("tp_et_vs_ieta_matched_rh","",57,-28.5,28.5,41,-0.25,20.25);
  tp_et_vs_ieta_matched_rh2=new TH2D("tp_et_vs_ieta_matched_rh2","",57,-28.5,28.5,41,-0.25,20.25);

  tp_et_vs_ieta_unmatched=new TH2D("tp_et_vs_ieta_unmatched","",57,-28.5,28.5,41,-0.25,20.25);



  ctminustp_et_vs_ieta_matched=new TH2D("ctminustp_et_vs_ieta_matched","",57,-28.5,28.5,31,-3.1,3.1);
  ctminustp_et_vs_ieta_matched2=new TH2D("ctminustp_et_vs_ieta_matched2","",57,-28.5,28.5,31,-3.1,3.1);

  rhminustp_et_vs_ieta_matched_rh=new TH2D("rhminustp_et_vs_ieta_matched_rh","",57,-28.5,28.5,31,-3.1,3.1);
  rhminustp_et_vs_ieta_matched_rh2=new TH2D("rhminustp_et_vs_ieta_matched_rh2","",57,-28.5,28.5,31,-3.1,3.1);


  
 


 

 
  tp_et_vs_ct_et_matched=new TH2D("tp_et_vs_ct_et_matched","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_ct_et_matched2=new TH2D("tp_et_vs_ct_et_matched2","",41,-0.25,20.25,41,-0.25,20.25);

  tp_et_vs_ct_et_matched_eb=new TH2D("tp_et_vs_ct_et_matched_eb","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_ct_et_matched_ee=new TH2D("tp_et_vs_ct_et_matched_ee","",41,-0.25,20.25,41,-0.25,20.25);
  
  tp_et_vs_ct_et_matched_eb2=new TH2D("tp_et_vs_ct_et_matched_eb2","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_ct_et_matched_ee2=new TH2D("tp_et_vs_ct_et_matched_ee2","",41,-0.25,20.25,41,-0.25,20.25);


  tp_et_vs_rh_et_matched_rh=new TH2D("tp_et_vs_rh_et_matched_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_rh2=new TH2D("tp_et_vs_rh_et_matched_rh2","",41,-0.25,20.25,41,-0.25,20.25);
 
  tp_et_vs_rh_et_matched_eb_rh=new TH2D("tp_et_vs_rh_et_matched_eb_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_rh=new TH2D("tp_et_vs_rh_et_matched_ee_rh","",41,-0.25,20.25,41,-0.25,20.25);

  tp_et_vs_rh_et_matched_eb_rh2=new TH2D("tp_et_vs_rh_et_matched_eb_rh2","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_rh2=new TH2D("tp_et_vs_rh_et_matched_ee_rh2","",41,-0.25,20.25,41,-0.25,20.25);
 

  tp_et_vs_rh_et_matched_ee_ieta1820_rh=new TH2D("tp_et_vs_rh_et_matched_ee_ieta1820_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_ieta2122_rh=new TH2D("tp_et_vs_rh_et_matched_ee_ieta2122_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_ieta2324_rh=new TH2D("tp_et_vs_rh_et_matched_ee_ieta2324_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_ieta2526_rh=new TH2D("tp_et_vs_rh_et_matched_ee_ieta2526_rh","",41,-0.25,20.25,41,-0.25,20.25);
  tp_et_vs_rh_et_matched_ee_ieta2728_rh=new TH2D("tp_et_vs_rh_et_matched_ee_ieta2728_rh","",41,-0.25,20.25,41,-0.25,20.25);


  eb_digi_01=new TH2D("eb_digi_01","",10,0,10,1000,0,1000);
  ee_digi_01=new TH2D("ee_digi_01","",10,0,10,1000,0,1000);

  eb_digi_05=new TH2D("eb_digi_05","",10,0,10,1000,0,1000);
  ee_digi_05=new TH2D("ee_digi_05","",10,0,10,1000,0,1000);

  eb_digi_30=new TH2D("eb_digi_30","",10,0,10,1000,0,1000);
  ee_digi_30=new TH2D("ee_digi_30","",10,0,10,1000,0,1000);
  

  eb_digi_0105=new TH2D("eb_digi_0105","",10,0,10,100,190,290);
  ee_digi_0105=new TH2D("ee_digi_0105","",10,0,10,100,190,290);

  eb_digi_0530=new TH2D("eb_digi_0530","",10,0,10,200,190,390);
  ee_digi_0530=new TH2D("ee_digi_0530","",10,0,10,200,190,390);


  eb_digi_0105_vs_time=new TH2D("eb_digi_0105_vs_time","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time=new TH2D("ee_digi_0105_vs_time","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time=new TH2D("eb_digi_0530_vs_time","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time=new TH2D("ee_digi_0530_vs_time","",120,-60,60,10,0,10);


  eb_digi_0105_vs_time_norm=new TH2D("eb_digi_0105_vs_time_norm","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_norm=new TH2D("ee_digi_0105_vs_time_norm","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_norm=new TH2D("eb_digi_0530_vs_time_norm","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_norm=new TH2D("ee_digi_0530_vs_time_norm","",120,-60,60,10,0,10);
  

  eb_digi_0105_vs_bxtrain=new TH2D("eb_digi_0105_vs_bxtrain","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain=new TH2D("ee_digi_0105_vs_bxtrain","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain=new TH2D("eb_digi_0530_vs_bxtrain","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain=new TH2D("ee_digi_0530_vs_bxtrain","",40,-2,38,10,0,10);



  eb_digi_0105_vs_bxtrain_norm=new TH2D("eb_digi_0105_vs_bxtrain_norm","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_norm=new TH2D("ee_digi_0105_vs_bxtrain_norm","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_norm=new TH2D("eb_digi_0530_vs_bxtrain_norm","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_norm=new TH2D("ee_digi_0530_vs_bxtrain_norm","",40,-2,38,10,0,10);



  eb_digi_0105_vs_time_eta15=new TH2D("eb_digi_0105_vs_time_eta15","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_eta15=new TH2D("ee_digi_0105_vs_time_eta15","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_eta15=new TH2D("eb_digi_0530_vs_time_eta15","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_eta15=new TH2D("ee_digi_0530_vs_time_eta15","",120,-60,60,10,0,10);


  eb_digi_0105_vs_time_norm_eta15=new TH2D("eb_digi_0105_vs_time_norm_eta15","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_norm_eta15=new TH2D("ee_digi_0105_vs_time_norm_eta15","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_norm_eta15=new TH2D("eb_digi_0530_vs_time_norm_eta15","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_norm_eta15=new TH2D("ee_digi_0530_vs_time_norm_eta15","",120,-60,60,10,0,10);
  

  eb_digi_0105_vs_bxtrain_eta15=new TH2D("eb_digi_0105_vs_bxtrain_eta15","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_eta15=new TH2D("ee_digi_0105_vs_bxtrain_eta15","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_eta15=new TH2D("eb_digi_0530_vs_bxtrain_eta15","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_eta15=new TH2D("ee_digi_0530_vs_bxtrain_eta15","",40,-2,38,10,0,10);



  eb_digi_0105_vs_bxtrain_norm_eta15=new TH2D("eb_digi_0105_vs_bxtrain_norm_eta15","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_norm_eta15=new TH2D("ee_digi_0105_vs_bxtrain_norm_eta15","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_norm_eta15=new TH2D("eb_digi_0530_vs_bxtrain_norm_eta15","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_norm_eta15=new TH2D("ee_digi_0530_vs_bxtrain_norm_eta15","",40,-2,38,10,0,10);



  eb_digi_0105_vs_time_eta20=new TH2D("eb_digi_0105_vs_time_eta20","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_eta20=new TH2D("ee_digi_0105_vs_time_eta20","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_eta20=new TH2D("eb_digi_0530_vs_time_eta20","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_eta20=new TH2D("ee_digi_0530_vs_time_eta20","",120,-60,60,10,0,10);


  eb_digi_0105_vs_time_norm_eta20=new TH2D("eb_digi_0105_vs_time_norm_eta20","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_norm_eta20=new TH2D("ee_digi_0105_vs_time_norm_eta20","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_norm_eta20=new TH2D("eb_digi_0530_vs_time_norm_eta20","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_norm_eta20=new TH2D("ee_digi_0530_vs_time_norm_eta20","",120,-60,60,10,0,10);
  

  eb_digi_0105_vs_bxtrain_eta20=new TH2D("eb_digi_0105_vs_bxtrain_eta20","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_eta20=new TH2D("ee_digi_0105_vs_bxtrain_eta20","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_eta20=new TH2D("eb_digi_0530_vs_bxtrain_eta20","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_eta20=new TH2D("ee_digi_0530_vs_bxtrain_eta20","",40,-2,38,10,0,10);



  eb_digi_0105_vs_bxtrain_norm_eta20=new TH2D("eb_digi_0105_vs_bxtrain_norm_eta20","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_norm_eta20=new TH2D("ee_digi_0105_vs_bxtrain_norm_eta20","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_norm_eta20=new TH2D("eb_digi_0530_vs_bxtrain_norm_eta20","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_norm_eta20=new TH2D("ee_digi_0530_vs_bxtrain_norm_eta20","",40,-2,38,10,0,10);




  eb_digi_0105_vs_time_eta25=new TH2D("eb_digi_0105_vs_time_eta25","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_eta25=new TH2D("ee_digi_0105_vs_time_eta25","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_eta25=new TH2D("eb_digi_0530_vs_time_eta25","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_eta25=new TH2D("ee_digi_0530_vs_time_eta25","",120,-60,60,10,0,10);


  eb_digi_0105_vs_time_norm_eta25=new TH2D("eb_digi_0105_vs_time_norm_eta25","",120,-60,60,10,0,10);
  ee_digi_0105_vs_time_norm_eta25=new TH2D("ee_digi_0105_vs_time_norm_eta25","",120,-60,60,10,0,10);

  eb_digi_0530_vs_time_norm_eta25=new TH2D("eb_digi_0530_vs_time_norm_eta25","",120,-60,60,10,0,10);
  ee_digi_0530_vs_time_norm_eta25=new TH2D("ee_digi_0530_vs_time_norm_eta25","",120,-60,60,10,0,10);
  

  eb_digi_0105_vs_bxtrain_eta25=new TH2D("eb_digi_0105_vs_bxtrain_eta25","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_eta25=new TH2D("ee_digi_0105_vs_bxtrain_eta25","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_eta25=new TH2D("eb_digi_0530_vs_bxtrain_eta25","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_eta25=new TH2D("ee_digi_0530_vs_bxtrain_eta25","",40,-2,38,10,0,10);



  eb_digi_0105_vs_bxtrain_norm_eta25=new TH2D("eb_digi_0105_vs_bxtrain_norm_eta25","",40,-2,38,10,0,10);
  ee_digi_0105_vs_bxtrain_norm_eta25=new TH2D("ee_digi_0105_vs_bxtrain_norm_eta25","",40,-2,38,10,0,10);

  eb_digi_0530_vs_bxtrain_norm_eta25=new TH2D("eb_digi_0530_vs_bxtrain_norm_eta25","",40,-2,38,10,0,10);
  ee_digi_0530_vs_bxtrain_norm_eta25=new TH2D("ee_digi_0530_vs_bxtrain_norm_eta25","",40,-2,38,10,0,10);


  //modif-alex
  h_ttonline_et  = new TH1F( "h_ttonline_et", "h_ttonline_et", 256,0,256);
  h_ttoffline_et = new TH1F( "h_ttoffline_et", "h_ttoffline_et", 256,0,256);


  dataTree_      = new TTree("dataTree","dataTree");

  // physics declared and technical trigger bits

  dataTree_->Branch("physdeclared",&physdeclared,"physdeclared/I");  
  dataTree_->Branch("bit36",       &bit36,   	 "bit36/I");
  dataTree_->Branch("bit37",       &bit37,   	 "bit37/I");
  dataTree_->Branch("bit38",       &bit38,   	 "bit38/I");
  dataTree_->Branch("bit39",       &bit39,   	 "bit39/I");
  dataTree_->Branch("bit40",       &bit40,   	 "bit40/I");
  dataTree_->Branch("bit41",       &bit41,   	 "bit41/I");
  dataTree_->Branch("bit3",        &bit3,        "bit3/I");
  dataTree_->Branch("bit4",        &bit4,        "bit4/I");
  dataTree_->Branch("bit9",        &bit9,        "bit9/I");
  
  dataTree_->Branch("bit0",        &bit0,        "bit0/I");


  dataTree_->Branch("eg1",         &eg1,         "eg1/I");
  dataTree_->Branch("eg2",         &eg2,         "eg2/I");
  dataTree_->Branch("eg5",         &eg5,         "eg5/I");
  dataTree_->Branch("algo124",     &algo124,     "algo124/I");


  // muon triggers

  dataTree_->Branch("algo54",     &algo54,     "algo54/I");
  dataTree_->Branch("algo55",     &algo55,     "algo55/I");
  dataTree_->Branch("algo56",     &algo56,     "algo56/I");
  dataTree_->Branch("algo57",     &algo57,     "algo57/I");
  dataTree_->Branch("algo58",     &algo58,     "algo58/I");
  dataTree_->Branch("algo59",     &algo59,     "algo59/I");
  dataTree_->Branch("algo60",     &algo60,     "algo60/I");
  dataTree_->Branch("algo61",     &algo61,     "algo61/I");
  dataTree_->Branch("algo62",     &algo62,     "algo62/I");

  dataTree_->Branch("algo106",    &algo106,    "algo106/I");
  dataTree_->Branch("algo107",    &algo107,    "algo107/I");


  // run/event info
  
  dataTree_->Branch("run",         &run,     	 "run/I");  
  dataTree_->Branch("even",        &even,   	 "even/I");
  dataTree_->Branch("lumi",        &lumi,    	 "lumi/I");
  dataTree_->Branch("bx",          &bx,      	 "bx/I");   
  dataTree_->Branch("orbit",       &orbit,   	 "orbit/I");
  dataTree_->Branch("time",        &time,    	 "time/F");
  
 
  // tracks - for monster event cut

  dataTree_->Branch("ntrk",        &ntrk,	  "ntrk/I");
  dataTree_->Branch("goodtrk",     &goodtrk,	  "goodtrk/I");


  //  primary vertex info

  dataTree_->Branch("numvtx",      &numvtx,	  "numvtx/I");
  dataTree_->Branch("numgoodvtx",  &numgoodvtx,	  "numgoodvtx/I");

  dataTree_->Branch("vtx_x",       &vtx_x,	  "vtx_x/F");
  dataTree_->Branch("vtx_y",       &vtx_y,	  "vtx_y/F");
  dataTree_->Branch("vtx_z",       &vtx_z,	  "vtx_z/F");

  dataTree_->Branch("vtx_x_err",   &vtx_x_err,	  "vtx_x_err/F");
  dataTree_->Branch("vtx_y_err",   &vtx_y_err,	  "vtx_y_err/F");
  dataTree_->Branch("vtx_z_err",   &vtx_z_err,	  "vtx_z_err/F");

  dataTree_->Branch("vtx_chi2",    &vtx_chi2,	  "vtx_chi2/F");
  dataTree_->Branch("vtx_ndof",    &vtx_ndof,	  "vtx_ndof/F");

  dataTree_->Branch("vtx_ntracks", &vtx_ntracks,  "vtx_ntracks/I");
  dataTree_->Branch("vtx_isfake",  &vtx_isfake,   "vtx_isfake/I");

  dataTree_->Branch("vtx_good",    &vtx_good,     "vtx_good/I");


  // PF jets

//   dataTree_->Branch("rank",        &rank_,	  "rank_/I");
//   dataTree_->Branch("ncr",         &ncr_,         "ncr_/I");
//   dataTree_->Branch("ptJet",       &ptJet_,       "ptJet_/F");  
//   dataTree_->Branch("etaJet",      &etaJet_,      "etaJet_/F");  
//   dataTree_->Branch("phiJet",      &phiJet_,      "phiJet_/F");
  
//   dataTree_->Branch("nrjets",      &nrjets_,      "nrjets_/I"); 
//   dataTree_->Branch("chfJet",      &chfJet_,      "chfJet_/F");
//   dataTree_->Branch("nhfJet",      &nhfJet_,      "nhfJet_/F");  
//   dataTree_->Branch("cemfJet",     &cemfJet_,     "cemfJet_/F");  
//   dataTree_->Branch("nemfJet",     &nemfJet_,     "nemfJet_/F");  
//   dataTree_->Branch("cmultiJet",   &cmultiJet_,   "cmultiJet_/I"); 
//   dataTree_->Branch("nmultiJet",   &nmultiJet_,   "nmultiJet_/I");
  
//   dataTree_->Branch("energy_pf",   &energy_pf_,   "energy_pf_/F");
//   dataTree_->Branch("energyc_pf",  &energyc_pf_,  "energyc_pf_/F");
//   dataTree_->Branch("energyn_pf",  &energyn_pf_,  "energyn_pf_/F");
//   dataTree_->Branch("energyg_pf",  &energyg_pf_,  "energyg_pf_/F");

//   dataTree_->Branch("energy_ecal",  &energy_ecal,  "energy_ecal/F");
//   dataTree_->Branch("energy_hcal",  &energy_hcal,  "energy_hcal/F");


  dataTree_->Branch("scale",        &scale,        "scale/F");

  dataTree_->Branch("ebmax",        &ebmax,        "ebmax/F");
  dataTree_->Branch("ebmaxet",      &ebmaxet,      "ebmaxet/F");
  dataTree_->Branch("ebtime",       &ebtime,       "ebtime/F");
  dataTree_->Branch("ebflags",      &ebflags,      "ebflags/I");
  dataTree_->Branch("eb_ieta",      &eb_ieta,      "eb_ieta/I");
  dataTree_->Branch("eb_iphi",      &eb_iphi,      "eb_iphi/I");  
  dataTree_->Branch("eb_eta",       &eb_eta,       "eb_eta/F");
  dataTree_->Branch("eb_phi",       &eb_phi,       "eb_phi/F");
  dataTree_->Branch("ebhits",       &ebhits,       "ebhits/I");
  dataTree_->Branch("ebhits1GeV",   &ebhits1GeV,   "ebhits1GeV/I");
  dataTree_->Branch("ebhits2GeV",   &ebhits2GeV,   "ebhits2GeV/I");
  dataTree_->Branch("ebhits4GeV",   &ebhits4GeV,   "ebhits4GeV/I");
  dataTree_->Branch("ebhits1GeVet",   &ebhits1GeVet,   "ebhits1GeVet/I");
  dataTree_->Branch("ebhits2GeVet",   &ebhits2GeVet,   "ebhits2GeVet/I");
  dataTree_->Branch("ebhits4GeVet",   &ebhits4GeVet,   "ebhits4GeVet/I");
  dataTree_->Branch("eb_r9",        &eb_r9,        "eb_r9/F");
  dataTree_->Branch("eb_r4",        &eb_r4,        "eb_r4/F");


  dataTree_->Branch("eb_e9",        &eb_e9,        "eb_e9/F");
  dataTree_->Branch("eb_e25",       &eb_e25,       "eb_e25/F");


  dataTree_->Branch("ebmax2",        &ebmax2,        "ebmax2/F");
  dataTree_->Branch("ebmaxet2",      &ebmaxet2,      "ebmaxet2/F");
  dataTree_->Branch("ebtime2",       &ebtime2,       "ebtime2/F");
  dataTree_->Branch("ebflags2",      &ebflags2,      "ebflags2/I");
  dataTree_->Branch("eb_ieta2",      &eb_ieta2,      "eb_ieta2/I");
  dataTree_->Branch("eb_iphi2",      &eb_iphi2,      "eb_iphi2/I");  
  dataTree_->Branch("eb_eta2",       &eb_eta2,       "eb_eta2/F");
  dataTree_->Branch("eb_phi2",       &eb_phi2,       "eb_phi2/F");
  dataTree_->Branch("ebhits1GeV2",   &ebhits1GeV2,   "ebhits1GeV2/I");
  dataTree_->Branch("eb_r92",        &eb_r92,        "eb_r92/F");
  dataTree_->Branch("eb_r42",        &eb_r42,        "eb_r42/F");


  dataTree_->Branch("ebchi2",        &ebchi2,        "ebchi2/F");
  dataTree_->Branch("ebchi2oot",     &ebchi2oot,     "ebchi2oot/F");
  dataTree_->Branch("eb2chi2",       &eb2chi2,       "eb2chi2/F");
  dataTree_->Branch("eb2chi2oot",    &eb2chi2oot,    "eb2chi2oot/F");


  dataTree_->Branch("ebsum_gt1",    &ebsum_gt1,    "ebsum_gt1/F");
  dataTree_->Branch("ebsum_gt2",    &ebsum_gt2,    "ebsum_gt2/F");
  dataTree_->Branch("ebsum_gt4",    &ebsum_gt4,    "ebsum_gt4/F");


  dataTree_->Branch("ebsum_gt1et",    &ebsum_gt1et,    "ebsum_gt1et/F");
  dataTree_->Branch("ebsum_gt2et",    &ebsum_gt2et,    "ebsum_gt2et/F");
  dataTree_->Branch("ebsum_gt4et",    &ebsum_gt4et,    "ebsum_gt4et/F");



  dataTree_->Branch("eesum_gt1",    &eesum_gt1,    "eesum_gt1/F");
  dataTree_->Branch("eesum_gt2",    &eesum_gt2,    "eesum_gt2/F");
  dataTree_->Branch("eesum_gt4",    &eesum_gt4,    "eesum_gt4/F");


  dataTree_->Branch("eesum_gt1et",    &eesum_gt1et,    "eesum_gt1et/F");
  dataTree_->Branch("eesum_gt2et",    &eesum_gt2et,    "eesum_gt2et/F");
  dataTree_->Branch("eesum_gt4et",    &eesum_gt4et,    "eesum_gt4et/F");


  dataTree_->Branch("ebflag_kgood",      &ebflag_kgood,       "ebflag_kgood/I");
  dataTree_->Branch("ebflag_kpoorreco", &ebflag_kpoorreco,  "ebflag_kpoorreco/I");
  dataTree_->Branch("ebflag_koutoftime", &ebflag_koutoftime,  "ebflag_koutoftime/I");
  dataTree_->Branch("ebflag_kfake",      &ebflag_kfake,       "ebflag_kfake/I");


  dataTree_->Branch("eemax",        &eemax,        "eemax/F");
  dataTree_->Branch("eemaxet",      &eemaxet,      "eemaxet/F");
  dataTree_->Branch("eetime",       &eetime,       "eetime/F");
  dataTree_->Branch("eeflags",      &eeflags,      "eeflags/I");
  dataTree_->Branch("eeix",         &eeix,         "eeix/I");
  dataTree_->Branch("eeiy",         &eeiy,         "eeiy/I");
  dataTree_->Branch("eeiz",         &eeiz,         "eeiz/I");
  dataTree_->Branch("ee_eta",       &ee_eta,       "ee_eta/F");
  dataTree_->Branch("ee_phi",       &ee_phi,       "ee_phi/F");
  dataTree_->Branch("eehits",       &eehits,       "eehits/I");
  dataTree_->Branch("eehits1GeV",   &eehits1GeV,   "eehits1GeV/I");
  dataTree_->Branch("eehits2GeV",   &eehits2GeV,   "eehits2GeV/I");
  dataTree_->Branch("eehits4GeV",   &eehits4GeV,   "eehits4GeV/I");
  dataTree_->Branch("eehits1GeVet",   &eehits1GeVet,   "eehits1GeVet/I");
  dataTree_->Branch("eehits2GeVet",   &eehits2GeVet,   "eehits2GeVet/I");
  dataTree_->Branch("eehits4GeVet",   &eehits4GeVet,   "eehits4GeVet/I");
  dataTree_->Branch("ee_r9",        &ee_r9,        "ee_r9/F");
  dataTree_->Branch("eephits",      &eephits,      "eephits/I");
  dataTree_->Branch("eemhits",      &eemhits,      "eemhits/I");


  dataTree_->Branch("eemax2",        &eemax2,        "eemax2/F");
  dataTree_->Branch("eemaxet2",      &eemaxet2,      "eemaxet2/F");
  dataTree_->Branch("eetime2",       &eetime2,       "eetime2/F");
  dataTree_->Branch("eeflags2",      &eeflags2,      "eeflags2/I");
  dataTree_->Branch("eeix2",         &eeix2,         "eeix2/I");
  dataTree_->Branch("eeiy2",         &eeiy2,         "eeiy2/I");
  dataTree_->Branch("eeiz2",         &eeiz2,         "eeiz2/I");
  dataTree_->Branch("ee_eta2",       &ee_eta2,       "ee_eta2/F");
  dataTree_->Branch("ee_phi2",       &ee_phi2,       "ee_phi2/F");
  dataTree_->Branch("eehits1GeV2",   &eehits1GeV2,   "eehits1GeV2/I");
  dataTree_->Branch("ee_r92",        &ee_r92,        "ee_r92/F");


  dataTree_->Branch("tmean_en",      &tmean_en,      "tmean_en/F");
  dataTree_->Branch("terr_en",       &terr_en,       "terr_en/F");
  dataTree_->Branch("tmean_sig",     &tmean_sig,     "tmean_sig/F");
  dataTree_->Branch("terr_sig",      &terr_sig,      "terr_sig/F");


  dataTree_->Branch("r4count",      &r4count,      "r4count/I");

  dataTree_->Branch("e2e9count_thresh0",  &e2e9count_thresh0,  "e2e9count_thresh0/I");
  dataTree_->Branch("e2e25count_thresh0",  &e2e25count_thresh0,  "e2e25count_thresh0/I");
  dataTree_->Branch("e2e9count_thresh1",  &e2e9count_thresh1,  "e2e9count_thresh1/I");
  dataTree_->Branch("e2e25count_thresh1",  &e2e25count_thresh1,  "e2e25count_thresh1/I");

  dataTree_->Branch("e2e9count_thresh0_nor4",  &e2e9count_thresh0_nor4,  "e2e9count_thresh0_nor4/I");
  dataTree_->Branch("e2e25count_thresh0_nor4",  &e2e25count_thresh0_nor4,  "e2e25count_thresh0_nor4/I");
  dataTree_->Branch("e2e9count_thresh1_nor4",  &e2e9count_thresh1_nor4,  "e2e9count_thresh1_nor4/I");
  dataTree_->Branch("e2e25count_thresh1_nor4",  &e2e25count_thresh1_nor4,  "e2e25count_thresh1_nor4/I");

  dataTree_->Branch("r4_algo_count",  &r4_algo_count,  "r4_algo_count/I");

  dataTree_->Branch("e2e9_algo_count",  &e2e9_algo_count,  "e2e9_algo_count/I");
  dataTree_->Branch("e2e9_algo_count_5_1",  &e2e9_algo_count_5_1,  "e2e9_algo_count_5_1/I");
  dataTree_->Branch("e2e9_algo_count_5_0",  &e2e9_algo_count_5_0,  "e2e9_algo_count_5_0/I");
 
  dataTree_->Branch("swisscross_algo",  &swisscross_algo,  "swisscross_algo/F");  dataTree_->Branch("e2e9_algo",  &e2e9_algo,  "e2e9_algo/F");
 


  // calotowers


  dataTree_->Branch("ncalotower", &ncalotower, "ncalotower/I");
 

  dataTree_->Branch("ncalotowereb", &ncalotowereb,"ncalotowereb/I");
  dataTree_->Branch("ncalotoweree", &ncalotoweree,"ncalotoweree/I");
  dataTree_->Branch("ncalotowerhf", &ncalotowerhf,"ncalotowerhf/I");



  dataTree_->Branch("ncalotowerebgt1",  &ncalotowerebgt1,  "ncalotowerebgt1/I");
  dataTree_->Branch("ncalotowerebgt2",  &ncalotowerebgt2,  "ncalotowerebgt2/I");
  dataTree_->Branch("ncalotowerebgt5",  &ncalotowerebgt5,  "ncalotowerebgt5/I");
  dataTree_->Branch("ncalotowerebgt10", &ncalotowerebgt10, "ncalotowerebgt10/I");

  dataTree_->Branch("ncalotowereegt1",  &ncalotowereegt1,  "ncalotowereegt1/I");
  dataTree_->Branch("ncalotowereegt2",  &ncalotowereegt2,  "ncalotowereegt2/I");
  dataTree_->Branch("ncalotowereegt5",  &ncalotowereegt5,  "ncalotowereegt5/I");
  dataTree_->Branch("ncalotowereegt10", &ncalotowereegt10, "ncalotowereegt10/I");

  dataTree_->Branch("ncalotowerhfgt1",  &ncalotowerhfgt1,  "ncalotowerhfgt1/I");
  dataTree_->Branch("ncalotowerhfgt2",  &ncalotowerhfgt2,  "ncalotowerhfgt2/I");
  dataTree_->Branch("ncalotowerhfgt5",  &ncalotowerhfgt5,  "ncalotowerhfgt5/I");
  dataTree_->Branch("ncalotowerhfgt10", &ncalotowerhfgt10, "ncalotowerhfgt10/I");



  dataTree_->Branch("ncalotowerebgt1had",  &ncalotowerebgt1had,  "ncalotowerebgt1had/I");
  dataTree_->Branch("ncalotowerebgt2had",  &ncalotowerebgt2had,  "ncalotowerebgt2had/I");
  dataTree_->Branch("ncalotowerebgt5had",  &ncalotowerebgt5had,  "ncalotowerebgt5had/I");
  dataTree_->Branch("ncalotowerebgt10had", &ncalotowerebgt10had, "ncalotowerebgt10had/I");

  dataTree_->Branch("ncalotowereegt1had",  &ncalotowereegt1had,  "ncalotowereegt1had/I");
  dataTree_->Branch("ncalotowereegt2had",  &ncalotowereegt2had,  "ncalotowereegt2had/I");
  dataTree_->Branch("ncalotowereegt5had",  &ncalotowereegt5had,  "ncalotowereegt5had/I");
  dataTree_->Branch("ncalotowereegt10had", &ncalotowereegt10had, "ncalotowereegt10had/I");

  dataTree_->Branch("ncalotowerhfgt1had",  &ncalotowerhfgt1had,  "ncalotowerhfgt1had/I");
  dataTree_->Branch("ncalotowerhfgt2had",  &ncalotowerhfgt2had,  "ncalotowerhfgt2had/I");
  dataTree_->Branch("ncalotowerhfgt5had",  &ncalotowerhfgt5had,  "ncalotowerhfgt5had/I");
  dataTree_->Branch("ncalotowerhfgt10had", &ncalotowerhfgt10had, "ncalotowerhfgt10had/I");



  dataTree_->Branch("ncalotowerebgt1hof",  &ncalotowerebgt1hof,  "ncalotowerebgt1hof/I");
  dataTree_->Branch("ncalotowerebgt2hof",  &ncalotowerebgt2hof,  "ncalotowerebgt2hof/I");
  dataTree_->Branch("ncalotowerebgt5hof",  &ncalotowerebgt5hof,  "ncalotowerebgt5hof/I");
  dataTree_->Branch("ncalotowerebgt10hof", &ncalotowerebgt10hof, "ncalotowerebgt10hof/I");

  dataTree_->Branch("ncalotowereegt1hof",  &ncalotowereegt1hof,  "ncalotowereegt1hof/I");
  dataTree_->Branch("ncalotowereegt2hof",  &ncalotowereegt2hof,  "ncalotowereegt2hof/I");
  dataTree_->Branch("ncalotowereegt5hof",  &ncalotowereegt5hof,  "ncalotowereegt5hof/I");
  dataTree_->Branch("ncalotowereegt10hof", &ncalotowereegt10hof, "ncalotowereegt10hof/I");

  dataTree_->Branch("ncalotowerhfgt1hof",  &ncalotowerhfgt1hof,  "ncalotowerhfgt1hof/I");
  dataTree_->Branch("ncalotowerhfgt2hof",  &ncalotowerhfgt2hof,  "ncalotowerhfgt2hof/I");
  dataTree_->Branch("ncalotowerhfgt5hof",  &ncalotowerhfgt5hof,  "ncalotowerhfgt5hof/I");
  dataTree_->Branch("ncalotowerhfgt10hof", &ncalotowerhfgt10hof, "ncalotowerhfgt10hof/I");


  dataTree_->Branch("ncalotowerebgt1all",  &ncalotowerebgt1all,  "ncalotowerebgt1all/I");
  dataTree_->Branch("ncalotowerebgt2all",  &ncalotowerebgt2all,  "ncalotowerebgt2all/I");
  dataTree_->Branch("ncalotowerebgt5all",  &ncalotowerebgt5all,  "ncalotowerebgt5all/I");
  dataTree_->Branch("ncalotowerebgt10all", &ncalotowerebgt10all, "ncalotowerebgt10all/I");

  dataTree_->Branch("ncalotowereegt1all",  &ncalotowereegt1all,  "ncalotowereegt1all/I");
  dataTree_->Branch("ncalotowereegt2all",  &ncalotowereegt2all,  "ncalotowereegt2all/I");
  dataTree_->Branch("ncalotowereegt5all",  &ncalotowereegt5all,  "ncalotowereegt5all/I");
  dataTree_->Branch("ncalotowereegt10all", &ncalotowereegt10all, "ncalotowereegt10all/I");

  dataTree_->Branch("ncalotowerhfgt1all",  &ncalotowerhfgt1all,  "ncalotowerhfgt1all/I");
  dataTree_->Branch("ncalotowerhfgt2all",  &ncalotowerhfgt2all,  "ncalotowerhfgt2all/I");
  dataTree_->Branch("ncalotowerhfgt5all",  &ncalotowerhfgt5all,  "ncalotowerhfgt5all/I");
  dataTree_->Branch("ncalotowerhfgt10all", &ncalotowerhfgt10all, "ncalotowerhfgt10all/I");



  dataTree_->Branch("ncalotowerebgt05",    &ncalotowerebgt05,    "ncalotowerebgt05/I");
  dataTree_->Branch("ncalotowereegt05",    &ncalotowereegt05,    "ncalotowereegt05/I");
  dataTree_->Branch("ncalotowerhfgt05",    &ncalotowerhfgt05,    "ncalotowerhfgt05/I");


  dataTree_->Branch("ncalotowerebgt05had", &ncalotowerebgt05had, "ncalotowerebgt05had/I");
  dataTree_->Branch("ncalotowereegt05had", &ncalotowereegt05had, "ncalotowereegt05had/I");
  dataTree_->Branch("ncalotowerhfgt05had", &ncalotowerhfgt05had, "ncalotowerhfgt05had/I");

  dataTree_->Branch("ncalotowerebgt05hof", &ncalotowerebgt05hof, "ncalotowerebgt05hof/I");
  dataTree_->Branch("ncalotowereegt05hof", &ncalotowereegt05hof, "ncalotowereegt05hof/I");
  dataTree_->Branch("ncalotowerhfgt05hof", &ncalotowerhfgt05hof, "ncalotowerhfgt05hof/I");

  dataTree_->Branch("ncalotowerebgt05all", &ncalotowerebgt05all, "ncalotowerebgt05all/I");
  dataTree_->Branch("ncalotowereegt05all", &ncalotowereegt05all, "ncalotowereegt05all/I");
  dataTree_->Branch("ncalotowerhfgt05all", &ncalotowerhfgt05all, "ncalotowerhfgt05all/I");



  dataTree_->Branch("ctsumebgt1",  &ctsumebgt1,  "ctsumebgt1/F");
  dataTree_->Branch("ctsumebgt2",  &ctsumebgt2,  "ctsumebgt2/F");
  dataTree_->Branch("ctsumebgt5",  &ctsumebgt5,  "ctsumebgt5/F");
  dataTree_->Branch("ctsumebgt10", &ctsumebgt10, "ctsumebgt10/F");


  dataTree_->Branch("ctsumeegt1",  &ctsumeegt1,  "ctsumeegt1/F");
  dataTree_->Branch("ctsumeegt2",  &ctsumeegt2,  "ctsumeegt2/F");
  dataTree_->Branch("ctsumeegt5",  &ctsumeegt5,  "ctsumeegt5/F");
  dataTree_->Branch("ctsumeegt10", &ctsumeegt10, "ctsumeegt10/F");


  dataTree_->Branch("ctsumhfgt1",  &ctsumhfgt1,  "ctsumhfgt1/F");
  dataTree_->Branch("ctsumhfgt2",  &ctsumhfgt2,  "ctsumhfgt2/F");
  dataTree_->Branch("ctsumhfgt5",  &ctsumhfgt5,  "ctsumhfgt5/F");
  dataTree_->Branch("ctsumhfgt10", &ctsumhfgt10, "ctsumhfgt10/F");



  dataTree_->Branch("ctsumebgt05",       &ctsumebgt05,       "ctsumebgt05/F");
  dataTree_->Branch("ctsumeegt05",       &ctsumeegt05,       "ctsumeegt05/F");
  dataTree_->Branch("ctsumhfgt05",       &ctsumhfgt05,       "ctsumhfgt05/F");

  dataTree_->Branch("ctsumebgt05had",    &ctsumebgt05had,    "ctsumebgt05had/F");
  dataTree_->Branch("ctsumeegt05had",    &ctsumeegt05had,    "ctsumeegt05had/F");
  dataTree_->Branch("ctsumhfgt05had",    &ctsumhfgt05had,    "ctsumhfgt05had/F");

  dataTree_->Branch("ctsumebgt05hof",    &ctsumebgt05hof,    "ctsumebgt05hof/F");
  dataTree_->Branch("ctsumeegt05hof",    &ctsumeegt05hof,    "ctsumeegt05hof/F");
  dataTree_->Branch("ctsumhfgt05hof",    &ctsumhfgt05hof,    "ctsumhfgt05hof/F");

  dataTree_->Branch("ctsumebgt05all",    &ctsumebgt05all,    "ctsumebgt05all/F");
  dataTree_->Branch("ctsumeegt05all",    &ctsumeegt05all,    "ctsumeegt05all/F");
  dataTree_->Branch("ctsumhfgt05all",    &ctsumhfgt05all,    "ctsumhfgt05all/F");



  dataTree_->Branch("ctsumebgt1had",  &ctsumebgt1had,  "ctsumebgt1had/F");
  dataTree_->Branch("ctsumebgt2had",  &ctsumebgt2had,  "ctsumebgt2had/F");
  dataTree_->Branch("ctsumebgt5had",  &ctsumebgt5had,  "ctsumebgt5had/F");
  dataTree_->Branch("ctsumebgt10had", &ctsumebgt10had, "ctsumebgt10had/F");


  dataTree_->Branch("ctsumeegt1had",  &ctsumeegt1had,  "ctsumeegt1had/F");
  dataTree_->Branch("ctsumeegt2had",  &ctsumeegt2had,  "ctsumeegt2had/F");
  dataTree_->Branch("ctsumeegt5had",  &ctsumeegt5had,  "ctsumeegt5had/F");
  dataTree_->Branch("ctsumeegt10had", &ctsumeegt10had, "ctsumeegt10had/F");


  dataTree_->Branch("ctsumhfgt1had",  &ctsumhfgt1had,  "ctsumhfgt1had/F");
  dataTree_->Branch("ctsumhfgt2had",  &ctsumhfgt2had,  "ctsumhfgt2had/F");
  dataTree_->Branch("ctsumhfgt5had",  &ctsumhfgt5had,  "ctsumhfgt5had/F");
  dataTree_->Branch("ctsumhfgt10had", &ctsumhfgt10had, "ctsumhfgt10had/F");


  dataTree_->Branch("ctsumebgt1hof",  &ctsumebgt1hof,  "ctsumebgt1hof/F");
  dataTree_->Branch("ctsumebgt2hof",  &ctsumebgt2hof,  "ctsumebgt2hof/F");
  dataTree_->Branch("ctsumebgt5hof",  &ctsumebgt5hof,  "ctsumebgt5hof/F");
  dataTree_->Branch("ctsumebgt10hof", &ctsumebgt10hof, "ctsumebgt10hof/F");


  dataTree_->Branch("ctsumeegt1hof",  &ctsumeegt1hof,  "ctsumeegt1hof/F");
  dataTree_->Branch("ctsumeegt2hof",  &ctsumeegt2hof,  "ctsumeegt2hof/F");
  dataTree_->Branch("ctsumeegt5hof",  &ctsumeegt5hof,  "ctsumeegt5hof/F");
  dataTree_->Branch("ctsumeegt10hof", &ctsumeegt10hof, "ctsumeegt10hof/F");


  dataTree_->Branch("ctsumhfgt1hof",  &ctsumhfgt1hof,  "ctsumhfgt1hof/F");
  dataTree_->Branch("ctsumhfgt2hof",  &ctsumhfgt2hof,  "ctsumhfgt2hof/F");
  dataTree_->Branch("ctsumhfgt5hof",  &ctsumhfgt5hof,  "ctsumhfgt5hof/F");
  dataTree_->Branch("ctsumhfgt10hof", &ctsumhfgt10hof, "ctsumhfgt10hof/F");


  dataTree_->Branch("ctsumebgt1all",  &ctsumebgt1all,  "ctsumebgt1all/F");
  dataTree_->Branch("ctsumebgt2all",  &ctsumebgt2all,  "ctsumebgt2all/F");
  dataTree_->Branch("ctsumebgt5all",  &ctsumebgt5all,  "ctsumebgt5all/F");
  dataTree_->Branch("ctsumebgt10all", &ctsumebgt10all, "ctsumebgt10all/F");


  dataTree_->Branch("ctsumeegt1all",  &ctsumeegt1all,  "ctsumeegt1all/F");
  dataTree_->Branch("ctsumeegt2all",  &ctsumeegt2all,  "ctsumeegt2all/F");
  dataTree_->Branch("ctsumeegt5all",  &ctsumeegt5all,  "ctsumeegt5all/F");
  dataTree_->Branch("ctsumeegt10all", &ctsumeegt10all, "ctsumeegt10all/F");


  dataTree_->Branch("ctsumhfgt1all",  &ctsumhfgt1all,  "ctsumhfgt1all/F");
  dataTree_->Branch("ctsumhfgt2all",  &ctsumhfgt2all,  "ctsumhfgt2all/F");
  dataTree_->Branch("ctsumhfgt5all",  &ctsumhfgt5all,  "ctsumhfgt5all/F");
  dataTree_->Branch("ctsumhfgt10all", &ctsumhfgt10all, "ctsumhfgt10all/F");




  dataTree_->Branch("rechitsumet_eb_all", &rechitsumet_eb_all, "rechitsumet_eb_all/F");
  dataTree_->Branch("rechitsumet_eb_01", &rechitsumet_eb_01, "rechitsumet_eb_01/F");
  dataTree_->Branch("rechitsumet_eb_05", &rechitsumet_eb_05, "rechitsumet_eb_05/F");

  dataTree_->Branch("rechitsumet_eb_0105", &rechitsumet_eb_0105, "rechitsumet_eb_0105/F");
  dataTree_->Branch("rechitsumet_eb_0530", &rechitsumet_eb_0530, "rechitsumet_eb_0530/F");

  dataTree_->Branch("rechitsumet_ee_all", &rechitsumet_ee_all, "rechitsumet_ee_all/F");
  dataTree_->Branch("rechitsumet_ee_01", &rechitsumet_ee_01, "rechitsumet_ee_01/F");
  dataTree_->Branch("rechitsumet_ee_05", &rechitsumet_ee_05, "rechitsumet_ee_05/F");

  dataTree_->Branch("rechitsumet_ee_0105", &rechitsumet_ee_0105, "rechitsumet_ee_0105/F");
  dataTree_->Branch("rechitsumet_ee_0530", &rechitsumet_ee_0530, "rechitsumet_ee_0530/F");



  dataTree_->Branch("bunchintrain", &bunchintrain, "bunchintrain/I");


  dataTree_->Branch("ebnumsc_all", &ebnumsc_all, "ebnumsc_all/I");
  dataTree_->Branch("eenumsc_all", &eenumsc_all, "eenumsc_all/I");


  dataTree_->Branch("ebscsumet_all", &ebscsumet_all, "ebscsumet_all/F");
  dataTree_->Branch("eescsumet_all", &eescsumet_all, "eescsumet_all/F");
  dataTree_->Branch("ebscsumet_all_eta15", &ebscsumet_all_eta15, "ebscsumet_all_eta15/F");
  dataTree_->Branch("ebscsumet_all_eta20", &ebscsumet_all_eta20, "ebscsumet_all_eta20/F");
  dataTree_->Branch("ebscsumet_all_eta25", &ebscsumet_all_eta25, "ebscsumet_all_eta25/F");
  dataTree_->Branch("eescsumet_all_eta15", &eescsumet_all_eta15, "eescsumet_all_eta15/F");
  dataTree_->Branch("eescsumet_all_eta20", &eescsumet_all_eta20, "eescsumet_all_eta20/F");
  dataTree_->Branch("eescsumet_all_eta25", &eescsumet_all_eta25, "eescsumet_all_eta25/F");

  dataTree_->Branch("ebnumrechits_01", &ebnumrechits_01, "ebnumrechits_01/I");
  dataTree_->Branch("ebnumrechits_0105", &ebnumrechits_0105, "ebnumrechits_0105/I");
  dataTree_->Branch("ebnumrechits_05", &ebnumrechits_05, "ebnumrechits_05/I");
  dataTree_->Branch("ebnumrechits_0530", &ebnumrechits_0530, "ebnumrechits_0530/I");


  dataTree_->Branch("eenumrechits_01", &eenumrechits_01, "eenumrechits_01/I");
  dataTree_->Branch("eenumrechits_0105", &eenumrechits_0105, "eenumrechits_0105/I");
  dataTree_->Branch("eenumrechits_05", &eenumrechits_05, "eenumrechits_05/I");
  dataTree_->Branch("eenumrechits_0530", &eenumrechits_0530, "eenumrechits_0530/I");

  dataTree_->Branch("numtp_eb",            &numtp_eb,            "numtp_eb/I");
  dataTree_->Branch("numtp_ee",            &numtp_ee,            "numtp_ee/I");
  dataTree_->Branch("numtp_geq1_eb",       &numtp_geq1_eb,       "numtp_geq1_eb/I");
  dataTree_->Branch("numtp_geq1_ee",       &numtp_geq1_ee,       "numtp_geq1_ee/I");



  dataTree_->Branch("numtp_samp2_eb",            &numtp_samp2_eb,            "numtp_samp2_eb/I");
  dataTree_->Branch("numtp_samp2_ee",            &numtp_samp2_ee,            "numtp_samp2_ee/I");



  dataTree_->Branch("tpsumet_eb",          &tpsumet_eb,          "tpsumet_eb/F");
  dataTree_->Branch("tpsumet_ee",          &tpsumet_ee,          "tpsumet_ee/F");
 

  dataTree_->Branch("tpsumet_samp2_eb",          &tpsumet_samp2_eb,          "tpsumet_samp2_eb/F");
  dataTree_->Branch("tpsumet_samp2_ee",          &tpsumet_samp2_ee,          "tpsumet_samp2_ee/F");
 

  dataTree_->Branch("tpsumet_geq1_eb",          &tpsumet_geq1_eb,          "tpsumet_geq1_eb/F");
  dataTree_->Branch("tpsumet_geq1_ee",          &tpsumet_geq1_ee,          "tpsumet_geq1_ee/F");
 
  dataTree_->Branch("tpsumet_mod1_eb",          &tpsumet_mod1_eb,          "tpsumet_mod1_eb/F");
  dataTree_->Branch("tpsumet_mod2_eb",          &tpsumet_mod2_eb,          "tpsumet_mod2_eb/F");
  dataTree_->Branch("tpsumet_mod3_eb",          &tpsumet_mod3_eb,          "tpsumet_mod3_eb/F");
  dataTree_->Branch("tpsumet_mod4_eb",          &tpsumet_mod4_eb,          "tpsumet_mod4_eb/F");

  dataTree_->Branch("tpsumet_geq1_mod1_eb",          &tpsumet_geq1_mod1_eb,          "tpsumet_geq1_mod1_eb/F");
  dataTree_->Branch("tpsumet_geq1_mod2_eb",          &tpsumet_geq1_mod2_eb,          "tpsumet_geq1_mod2_eb/F");
  dataTree_->Branch("tpsumet_geq1_mod3_eb",          &tpsumet_geq1_mod3_eb,          "tpsumet_geq1_mod3_eb/F");
  dataTree_->Branch("tpsumet_geq1_mod4_eb",          &tpsumet_geq1_mod4_eb,          "tpsumet_geq1_mod4_eb/F");


  dataTree_->Branch("tpsumet_ring1_ee",          &tpsumet_ring1_ee,          "tpsumet_ring1_ee/F");
  dataTree_->Branch("tpsumet_ring2_ee",          &tpsumet_ring2_ee,          "tpsumet_ring2_ee/F");
  dataTree_->Branch("tpsumet_ring3_ee",          &tpsumet_ring3_ee,          "tpsumet_ring3_ee/F");

  dataTree_->Branch("tpsumet_geq1_ring1_ee",          &tpsumet_geq1_ring1_ee,          "tpsumet_geq1_ring1_ee/F");
  dataTree_->Branch("tpsumet_geq1_ring2_ee",          &tpsumet_geq1_ring2_ee,          "tpsumet_geq1_ring2_ee/F");
  dataTree_->Branch("tpsumet_geq1_ring3_ee",          &tpsumet_geq1_ring3_ee,          "tpsumet_geq1_ring3_ee/F");




  dataTree_->Branch("numcalotower",       &numcalotower,      "numcalotower/I");
  dataTree_->Branch("numtrigprim",        &numtrigprim,       "numtrigprim/I");
  dataTree_->Branch("numrechittower",     &numrechittower,    "numrechittower/I");

  dataTree_->Branch("numcalotower_eb",    &numcalotower_eb,   "numcalotower_eb/I");
  dataTree_->Branch("numtrigprim_eb",     &numtrigprim_eb,    "numtrigprim_eb/I");
  dataTree_->Branch("numrechittower_eb",  &numrechittower_eb, "numrechittower_eb/I");

  dataTree_->Branch("numcalotower_ee",    &numcalotower_ee,   "numcalotower_ee/I");
  dataTree_->Branch("numtrigprim_ee",     &numtrigprim_ee,    "numtrigprim_ee/I");
  dataTree_->Branch("numtrigprim_ee24",     &numtrigprim_ee24,    "numtrigprim_ee24/I");
  dataTree_->Branch("numrechittower_ee",  &numrechittower_ee, "numrechittower_ee/I");
  dataTree_->Branch("numrechittower_ee24",  &numrechittower_ee24, "numrechittower_ee24/I");

  
  dataTree_->Branch("numcalotower_notmatched1",    &numcalotower_notmatched1,    "numcalotower_notmatched1/I");
  dataTree_->Branch("numtrigprim_notmatched1",     &numtrigprim_notmatched1,     "numtrigprim_notmatched1/I");
  dataTree_->Branch("numcalotower_notmatched1_eb", &numcalotower_notmatched1_eb, "numcalotower_notmatched1_eb/I");
  dataTree_->Branch("numtrigprim_notmatched1_eb",  &numtrigprim_notmatched1_eb,  "numtrigprim_notmatched1_eb/I");
  dataTree_->Branch("numcalotower_notmatched1_ee", &numcalotower_notmatched1_ee, "numtrigprim_notmatched1_ee/I");
  dataTree_->Branch("numtrigprim_notmatched1_ee",  &numtrigprim_notmatched1_ee,  "numcalotower_notmatched1_ee/I");
 

  dataTree_->Branch("numcalotower_notmatched2",    &numcalotower_notmatched2,    "numcalotower_notmatched2/I");
  dataTree_->Branch("numtrigprim_notmatched2",     &numtrigprim_notmatched2,     "numtrigprim_notmatched2/I");
  dataTree_->Branch("numcalotower_notmatched2_eb", &numcalotower_notmatched2_eb, "numcalotower_notmatched2_eb/I");
  dataTree_->Branch("numtrigprim_notmatched2_eb",  &numtrigprim_notmatched2_eb,  "numtrigprim_notmatched2_eb/I");
  dataTree_->Branch("numcalotower_notmatched2_ee", &numcalotower_notmatched2_ee, "numtrigprim_notmatched2_ee/I");
  dataTree_->Branch("numtrigprim_notmatched2_ee",  &numtrigprim_notmatched2_ee,  "numcalotower_notmatched2_ee/I");
 
  dataTree_->Branch("numcalotower_matched",    &numcalotower_matched,    "numcalotower_matched/I");
  dataTree_->Branch("numtrigprim_matched",     &numtrigprim_matched,     "numtrigprim_matched/I");
  dataTree_->Branch("numcalotower_matched_eb", &numcalotower_matched_eb, "numcalotower_matched_eb/I");
  dataTree_->Branch("numtrigprim_matched_eb",  &numtrigprim_matched_eb,  "numtrigprim_matched_eb/I");
  dataTree_->Branch("numcalotower_matched_ee", &numcalotower_matched_ee, "numcalotower_matched_ee/I");
  dataTree_->Branch("numtrigprim_matched_ee",  &numtrigprim_matched_ee,  "numtrigprim_matched_ee/I");
 



  dataTree_->Branch("numcalotower_matched2",    &numcalotower_matched2,    "numcalotower_matched2/I");
  dataTree_->Branch("numtrigprim_matched2",     &numtrigprim_matched2,     "numtrigprim_matched2/I");
  dataTree_->Branch("numcalotower_matched_eb2", &numcalotower_matched_eb2, "numcalotower_matched_eb2/I");
  dataTree_->Branch("numtrigprim_matched_eb2",  &numtrigprim_matched_eb2,  "numtrigprim_matched_eb2/I");
  dataTree_->Branch("numcalotower_matched_ee2", &numcalotower_matched_ee2, "numcalotower_matched_ee2/I");
  dataTree_->Branch("numtrigprim_matched_ee2",  &numtrigprim_matched_ee2,  "numtrigprim_matched_ee2/I");
 
  dataTree_->Branch("numrechittower_matched_rh",    &numrechittower_matched_rh,    "numrechittower_matched_rh/I");
  dataTree_->Branch("numtrigprim_matched_rh",     &numtrigprim_matched_rh,     "numtrigprim_matched_rh/I");
  dataTree_->Branch("numrechittower_matched_eb_rh", &numrechittower_matched_eb_rh, "numrechittower_matched_eb_rh/I");
  dataTree_->Branch("numtrigprim_matched_eb_rh",  &numtrigprim_matched_eb_rh,  "numtrigprim_matched_eb_rh/I");
  dataTree_->Branch("numrechittower_matched_ee_rh", &numrechittower_matched_ee_rh, "numrechit_matched_ee_rh/I");
  dataTree_->Branch("numtrigprim_matched_ee_rh",  &numtrigprim_matched_ee_rh,  "numctrigprim_matched_ee_rh/I");

  dataTree_->Branch("numrechittower_matched_ee24_rh", &numrechittower_matched_ee24_rh, "numrechit_matched_ee24_rh/I");
  dataTree_->Branch("numtrigprim_matched_ee24_rh",  &numtrigprim_matched_ee24_rh,  "numctrigprim_matched_ee24_rh/I");
 


  dataTree_->Branch("numrechittower_matched_rh2",    &numrechittower_matched_rh,    "numrechittower_matched_rh2/I");
  dataTree_->Branch("numtrigprim_matched_rh2",     &numtrigprim_matched_rh,     "numtrigprim_matched_rh2/I");
  dataTree_->Branch("numrechittower_matched_eb_rh2", &numrechittower_matched_eb_rh, "numrechittower_matched_eb_rh2/I");
  dataTree_->Branch("numtrigprim_matched_eb_rh2",  &numtrigprim_matched_eb_rh,  "numtrigprim_matched_eb_rh2/I");
  dataTree_->Branch("numrechittower_matched_ee_rh2", &numrechittower_matched_ee_rh, "numrechit_matched_ee_rh2/I");
  dataTree_->Branch("numtrigprim_matched_ee_rh2",  &numtrigprim_matched_ee_rh,  "numctrigprim_matched_ee_rh2/I");
 



 

  dataTree_->Branch("sumcalotower",       &sumcalotower,       "sumcalotower/F");
  dataTree_->Branch("sumtrigprim",        &sumtrigprim,        "sumtrigprim/F");
  dataTree_->Branch("sumrechittower",     &sumrechittower,     "sumrechittower/F");
  dataTree_->Branch("sumcalotower_eb",    &sumcalotower_eb,    "sumcalotower_eb/F");
  dataTree_->Branch("sumtrigprim_eb",     &sumtrigprim_eb,     "sumtrigprim_eb/F");
  dataTree_->Branch("sumrechittower_eb",  &sumrechittower_eb,  "sumrechittower_eb/F");
  dataTree_->Branch("sumcalotower_ee",    &sumcalotower_ee,    "sumcalotower_ee/F");
  dataTree_->Branch("sumtrigprim_ee",     &sumtrigprim_ee,     "sumtrigprim_ee/F");
  dataTree_->Branch("sumrechittower_ee",  &sumrechittower_ee,  "sumrechittower_ee/F");
  dataTree_->Branch("sumtrigprim_ee24",     &sumtrigprim_ee24,     "sumtrigprim_ee24/F");
  dataTree_->Branch("sumrechittower_ee24",  &sumrechittower_ee24,  "sumrechittower_ee24/F");

  dataTree_->Branch("sumtrigprim_samp2_eb",     &sumtrigprim_samp2_eb,     "sumtrigprim_samp2_eb/F");
  dataTree_->Branch("sumtrigprim_samp2_ee",     &sumtrigprim_samp2_ee,     "sumtrigprim_samp2_ee/F");

  
  dataTree_->Branch("sumcalotower_notmatched1",    &sumcalotower_notmatched1,    "sumcalotower_notmatched1/F");
  dataTree_->Branch("sumtrigprim_notmatched1",     &sumtrigprim_notmatched1,     "sumtrigprim_notmatched1/F");
  dataTree_->Branch("sumcalotower_notmatched1_eb", &sumcalotower_notmatched1_eb, "sumcalotower_notmatched1_eb/F");
  dataTree_->Branch("sumtrigprim_notmatched1_eb",  &sumtrigprim_notmatched1_eb,  "sumtrigprim_notmatched1_eb/F");
  dataTree_->Branch("sumcalotower_notmatched1_ee", &sumcalotower_notmatched1_ee, "sumtrigprim_notmatched1_ee/F");
  dataTree_->Branch("sumtrigprim_notmatched1_ee",  &sumtrigprim_notmatched1_ee,  "sumcalotower_notmatched1_ee/F");
 

  dataTree_->Branch("sumcalotower_notmatched2",    &sumcalotower_notmatched2,    "sumcalotower_notmatched2/F");
  dataTree_->Branch("sumtrigprim_notmatched2",     &sumtrigprim_notmatched2,     "sumtrigprim_notmatched2/F");
  dataTree_->Branch("sumcalotower_notmatched2_eb", &sumcalotower_notmatched2_eb, "sumcalotower_notmatched2_eb/F");
  dataTree_->Branch("sumtrigprim_notmatched2_eb",  &sumtrigprim_notmatched2_eb,  "sumtrigprim_notmatched2_eb/F");
  dataTree_->Branch("sumcalotower_notmatched2_ee", &sumcalotower_notmatched2_ee, "sumtrigprim_notmatched2_ee/F");
  dataTree_->Branch("sumtrigprim_notmatched2_ee",  &sumtrigprim_notmatched2_ee,  "sumcalotower_notmatched2_ee/F");
 
  dataTree_->Branch("sumcalotower_matched",    &sumcalotower_matched,    "sumcalotower_matched/F");
  dataTree_->Branch("sumtrigprim_matched",     &sumtrigprim_matched,     "sumtrigprim_matched/F");
  dataTree_->Branch("sumcalotower_matched_eb", &sumcalotower_matched_eb, "sumcalotower_matched_eb/F");
  dataTree_->Branch("sumtrigprim_matched_eb",  &sumtrigprim_matched_eb,  "sumtrigprim_matched_eb/F");
  dataTree_->Branch("sumcalotower_matched_ee", &sumcalotower_matched_ee, "sumtrigprim_matched_ee/F");
  dataTree_->Branch("sumtrigprim_matched_ee",  &sumtrigprim_matched_ee,  "sumcalotower_matched_ee/F");
 

 


   dataTree_->Branch("sumcalotower_matched2", &sumcalotower_matched2, "sumcalotower_matched2/F");
   dataTree_->Branch("sumtrigprim_matched2", &sumtrigprim_matched2, "sumtrigprim_matched2/F");
   dataTree_->Branch("sumcalotower_matched_eb2", &sumcalotower_matched_eb2, "sumcalotower_matched_eb/F");
   dataTree_->Branch("sumtrigprim_matched_eb2", &sumtrigprim_matched_eb2, "sumtrigprim_matched_eb2/F");
   dataTree_->Branch("sumcalotower_matched_ee2", &sumcalotower_matched_ee2, "sumcalotower_matched_ee2/F");
   dataTree_->Branch("sumtrigprim_matched_ee2", &sumtrigprim_matched_ee2, "sumtrigprim_matched_ee2/F");
 

   dataTree_->Branch("sumrechittower_matched_rh", &sumrechittower_matched_rh, "sumrechittower_matched_rh/F");
   dataTree_->Branch("sumtrigprim_matched_rh", &sumtrigprim_matched_rh, "sumtrigprim_matched_rh/F");
   dataTree_->Branch("sumrechittower_matched_eb_rh", &sumrechittower_matched_eb_rh, "sumrechittower_matched_eb_rh/F");
   dataTree_->Branch("sumtrigprim_matched_eb_rh", &sumtrigprim_matched_eb_rh, "sumtrigprim_matched_eb_rh/F");
   dataTree_->Branch("sumrechittower_matched_ee_rh", &sumrechittower_matched_ee_rh, "sumrechittower_matched_ee_rh/F");
   dataTree_->Branch("sumtrigprim_matched_ee_rh", &sumtrigprim_matched_ee_rh, "sumtrigprim_matched_ee_rh/F");
 
   dataTree_->Branch("sumrechittower_matched_ee24_rh", &sumrechittower_matched_ee24_rh, "sumrechittower_matched_ee24_rh/F");
   dataTree_->Branch("sumtrigprim_matched_ee24_rh", &sumtrigprim_matched_ee24_rh, "sumtrigprim_matched_ee24_rh/F");
 
   dataTree_->Branch("sumrechittower_matched_rh2", &sumrechittower_matched_rh2, "sumrechittower_matched_rh2/F");
   dataTree_->Branch("sumtrigprim_matched_rh2", &sumtrigprim_matched_rh2, "sumtrigprim_matched_rh2/F");
   dataTree_->Branch("sumrechittower_matched_eb_rh2", &sumrechittower_matched_eb_rh2, "sumrechittower_matched_eb_rh2/F");
   dataTree_->Branch("sumtrigprim_matched_eb_rh2", &sumtrigprim_matched_eb_rh2, "sumtrigprim_matched_eb_rh2/F");
   dataTree_->Branch("sumrechittower_matched_ee_rh2", &sumrechittower_matched_ee_rh2, "sumrechittower_matched_ee_rh2/F");
   dataTree_->Branch("sumtrigprim_matched_ee_rh2", &sumtrigprim_matched_ee_rh2, "sumtrigprim_matched_ee_rh2/F");



   dataTree_->Branch("sumtrigprim_ee_eta1820", &sumtrigprim_ee_eta1820, "sumtrigprim_ee_eta1820/F");
   dataTree_->Branch("sumtrigprim_ee_eta2122", &sumtrigprim_ee_eta2122, "sumtrigprim_ee_eta2122/F");
   dataTree_->Branch("sumtrigprim_ee_eta2324", &sumtrigprim_ee_eta2324, "sumtrigprim_ee_eta2324/F");
   dataTree_->Branch("sumtrigprim_ee_eta2526", &sumtrigprim_ee_eta2526, "sumtrigprim_ee_eta2526/F");
   dataTree_->Branch("sumtrigprim_ee_eta2728", &sumtrigprim_ee_eta2728, "sumtrigprim_ee_eta2728/F");

   dataTree_->Branch("sumtrigprim_ee_matched_eta1820", &sumtrigprim_ee_matched_eta1820, "sumtrigprim_ee_matched_eta1820/F");
   dataTree_->Branch("sumtrigprim_ee_matched_eta2122", &sumtrigprim_ee_matched_eta2122, "sumtrigprim_ee_matched_eta2122/F");
   dataTree_->Branch("sumtrigprim_ee_matched_eta2324", &sumtrigprim_ee_matched_eta2324, "sumtrigprim_ee_matched_eta2324/F");
   dataTree_->Branch("sumtrigprim_ee_matched_eta2526", &sumtrigprim_ee_matched_eta2526, "sumtrigprim_ee_matched_eta2526/F");
   dataTree_->Branch("sumtrigprim_ee_matched_eta2728", &sumtrigprim_ee_matched_eta2728, "sumtrigprim_ee_matched_eta2728/F");

   dataTree_->Branch("sumrechittower_ee_matched_eta1820", &sumrechittower_ee_matched_eta1820, "sumrechittower_ee_matched_eta1820/F");
   dataTree_->Branch("sumrechittower_ee_matched_eta2122", &sumrechittower_ee_matched_eta2122, "sumrechittower_ee_matched_eta2122/F");
   dataTree_->Branch("sumrechittower_ee_matched_eta2324", &sumrechittower_ee_matched_eta2324, "sumrechittower_ee_matched_eta2324/F");
   dataTree_->Branch("sumrechittower_ee_matched_eta2526", &sumrechittower_ee_matched_eta2526, "sumrechittower_ee_matched_eta2526/F");
   dataTree_->Branch("sumrechittower_ee_matched_eta2728", &sumrechittower_ee_matched_eta2728, "sumrechittower_ee_matched_eta2728/F");


   dataTree_->Branch("sumtrigprim_ee_bxm2", &sumtrigprim_ee_bxm2, "sumtrigprim_ee_bxm2/F");
   dataTree_->Branch("sumtrigprim_ee_bxm1", &sumtrigprim_ee_bxm1, "sumtrigprim_ee_bxm1/F");
   dataTree_->Branch("sumtrigprim_ee_bx0",  &sumtrigprim_ee_bx0,  "sumtrigprim_ee_bx0/F");
   dataTree_->Branch("sumtrigprim_ee_bxp1", &sumtrigprim_ee_bxp1, "sumtrigprim_ee_bxp1/F");
   dataTree_->Branch("sumtrigprim_ee_bxp2", &sumtrigprim_ee_bxp2, "sumtrigprim_ee_bxp2/F");



















































    
 





}
//////////////////////////////////////////////////////////////////////////////////////////
void PFDataTreeProducer::endJob() 
{
  if (file_ !=0) 
    {
      file_->cd();

      h_ttonline_et->Write();
      h_ttoffline_et->Write();

      histo_event_->Write();
      eb_rechitenergy_->Write();
      ee_rechitenergy_->Write();

      ee_rechitenergy_notypeb_->Write();

      calotower_map_cumul->Write();
      tpmap_cumul->Write();
      tpmap_cumul2->Write();
      rechitmap_cumul->Write();

      tower_map_cumul_matched->Write();
      tower_map_cumul_matched2->Write();
      tower_map_cumul_matched_rh->Write();
      tower_map_cumul_matched_rh2->Write();

      tower_map_cumul_unmatched1->Write();
      tower_map_cumul_unmatched2->Write();


      ct_et_all->Write();
      ct_et_matched->Write();
      ct_et_matched2->Write();
      ct_et_unmatched->Write();

    rh_et_matched_rh->Write();
    rh_et_matched_rh2->Write();

    rh_et_all_ee->Write();
    rh_et_all_ee24->Write();

    tp_et_all_ee->Write();
    tp_et_all_ee24->Write();


    tp_et_all->Write();
    tp_et_matched->Write();
    tp_et_matched2->Write();
 
    tp_et_matched_rh->Write();
    tp_et_matched_rh2->Write();
    tp_et_unmatched->Write();


    ctminustp_et_matched->Write();
    ctminustp_et_matched2->Write();

    rhminustp_et_matched_rh->Write();
    rhminustp_et_matched_rh2->Write();



    ct_et_all_eb->Write();
    ct_et_matched_eb->Write();
    ct_et_matched_eb2->Write();


    rh_et_matched_eb_rh->Write();
    rh_et_matched_eb_rh2->Write();


    rh_et_matched_ee24_rh->Write();


    ct_et_unmatched_eb->Write();

    tp_et_all_eb->Write();
    tp_et_matched_eb->Write();
    tp_et_matched_eb2->Write();
    tp_et_matched_eb_rh->Write();
    tp_et_matched_eb_rh2->Write();
    tp_et_unmatched_eb->Write();


    ctminustp_et_matched_eb->Write();
    ctminustp_et_matched_eb2->Write();

    rhminustp_et_matched_eb_rh->Write();
    rhminustp_et_matched_eb_rh2->Write();


    ct_et_all_ee->Write();
    ct_et_matched_ee->Write();
    ct_et_matched_ee2->Write();
    ct_et_unmatched_ee->Write();

    rh_et_matched_ee_rh->Write();
    rh_et_matched_ee_rh2->Write();



    tp_et_all_ee->Write();
    tp_et_matched_ee->Write();
    tp_et_matched_ee2->Write();
    tp_et_matched_ee_rh->Write();
    tp_et_matched_ee24_rh->Write();
    tp_et_matched_ee_rh2->Write();
    tp_et_unmatched_ee->Write();


    ctminustp_et_matched_ee->Write();
    ctminustp_et_matched_ee2->Write();

    rhminustp_et_matched_ee_rh->Write();
    rhminustp_et_matched_ee_rh2->Write();



    ct_et_vs_ieta_all->Write();
    ct_et_vs_ieta_matched->Write();
    ct_et_vs_ieta_matched2->Write();
    ct_et_vs_ieta_unmatched->Write();

    rh_et_vs_ieta_matched_rh->Write();
    rh_et_vs_ieta_matched_rh2->Write();



    matched_timevset->Write();
    unmatched_timevset->Write();

    matched_timevsieta->Write();
    unmatched_timevsieta->Write();

    matched_timevsieta_etweighted->Write();
    unmatched_timevsieta_etweighted->Write();
 
    matched_timevsetvsieta->Write();
    unmatched_timevsetvsieta->Write();

    ietamatched->Write();
    ietatotal->Write();


    tp_et_vs_ieta_all->Write();
    tp_et_vs_ieta_matched->Write();
    tp_et_vs_ieta_matched2->Write();
    tp_et_vs_ieta_matched_rh->Write();
    tp_et_vs_ieta_matched_rh2->Write();

    tp_et_vs_ieta_unmatched->Write();


    ctminustp_et_vs_ieta_matched->Write();
    ctminustp_et_vs_ieta_matched2->Write();


    rhminustp_et_vs_ieta_matched_rh->Write();
    rhminustp_et_vs_ieta_matched_rh2->Write();

    tp_et_vs_ct_et_matched->Write();
    tp_et_vs_ct_et_matched2->Write();

    tp_et_vs_ct_et_matched_eb->Write();
    tp_et_vs_ct_et_matched_ee->Write();
 
    tp_et_vs_ct_et_matched_eb2->Write();
    tp_et_vs_ct_et_matched_ee2->Write();




    tp_et_vs_rh_et_matched_rh->Write();
    tp_et_vs_rh_et_matched_rh2->Write();

    tp_et_vs_rh_et_matched_eb_rh->Write();
    tp_et_vs_rh_et_matched_ee_rh->Write();
 
    tp_et_vs_rh_et_matched_eb_rh2->Write();
    tp_et_vs_rh_et_matched_ee_rh2->Write();

    tp_et_vs_rh_et_matched_ee_ieta1820_rh->Write();
    tp_et_vs_rh_et_matched_ee_ieta2122_rh->Write();
    tp_et_vs_rh_et_matched_ee_ieta2324_rh->Write();
    tp_et_vs_rh_et_matched_ee_ieta2526_rh->Write();
    tp_et_vs_rh_et_matched_ee_ieta2728_rh->Write();





      eb_rechitenergy_02->Write();
      eb_rechitenergy_04->Write();
      eb_rechitenergy_06->Write();
      eb_rechitenergy_08->Write();
      eb_rechitenergy_10->Write();
      eb_rechitenergy_12->Write();
      eb_rechitenergy_14->Write();
      eb_rechitenergy_148->Write();
      
      
      
      ee_rechitet_16->Write();
      ee_rechitet_18->Write();
      ee_rechitet_20->Write();
      ee_rechitet_22->Write();
      ee_rechitet_24->Write();
      ee_rechitet_26->Write();
      ee_rechitet_28->Write();
      ee_rechitet_30->Write();
      



      eb_rechitet_02->Write();
      eb_rechitet_04->Write();
      eb_rechitet_06->Write();
      eb_rechitet_08->Write();
      eb_rechitet_10->Write();
      eb_rechitet_12->Write();
      eb_rechitet_14->Write();
      eb_rechitet_148->Write();
      
      
      
      ee_rechitenergy_16->Write();
      ee_rechitenergy_18->Write();
      ee_rechitenergy_20->Write();
      ee_rechitenergy_22->Write();
      ee_rechitenergy_24->Write();
      ee_rechitenergy_26->Write();
      ee_rechitenergy_28->Write();
      ee_rechitenergy_30->Write();
      



      eb_rechitetvspu_05->Write();
      eb_rechitetvspu_10->Write();
      eb_rechitetvspu_15->Write();

      ee_rechitetvspu_20->Write();
      ee_rechitetvspu_25->Write();
      ee_rechitetvspu_30->Write();

      eb_rechitet_->Write();
      ee_rechitet_->Write();

      eb_rechiten_vs_eta->Write();
  
      eb_rechitet_vs_eta->Write();
  



      eep_rechiten_vs_eta->Write();
      eep_rechiten_vs_phi->Write();

      eem_rechiten_vs_eta->Write();
      eem_rechiten_vs_phi->Write();

      eep_rechitet_vs_eta->Write();
      eep_rechitet_vs_phi->Write();

      eem_rechitet_vs_eta->Write();
      eem_rechitet_vs_phi->Write();


      ebocc->Write();
      eboccgt1->Write();
      eboccgt1et->Write();
      eboccet->Write();
      eboccetgt1et->Write();
      eboccen->Write();
      eboccengt1->Write();
 
      eeocc->Write();
      eeoccgt1->Write();
      eeoccgt1et->Write();
      eeoccet->Write();
      eeoccetgt1et->Write();
      eeoccen->Write();
      eeoccengt1->Write();


      eb_timing_0->Write();
      eb_timing_200->Write();
      eb_timing_400->Write();
      eb_timing_600->Write();
      eb_timing_800->Write();
      eb_timing_1000->Write();
      eb_timing_2000->Write();
      eb_timing_3000->Write();
      eb_timing_5000->Write();

      eb_r4_0->Write();
      eb_r4_200->Write();
      eb_r4_400->Write();
      eb_r4_600->Write();
      eb_r4_800->Write();
      eb_r4_1000->Write();
      eb_r4_2000->Write();
      eb_r4_3000->Write();
      eb_r4_5000->Write();


      eb_timing_r4_0->Write();
      eb_timing_r4_200->Write();
      eb_timing_r4_400->Write();
      eb_timing_r4_600->Write();
      eb_timing_r4_800->Write();
      eb_timing_r4_1000->Write();
      eb_timing_r4_2000->Write();
      eb_timing_r4_3000->Write();
      eb_timing_r4_5000->Write();

      eb_timing_vs_r4_0->Write();
      eb_timing_vs_r4_200->Write();
      eb_timing_vs_r4_400->Write();
      eb_timing_vs_r4_600->Write();
      eb_timing_vs_r4_800->Write();
      eb_timing_vs_r4_1000->Write();
      eb_timing_vs_r4_2000->Write();
      eb_timing_vs_r4_3000->Write();
      eb_timing_vs_r4_5000->Write();



      rechiteta_all->Write();
      rechiteta_gt1et->Write();
      rechiteta_etweight->Write();
      rechiteta_gt1et_pu10->Write();
      rechiteta_gt1et_pu20->Write();
      rechiteta_gt1et_pu30->Write();


      calotowereta_all->Write();
      calotowereta_gt1et->Write();
      calotowereta_etweight->Write();
      calotowereta_gt1et_pu10->Write();
      calotowereta_gt1et_pu20->Write();
      calotowereta_gt1et_pu30->Write();

      sceta_all->Write();
      sceta_severity0->Write();
      sceta_koot0->Write();



      sceta_all_pueq01->Write();
      sceta_severity0_pueq01->Write();

      sceta_all_pueq02->Write();
      sceta_severity0_pueq02->Write();

      sceta_all_pueq03->Write();
      sceta_severity0_pueq03->Write();

      sceta_all_pueq04->Write();
      sceta_severity0_pueq04->Write();

      sceta_all_pueq05->Write();
      sceta_severity0_pueq05->Write();

      sceta_all_pueq06->Write();
      sceta_severity0_pueq06->Write();

      sceta_all_pueq07->Write();
      sceta_severity0_pueq07->Write();
      
      sceta_all_pueq08->Write();
      sceta_severity0_pueq08->Write();

      sceta_all_pueq09->Write();
      sceta_severity0_pueq09->Write();


      sceta_all_gt2->Write();
      sceta_severity0_gt2->Write();
      sceta_koot0_gt2->Write();


      sceta_all_gt5->Write();
      sceta_severity0_gt5->Write();
      sceta_koot0_gt5->Write();


      sceta_all_gt10->Write();
      sceta_severity0_gt10->Write();
      sceta_koot0_gt10->Write();


      scet_eb_all->Write();
      scet_eb_severity0->Write();
      scet_eb_koot0->Write();

      scet_ee_all->Write();
      scet_ee_severity0->Write();
      scet_ee_koot0->Write();

      
      scsumet_eb_all->Write();
      scsumet_eb_severity0->Write();
      scsumet_eb_koot0->Write();
   
      scsumet_ee_all->Write();
      scsumet_ee_severity0->Write();
      scsumet_ee_koot0->Write();


      scsumet_eb_all_gt2->Write();
      scsumet_eb_severity0_gt2->Write();
      scsumet_eb_koot0_gt2->Write();
   
      scsumet_ee_all_gt2->Write();
      scsumet_ee_severity0_gt2->Write();
      scsumet_ee_koot0_gt2->Write();

      scsumet_eb_all_gt5->Write();
      scsumet_eb_severity0_gt5->Write();
      scsumet_eb_koot0_gt5->Write();
   
      scsumet_ee_all_gt5->Write();
      scsumet_ee_severity0_gt5->Write();
      scsumet_ee_koot0_gt5->Write();

      scsumet_eb_all_gt10->Write();
      scsumet_eb_severity0_gt10->Write();
      scsumet_eb_koot0_gt10->Write();
   
      scsumet_ee_all_gt10->Write();
      scsumet_ee_severity0_gt10->Write();
      scsumet_ee_koot0_gt10->Write();



      scocc_eb_gt50->Write();
      scocc_ee_gt50->Write();




      scet_eb_all_eta15->Write();
    scet_eb_all_eta20->Write();
    scet_eb_all_eta25->Write();
 
    scet_ee_all_eta15->Write();
    scet_ee_all_eta20->Write();
    scet_ee_all_eta25->Write();



    scet_eb_all_eta15_pu10->Write();
    scet_eb_all_eta20_pu10->Write();
    scet_eb_all_eta25_pu10->Write();
 
    scet_ee_all_eta15_pu10->Write();
    scet_ee_all_eta20_pu10->Write();
    scet_ee_all_eta25_pu10->Write();

    scet_eb_all_eta15_pu20->Write();
    scet_eb_all_eta20_pu20->Write();
    scet_eb_all_eta25_pu20->Write();
 
    scet_ee_all_eta15_pu20->Write();
    scet_ee_all_eta20_pu20->Write();
    scet_ee_all_eta25_pu20->Write();

    scet_eb_all_eta15_pu30->Write();
    scet_eb_all_eta20_pu30->Write();
    scet_eb_all_eta25_pu30->Write();
 
    scet_ee_all_eta15_pu30->Write();
    scet_ee_all_eta20_pu30->Write();
    scet_ee_all_eta25_pu30->Write();


    scet_eb_all_eta15_pueq10->Write();
    scet_eb_all_eta20_pueq10->Write();
    scet_eb_all_eta25_pueq10->Write();
 
    scet_ee_all_eta15_pueq10->Write();
    scet_ee_all_eta20_pueq10->Write();
    scet_ee_all_eta25_pueq10->Write();

    scet_eb_all_eta15_pueq20->Write();
    scet_eb_all_eta20_pueq20->Write();
    scet_eb_all_eta25_pueq20->Write();
 
    scet_ee_all_eta15_pueq20->Write();
    scet_ee_all_eta20_pueq20->Write();
    scet_ee_all_eta25_pueq20->Write();
 



    scsumet_eb_all_eta15->Write();
    scsumet_eb_all_eta20->Write();
    scsumet_eb_all_eta25->Write();
 
    scsumet_ee_all_eta15->Write();
    scsumet_ee_all_eta20->Write();
    scsumet_ee_all_eta25->Write();



    scsumet_eb_all_eta15_pu10->Write();
    scsumet_eb_all_eta20_pu10->Write();
    scsumet_eb_all_eta25_pu10->Write();
 
    scsumet_ee_all_eta15_pu10->Write();
    scsumet_ee_all_eta20_pu10->Write();
    scsumet_ee_all_eta25_pu10->Write();

    scsumet_eb_all_eta15_pu20->Write();
    scsumet_eb_all_eta20_pu20->Write();
    scsumet_eb_all_eta25_pu20->Write();
 
    scsumet_ee_all_eta15_pu20->Write();
    scsumet_ee_all_eta20_pu20->Write();
    scsumet_ee_all_eta25_pu20->Write();

    scsumet_eb_all_eta15_pu30->Write();
    scsumet_eb_all_eta20_pu30->Write();
    scsumet_eb_all_eta25_pu30->Write();
 
    scsumet_ee_all_eta15_pu30->Write();
    scsumet_ee_all_eta20_pu30->Write();
    scsumet_ee_all_eta25_pu30->Write();


    scsumet_eb_all_eta15_pueq10->Write();
    scsumet_eb_all_eta20_pueq10->Write();
    scsumet_eb_all_eta25_pueq10->Write();
 
    scsumet_ee_all_eta15_pueq10->Write();
    scsumet_ee_all_eta20_pueq10->Write();
    scsumet_ee_all_eta25_pueq10->Write();

    scsumet_eb_all_eta15_pueq20->Write();
    scsumet_eb_all_eta20_pueq20->Write();
    scsumet_eb_all_eta25_pueq20->Write();
 
    scsumet_ee_all_eta15_pueq20->Write();
    scsumet_ee_all_eta20_pueq20->Write();
    scsumet_ee_all_eta25_pueq20->Write();
 


    scet_eb_all_eta15_pueq01->Write();
    scet_eb_all_eta20_pueq01->Write();
    scet_eb_all_eta25_pueq01->Write();
 
    scet_ee_all_eta15_pueq01->Write();
    scet_ee_all_eta20_pueq01->Write();
    scet_ee_all_eta25_pueq01->Write();

    scsumet_eb_all_eta15_pueq01->Write();
    scsumet_eb_all_eta20_pueq01->Write();
    scsumet_eb_all_eta25_pueq01->Write();
 
    scsumet_ee_all_eta15_pueq01->Write();
    scsumet_ee_all_eta20_pueq01->Write();
    scsumet_ee_all_eta25_pueq01->Write();



    scet_eb_all_eta15_pueq02->Write();
    scet_eb_all_eta20_pueq02->Write();
    scet_eb_all_eta25_pueq02->Write();
 
    scet_ee_all_eta15_pueq02->Write();
    scet_ee_all_eta20_pueq02->Write();
    scet_ee_all_eta25_pueq02->Write();

    scsumet_eb_all_eta15_pueq02->Write();
    scsumet_eb_all_eta20_pueq02->Write();
    scsumet_eb_all_eta25_pueq02->Write();
 
    scsumet_ee_all_eta15_pueq02->Write();
    scsumet_ee_all_eta20_pueq02->Write();
    scsumet_ee_all_eta25_pueq02->Write();



    scet_eb_all_eta15_pueq03->Write();
    scet_eb_all_eta20_pueq03->Write();
    scet_eb_all_eta25_pueq03->Write();
 
    scet_ee_all_eta15_pueq03->Write();
    scet_ee_all_eta20_pueq03->Write();
    scet_ee_all_eta25_pueq03->Write();

    scsumet_eb_all_eta15_pueq03->Write();
    scsumet_eb_all_eta20_pueq03->Write();
    scsumet_eb_all_eta25_pueq03->Write();
 
    scsumet_ee_all_eta15_pueq03->Write();
    scsumet_ee_all_eta20_pueq03->Write();
    scsumet_ee_all_eta25_pueq03->Write();



    scet_eb_all_eta15_pueq04->Write();
    scet_eb_all_eta20_pueq04->Write();
    scet_eb_all_eta25_pueq04->Write();
 
    scet_ee_all_eta15_pueq04->Write();
    scet_ee_all_eta20_pueq04->Write();
    scet_ee_all_eta25_pueq04->Write();

    scsumet_eb_all_eta15_pueq04->Write();
    scsumet_eb_all_eta20_pueq04->Write();
    scsumet_eb_all_eta25_pueq04->Write();
 
    scsumet_ee_all_eta15_pueq04->Write();
    scsumet_ee_all_eta20_pueq04->Write();
    scsumet_ee_all_eta25_pueq04->Write();



    scet_eb_all_eta15_pueq05->Write();
    scet_eb_all_eta20_pueq05->Write();
    scet_eb_all_eta25_pueq05->Write();
 
    scet_ee_all_eta15_pueq05->Write();
    scet_ee_all_eta20_pueq05->Write();
    scet_ee_all_eta25_pueq05->Write();

    scsumet_eb_all_eta15_pueq05->Write();
    scsumet_eb_all_eta20_pueq05->Write();
    scsumet_eb_all_eta25_pueq05->Write();
 
    scsumet_ee_all_eta15_pueq05->Write();
    scsumet_ee_all_eta20_pueq05->Write();
    scsumet_ee_all_eta25_pueq05->Write();



    scet_eb_all_eta15_pueq06->Write();
    scet_eb_all_eta20_pueq06->Write();
    scet_eb_all_eta25_pueq06->Write();
 
    scet_ee_all_eta15_pueq06->Write();
    scet_ee_all_eta20_pueq06->Write();
    scet_ee_all_eta25_pueq06->Write();

    scsumet_eb_all_eta15_pueq06->Write();
    scsumet_eb_all_eta20_pueq06->Write();
    scsumet_eb_all_eta25_pueq06->Write();
 
    scsumet_ee_all_eta15_pueq06->Write();
    scsumet_ee_all_eta20_pueq06->Write();
    scsumet_ee_all_eta25_pueq06->Write();



    scet_eb_all_eta15_pueq07->Write();
    scet_eb_all_eta20_pueq07->Write();
    scet_eb_all_eta25_pueq07->Write();
 
    scet_ee_all_eta15_pueq07->Write();
    scet_ee_all_eta20_pueq07->Write();
    scet_ee_all_eta25_pueq07->Write();

    scsumet_eb_all_eta15_pueq07->Write();
    scsumet_eb_all_eta20_pueq07->Write();
    scsumet_eb_all_eta25_pueq07->Write();
 
    scsumet_ee_all_eta15_pueq07->Write();
    scsumet_ee_all_eta20_pueq07->Write();
    scsumet_ee_all_eta25_pueq07->Write();



    scet_eb_all_eta15_pueq08->Write();
    scet_eb_all_eta20_pueq08->Write();
    scet_eb_all_eta25_pueq08->Write();
 
    scet_ee_all_eta15_pueq08->Write();
    scet_ee_all_eta20_pueq08->Write();
    scet_ee_all_eta25_pueq08->Write();

    scsumet_eb_all_eta15_pueq08->Write();
    scsumet_eb_all_eta20_pueq08->Write();
    scsumet_eb_all_eta25_pueq08->Write();
 
    scsumet_ee_all_eta15_pueq08->Write();
    scsumet_ee_all_eta20_pueq08->Write();
    scsumet_ee_all_eta25_pueq08->Write();



    scet_eb_all_eta15_pueq09->Write();
    scet_eb_all_eta20_pueq09->Write();
    scet_eb_all_eta25_pueq09->Write();
 
    scet_ee_all_eta15_pueq09->Write();
    scet_ee_all_eta20_pueq09->Write();
    scet_ee_all_eta25_pueq09->Write();

    scsumet_eb_all_eta15_pueq09->Write();
    scsumet_eb_all_eta20_pueq09->Write();
    scsumet_eb_all_eta25_pueq09->Write();
 
    scsumet_ee_all_eta15_pueq09->Write();
    scsumet_ee_all_eta20_pueq09->Write();
    scsumet_ee_all_eta25_pueq09->Write();


    ebtime_vs_bxtrain_01->Write();
    ebtime_vs_bxtrain_05->Write();
    eetime_vs_bxtrain_01->Write();
    eetime_vs_bxtrain_05->Write();



    rechiteta_vs_bxtrain_01->Write();
    rechiteta_vs_bxtrain_05->Write();
    sceta_vs_bxtrain->Write();



    eb_digi_01->Write();
    eb_digi_05->Write();
    eb_digi_0105->Write();
    eb_digi_0530->Write();

    ee_digi_01->Write();
    ee_digi_05->Write();
    ee_digi_0105->Write();
    ee_digi_0530->Write();


    numtp_vs_ieta->Write();
    numtp_vs_ieta_samp2->Write();
    numtp_geq1_vs_ieta->Write();
    numtp_etweighted_vs_ieta->Write();

    
    // normalise 3d digi plots

    for (Int_t i=0;i<120;i++) {
      for (Int_t j=0;j<10;j++) {

	float r=0;

	float a=eb_digi_0105_vs_time->GetBinContent(i+1,j+1);
	float b=eb_digi_0105_vs_time_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_time->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_time->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_time_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_time->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_time->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_time_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_time->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_time->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_time_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_time->SetBinContent(i+1,j+1,r);




	r=0;

	a=eb_digi_0105_vs_time_eta15->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_time_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_time_eta15->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_time_eta15->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_time_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_time_eta15->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_time_eta15->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_time_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_time_eta15->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_time_eta15->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_time_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_time_eta15->SetBinContent(i+1,j+1,r);







	r=0;

	a=eb_digi_0105_vs_time_eta20->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_time_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_time_eta20->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_time_eta20->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_time_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_time_eta20->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_time_eta20->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_time_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_time_eta20->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_time_eta20->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_time_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_time_eta20->SetBinContent(i+1,j+1,r);



	r=0;

	a=eb_digi_0105_vs_time_eta25->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_time_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_time_eta25->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_time_eta25->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_time_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_time_eta25->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_time_eta25->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_time_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_time_eta25->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_time_eta25->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_time_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_time_eta25->SetBinContent(i+1,j+1,r);










      }
    }




    for (Int_t i=0;i<40;i++) {
      for (Int_t j=0;j<10;j++) {

	float r=0;

	float a=eb_digi_0105_vs_bxtrain->GetBinContent(i+1,j+1);
	float b=eb_digi_0105_vs_bxtrain_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_bxtrain->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_bxtrain->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_bxtrain_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_bxtrain->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_bxtrain->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_bxtrain_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_bxtrain->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_bxtrain->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_bxtrain_norm->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_bxtrain->SetBinContent(i+1,j+1,r);





	r=0;

	a=eb_digi_0105_vs_bxtrain_eta15->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_bxtrain_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_bxtrain_eta15->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_bxtrain_eta15->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_bxtrain_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_bxtrain_eta15->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_bxtrain_eta15->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_bxtrain_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_bxtrain_eta15->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_bxtrain_eta15->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_bxtrain_norm_eta15->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_bxtrain_eta15->SetBinContent(i+1,j+1,r);







	r=0;

	a=eb_digi_0105_vs_bxtrain_eta20->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_bxtrain_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_bxtrain_eta20->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_bxtrain_eta20->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_bxtrain_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_bxtrain_eta20->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_bxtrain_eta20->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_bxtrain_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_bxtrain_eta20->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_bxtrain_eta20->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_bxtrain_norm_eta20->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_bxtrain_eta20->SetBinContent(i+1,j+1,r);







	r=0;

	a=eb_digi_0105_vs_bxtrain_eta25->GetBinContent(i+1,j+1);
	b=eb_digi_0105_vs_bxtrain_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0105_vs_bxtrain_eta25->SetBinContent(i+1,j+1,r);


	r=0;

	a=eb_digi_0530_vs_bxtrain_eta25->GetBinContent(i+1,j+1);
	b=eb_digi_0530_vs_bxtrain_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	eb_digi_0530_vs_bxtrain_eta25->SetBinContent(i+1,j+1,r);



	r=0;

	a=ee_digi_0105_vs_bxtrain_eta25->GetBinContent(i+1,j+1);
	b=ee_digi_0105_vs_bxtrain_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0105_vs_bxtrain_eta25->SetBinContent(i+1,j+1,r);


	r=0;

	a=ee_digi_0530_vs_bxtrain_eta25->GetBinContent(i+1,j+1);
	b=ee_digi_0530_vs_bxtrain_norm_eta25->GetBinContent(i+1,j+1);

	if (b>0) r=a/b;

	ee_digi_0530_vs_bxtrain_eta25->SetBinContent(i+1,j+1,r);








      }
    }
    

    eb_digi_0105_vs_time->Write();
    eb_digi_0530_vs_time->Write();


    ee_digi_0105_vs_time->Write();
    ee_digi_0530_vs_time->Write();


    eb_digi_0105_vs_bxtrain->Write();
    eb_digi_0530_vs_bxtrain->Write();


    ee_digi_0105_vs_bxtrain->Write();
    ee_digi_0530_vs_bxtrain->Write();



    eb_digi_0105_vs_time_eta15->Write();
    eb_digi_0530_vs_time_eta15->Write();


    ee_digi_0105_vs_time_eta15->Write();
    ee_digi_0530_vs_time_eta15->Write();


    eb_digi_0105_vs_bxtrain_eta15->Write();
    eb_digi_0530_vs_bxtrain_eta15->Write();


    ee_digi_0105_vs_bxtrain_eta15->Write();
    ee_digi_0530_vs_bxtrain_eta15->Write();






    eb_digi_0105_vs_time_eta20->Write();
    eb_digi_0530_vs_time_eta20->Write();


    ee_digi_0105_vs_time_eta20->Write();
    ee_digi_0530_vs_time_eta20->Write();


    eb_digi_0105_vs_bxtrain_eta20->Write();
    eb_digi_0530_vs_bxtrain_eta20->Write();


    ee_digi_0105_vs_bxtrain_eta20->Write();
    ee_digi_0530_vs_bxtrain_eta20->Write();






    eb_digi_0105_vs_time_eta25->Write();
    eb_digi_0530_vs_time_eta25->Write();


    ee_digi_0105_vs_time_eta25->Write();
    ee_digi_0530_vs_time_eta25->Write();


    eb_digi_0105_vs_bxtrain_eta25->Write();
    eb_digi_0530_vs_bxtrain_eta25->Write();


    ee_digi_0105_vs_bxtrain_eta25->Write();
    ee_digi_0530_vs_bxtrain_eta25->Write();






    dataTree_   ->Write();
    }
  file_ = 0;
}


void PFDataTreeProducer::scan5x5(const DetId & det, const edm::Handle<EcalRecHitCollection> &hits, const edm::ESHandle<CaloTopology>  &caloTopo, const edm::ESHandle<CaloGeometry>  &geometry, int &nHits, float & totEt)
{
  nHits = 0; 
  totEt = 0;

  CaloNavigator<DetId> cursor = CaloNavigator<DetId>(det,caloTopo->getSubdetectorTopology(det));

  for(int j=side_/2; j>=-side_/2; --j)
    {
      for(int i=-side_/2; i<=side_/2; ++i)
	{
	  cursor.home();
	  cursor.offsetBy(i,j);
	  if(hits->find(*cursor)!=hits->end()) // if hit exists in the rechit collection
	    {
	      EcalRecHit tmpHit = *hits->find(*cursor); // get rechit with detID at cursor

	
	      const GlobalPoint p ( geometry->getPosition(*cursor) ) ; // calculate Et of the rechit
	      TVector3 hitPos(p.x(),p.y(),p.z());
	      hitPos *= 1.0/hitPos.Mag();
	      hitPos *= tmpHit.energy();
	      float rechitEt =  hitPos.Pt();
        
	      //--- return values
	      totEt += rechitEt;  
   

	      if(tmpHit.energy()>Emin_ && !tmpHit.checkFlag(EcalRecHit::kGood))nHits++;

	    }
	}
    }


}


//////////////////////////////////////////////////////////////////////////////////////////
void PFDataTreeProducer::analyze(edm::Event const& event, edm::EventSetup const& iSetup) 
{ 

  /*
   edm::Handle<L1GlobalTriggerReadoutRecord> gtRecord;
   event.getByLabel("gtDigis", gtRecord);
  */


  //  edm::ESHandle<EcalLaserDbService> laser;
  // iSetup.get<EcalLaserDbRecord>().get(laser);
  
  calotower_map->Reset();
  tp_map->Reset();
  tp_map_samp2->Reset();
  tp_map2->Reset();
  rechit_map->Reset();
  rechit_map_time->Reset();

  rechit_map_ixiy->Reset();
  rechit_map_time_ixiy->Reset();
   
   bit36=0;
   bit37=0;
   bit38=0;
   bit39=0;
   bit40=0;
   bit41=0;
   bit3 =0;
   bit4 =0;
   bit9 =0;
   bit0 =0;

   eg1=0;
   eg2=0;
   eg5=0;
   algo124=0;

   algo54=0;
   algo55=0;
   algo56=0;
   algo57=0;
   algo58=0;
   algo59=0;
   algo60=0;
   algo61=0;
   algo62=0;

   algo106=0;
   algo107=0;



   rank_=0;
   ntrk=0;
   goodtrk=0;

   numvtx=0;
   vtx_x=0;
   vtx_y=0;
   vtx_z=0;
   vtx_x_err=0;
   vtx_y_err=0;
   vtx_z_err=0;

   vtx_chi2=-1;
   vtx_ndof=0;
   vtx_ntracks=0;
   vtx_isfake=0;

   vtx_good=0;
   numgoodvtx=0;

   scale=0;

   energy_ecal=0;
   energy_hcal=0;


  eemax=0;
  eemaxet=0;
  eeix=0;
  eeiy=0;
  eeix=0;
  eetime=0;
  eeflags=0;
  eehits=0;
  eehits1GeV=0;
  ee_eta=0;
  ee_phi=0;
  ee_r9=0;
  eephits=0;
  eemhits=0;

  ebmax=0;
  eb_ieta=0;
  eb_iphi=0;
  eb_eta=0;
  eb_phi=0;
  ebtime=0;
  ebflags=0;
  ebhits=0;
  ebhits1GeV=0;
  ebmaxet=0;
  eb_r9=0;
  eb_r4=0;
  eb_e9=0;
  eb_e25=0;


  ebchi2=0;
  eb2chi2=0;
  ebchi2oot=0;
  eb2chi2oot=0;


  eemax2=0;
  eemaxet2=0;
  eeix2=0;
  eeiy2=0;
  eeix2=0;
  eetime2=0;
  eeflags2=0;
  eehits1GeV2=0;
  ee_eta2=0;
  ee_phi2=0;
  ee_r92=0;

  ebmax2=0;
  eb_ieta2=0;
  eb_iphi2=0;
  eb_eta2=0;
  eb_phi2=0;
  ebtime2=0;
  ebflags2=0;
  ebhits1GeV2=0;
  ebmaxet2=0;
  eb_r92=0;
  eb_r42=0;




  sumtrigprim_ee_eta1820=0;
  sumtrigprim_ee_eta2122=0;
  sumtrigprim_ee_eta2324=0;
  sumtrigprim_ee_eta2526=0;
  sumtrigprim_ee_eta2728=0;

  sumtrigprim_ee_matched_eta1820=0;
  sumtrigprim_ee_matched_eta2122=0;
  sumtrigprim_ee_matched_eta2324=0;
  sumtrigprim_ee_matched_eta2526=0;
  sumtrigprim_ee_matched_eta2728=0;


  sumrechittower_ee_matched_eta1820=0;
  sumrechittower_ee_matched_eta2122=0;
  sumrechittower_ee_matched_eta2324=0;
  sumrechittower_ee_matched_eta2526=0;
  sumrechittower_ee_matched_eta2728=0;

  sumtrigprim_ee_bxm2=0;
  sumtrigprim_ee_bxm1=0;
  sumtrigprim_ee_bx0=0;
  sumtrigprim_ee_bxp1=0;
  sumtrigprim_ee_bxp2=0;



  ebnumrechits_01=0;
  ebnumrechits_05=0;
  ebnumrechits_0105=0;
  ebnumrechits_0530=0;

  eenumrechits_01=0;
  eenumrechits_05=0;
  eenumrechits_0105=0;
  eenumrechits_0530=0;



  ebflag_kgood=0;
  ebflag_kpoorreco=0;
  ebflag_koutoftime=0;
  ebflag_kfake=0;


  tmean_en=0;
  terr_en=0;

  tmean_sig=0;
  terr_sig=0;


  r4count=0;
  e2e9count_thresh0=0;
  e2e25count_thresh0=0;
  e2e9count_thresh1=0;
  e2e25count_thresh1=0;

  e2e9count_thresh0_nor4=0;
  e2e25count_thresh0_nor4=0;
  e2e9count_thresh1_nor4=0;
  e2e25count_thresh1_nor4=0;

  r4_algo_count=0;
  e2e9_algo_count=0;
  e2e9_algo_count_5_1=0;
  e2e9_algo_count_5_0=0;

  swisscross_algo=0;
  e2e9_algo=0;

  scale=1.0;
  ncr_ = 0 ;
  
  energy_pf_     = 0;
  energyc_pf_    = 0;
  energyn_pf_    = 0;
  energyg_pf_    = 0;
  
  energy_ecal    = 0;
  energy_hcal    = 0;
  
  ptJet_      = 0;
  etaJet_     = 0;
  phiJet_     = 0;
  chfJet_     = 0;
  nhfJet_     = 0;
  cemfJet_    = 0;
  nemfJet_    = 0;
  cmultiJet_  = 0;
  nmultiJet_  = 0;
  
  
  nrjets_     = 1;

	
  eesum_gt1=0;
  eesum_gt2=0;
  eesum_gt4=0;

  ebsum_gt1=0;
  ebsum_gt2=0;
  ebsum_gt4=0;



  eesum_gt1et=0;
  eesum_gt2et=0;
  eesum_gt4et=0;

  ebsum_gt1et=0;
  ebsum_gt2et=0;
  ebsum_gt4et=0;

   
  ebhits1GeVet=0;
  ebhits2GeV=0;
  ebhits4GeV=0;
  eehits2GeV=0;
  eehits4GeV=0;

  eehits1GeVet=0;
  ebhits2GeVet=0;
  ebhits4GeVet=0;
  eehits2GeVet=0;
  eehits4GeVet=0;





  Bool_t trigbits=false;
  Bool_t physdec=false;
  Bool_t goodbx=false;
  Bool_t nomonster=false;


 
 

  run         = event.id().run();
  even        = event.id().event();
  lumi        = event.luminosityBlock();
  bx          = event.bunchCrossing();
  //  physdeclared= gtRecord->gtFdlWord().physicsDeclared();
  orbit       = event.orbitNumber(); 
  double sec  = event.time().value() >> 32 ;
  double usec = 0xFFFFFFFF & event.time().value();
  double conv = 1e6;
  time        = (sec+usec/conv);




  // get position in bunch train

  bunchintrain=-1;

  for (std::vector<int>::const_iterator bxit=bunchstartbx_.begin(); bxit!=bunchstartbx_.end(); ++bxit) {

    Int_t bxpos=bx - *bxit;


    // 50ns
    if (bxpos>=0 && bxpos<=70) bunchintrain=bxpos/2;


    // 25ns
    // if (bxpos>=0 && bxpos<=50) bunchintrain=bxpos;

  }


  Float_t toffset=0.0;
  


  Int_t numtracks=0;
  Int_t goodtracks=0;




 
  edm::Handle<VertexCollection> vertices;
  event.getByLabel(vertex_coll_,vertices);
  VertexCollection::const_iterator vit;

  numvtx=vertices->size();
  

  if (vertices->size()>0) {

    for (vit=vertices->begin();vit!=vertices->end();++vit) {

      vtx_x=vit->x();
      vtx_y=vit->y();
      vtx_z=vit->z();
      vtx_x_err=vit->xError();
      vtx_y_err=vit->yError();
      vtx_z_err=vit->zError();

      vtx_chi2=vit->chi2();
      vtx_ndof=vit->ndof();

      //      vector<TrackBaseRef>::const_iterator trstart=vit->tracks_begin();
      //      vector<TrackBaseRef>::const_iterator trend=vit->tracks_end();


      vtx_ntracks=vit->tracksSize();


      vtx_isfake=vit->isFake();

      if (vit->isValid() && vtx_isfake==0 && vtx_ndof>4 && vtx_chi2>0 && vtx_chi2<10000) {

        vtx_good=1;
        numgoodvtx++;

      }
     
    }

  } 

  
 

  physdec=true;



  trigbits=true;



 
  if (bit0) goodbx=true;

  if (isMC_) {
    goodbx=true;
    physdec=true;
  }

  Bool_t goodevent=false;

  goodevent= goodbx && physdec && trigbits && nomonster;







 


  // Rechit Collections

  Handle<EcalRecHitCollection> EBhits;
  event.getByLabel("ecalRecHit","EcalRecHitsEB",EBhits);

  //event.getByLabel("reducedEcalRecHitsEB","",EBhits);
  //event.getByLabel("alCaIsolatedElectrons","alCaRecHitsEB",EBhits);

  const EcalRecHitCollection *ebRecHits=EBhits.product();



  Handle<EcalRecHitCollection> EEhits;
  event.getByLabel("ecalRecHit","EcalRecHitsEE",EEhits);
 
  //event.getByLabel("reducedEcalRecHitsEE","",EEhits);
  // event.getByLabel("alCaIsolatedElectrons","alCaRecHitsEE",EEhits);

  const EcalRecHitCollection *eeRecHits=EEhits.product();

 

  // Digis, from RAW
 
  
  Handle<EBDigiCollection>  EBdigis;
  event.getByLabel("ecalEBunpacker","ebDigis",EBdigis);

  Handle<EEDigiCollection>  EEdigis;

  event.getByLabel("ecalEBunpacker","eeDigis",EEdigis);
  


  // use this if we don't have RAW files

  /*
  Handle<EBDigiCollection>  EBdigis;
  event.getByLabel("selectDigi","selectedEcalEBDigiCollection",EBdigis);

  Handle<EEDigiCollection>  EEdigis;

  event.getByLabel("selectDigi","selectedEcalEEDigiCollection",EEdigis);
  */
  




  // Calotower section

  edm::Handle<CaloTowerCollection> towers;
  event.getByLabel("towerMaker",towers);


  ncalotower = 0;

  ncalotowereb=0;
  ncalotoweree=0;
  ncalotowerhf=0;

  ncalotowerebgt05=0;
  ncalotowerebgt1=0;
  ncalotowerebgt2=0;
  ncalotowerebgt5=0;
  ncalotowerebgt10=0;

  ncalotowereegt05=0;
  ncalotowereegt1=0;
  ncalotowereegt2=0;
  ncalotowereegt5=0;
  ncalotowereegt10=0;

  ncalotowerhfgt05=0;
  ncalotowerhfgt1=0;
  ncalotowerhfgt2=0;
  ncalotowerhfgt5=0;
  ncalotowerhfgt10=0;


  ncalotowerebgt05had=0;
  ncalotowerebgt1had=0;
  ncalotowerebgt2had=0;
  ncalotowerebgt5had=0;
  ncalotowerebgt10had=0;

  ncalotowereegt05had=0;
  ncalotowereegt1had=0;
  ncalotowereegt2had=0;
  ncalotowereegt5had=0;
  ncalotowereegt10had=0;

  ncalotowerhfgt05had=0;
  ncalotowerhfgt1had=0;
  ncalotowerhfgt2had=0;
  ncalotowerhfgt5had=0;
  ncalotowerhfgt10had=0;

  ncalotowerebgt05hof=0;
  ncalotowerebgt1hof=0;
  ncalotowerebgt2hof=0;
  ncalotowerebgt5hof=0;
  ncalotowerebgt10hof=0;

  ncalotowereegt05hof=0;
  ncalotowereegt1hof=0;
  ncalotowereegt2hof=0;
  ncalotowereegt5hof=0;
  ncalotowereegt10hof=0;

  ncalotowerhfgt05hof=0;
  ncalotowerhfgt1hof=0;
  ncalotowerhfgt2hof=0;
  ncalotowerhfgt5hof=0;
  ncalotowerhfgt10hof=0;

  ncalotowerebgt05all=0;
  ncalotowerebgt1all=0;
  ncalotowerebgt2all=0;
  ncalotowerebgt5all=0;
  ncalotowerebgt10all=0;

  ncalotowereegt1all=0;
  ncalotowereegt2all=0;
  ncalotowereegt5all=0;
  ncalotowereegt10all=0;

  ncalotowerhfgt05all=0;
  ncalotowerhfgt1all=0;
  ncalotowerhfgt2all=0;
  ncalotowerhfgt5all=0;
  ncalotowerhfgt10all=0;



  ctsumebgt05=0;
  ctsumebgt1=0;
  ctsumebgt2=0;
  ctsumebgt5=0;
  ctsumebgt10=0;

  ctsumeegt05=0;
  ctsumeegt1=0;
  ctsumeegt2=0;
  ctsumeegt5=0;
  ctsumeegt10=0;

  ctsumhfgt05=0;
  ctsumhfgt1=0;
  ctsumhfgt2=0;
  ctsumhfgt5=0;
  ctsumhfgt10=0;


  ctsumebgt05had=0;
  ctsumebgt1had=0;
  ctsumebgt2had=0;
  ctsumebgt5had=0;
  ctsumebgt10had=0;

  ctsumeegt05had=0;
  ctsumeegt1had=0;
  ctsumeegt2had=0;
  ctsumeegt5had=0;
  ctsumeegt10had=0;

  ctsumhfgt05had=0;
  ctsumhfgt1had=0;
  ctsumhfgt2had=0;
  ctsumhfgt5had=0;
  ctsumhfgt10had=0;


  ctsumebgt05hof=0;
  ctsumebgt1hof=0;
  ctsumebgt2hof=0;
  ctsumebgt5hof=0;
  ctsumebgt10hof=0;

  ctsumeegt05hof=0;
  ctsumeegt1hof=0;
  ctsumeegt2hof=0;
  ctsumeegt5hof=0;
  ctsumeegt10hof=0;

  ctsumhfgt05hof=0;
  ctsumhfgt1hof=0;
  ctsumhfgt2hof=0;
  ctsumhfgt5hof=0;
  ctsumhfgt10hof=0;


  ctsumebgt05all=0;
  ctsumebgt1all=0;
  ctsumebgt2all=0;
  ctsumebgt5all=0;
  ctsumebgt10all=0;

  ctsumeegt05all=0;
  ctsumeegt1all=0;
  ctsumeegt2all=0;
  ctsumeegt5all=0;
  ctsumeegt10all=0;

  ctsumhfgt05all=0;
  ctsumhfgt1all=0;
  ctsumhfgt2all=0;
  ctsumhfgt5all=0;
  ctsumhfgt10all=0;




  if (towers.isValid() && towers->size()>0) {
    for(CaloTowerCollection::const_iterator it = towers->begin();it != towers->end(); it++) {
      
      double ct_eta=it->eta();
      double ct_phi = it->phi();
      double ct_et = it->et();
      double ct_emf=it->emEt();
      double ct_had=it->hadEt();
      double ct_hof=it->outerEt();


      int ct_ieta=it->ieta();
      int ct_iphi=it->iphi();

      if (ct_emf>0.5 && fabs(ct_ieta)<=20) {

       
	calotower_map->Fill(ct_iphi,ct_ieta,ct_emf);

      }

 

      if (ct_emf>1.0 && fabs(ct_ieta)<=28 && fabs(ct_ieta)>20) {

       
	calotower_map->Fill(ct_iphi,ct_ieta,ct_emf/2);
	calotower_map->Fill(ct_iphi+1,ct_ieta,ct_emf/2);

      }

      ncalotower++;



      if (fabs(ct_eta)<3.0 && ct_emf>0.0) {

	calotowereta_all->Fill(ct_eta);
	calotowereta_etweight->Fill(ct_eta,ct_emf);

	if (ct_emf>1.0) {

	  calotowereta_gt1et->Fill(ct_eta);
	  calotowereta_etweight_gt1et->Fill(ct_eta,ct_emf);


	  if (numgoodvtx<10) 	  calotowereta_gt1et_pu10->Fill(ct_eta);
	  if (numgoodvtx>=10 && numgoodvtx<20) 	  calotowereta_gt1et_pu20->Fill(ct_eta);
	  if (numgoodvtx>20) 	  calotowereta_gt1et_pu30->Fill(ct_eta);


	}	

      }


      if (fabs(ct_eta)<=1.48) {

	ncalotowereb++;

	if (ct_emf>0.5) {
	  ncalotowerebgt05++;
	  ctsumebgt05+=ct_emf;
	}

	if (ct_emf>1.0) {
	  ncalotowerebgt1++;
	  ctsumebgt1+=ct_emf;
	}

	if (ct_emf>2.0) {
	  ncalotowerebgt2++;
	  ctsumebgt2+=ct_emf;
	}

	if (ct_emf>5.0) {
	  ncalotowerebgt5++;
	  ctsumebgt5+=ct_emf;
	}

	if (ct_emf>10.0) {
	  ncalotowerebgt10++;
	  ctsumebgt10+=ct_emf;
	}


	if (ct_had>0.5) {
	  ncalotowerebgt05had++;
	  ctsumebgt05had+=ct_had;
	}

	if (ct_had>1.0) {
	  ncalotowerebgt1had++;
	  ctsumebgt1had+=ct_had;
	}

	if (ct_had>2.0) {
	  ncalotowerebgt2had++;
	  ctsumebgt2had+=ct_had;
	}

	if (ct_had>5.0) {
	  ncalotowerebgt5had++;
	  ctsumebgt5had+=ct_had;
	}

	if (ct_had>10.0) {
	  ncalotowerebgt10had++;
	  ctsumebgt10had+=ct_had;
	}


	if (ct_hof>0.5) {
	  ncalotowerebgt05hof++;
	  ctsumebgt05hof+=ct_hof;
	}

	if (ct_hof>1.0) {
	  ncalotowerebgt1hof++;
	  ctsumebgt1hof+=ct_hof;
	}

	if (ct_hof>2.0) {
	  ncalotowerebgt2hof++;
	  ctsumebgt2hof+=ct_hof;
	}

	if (ct_hof>5.0) {
	  ncalotowerebgt5hof++;
	  ctsumebgt5hof+=ct_hof;
	}

	if (ct_hof>10.0) {
	  ncalotowerebgt10hof++;
	  ctsumebgt10hof+=ct_hof;
	}


	if (ct_et>0.5) {
	  ncalotowerebgt05all++;
	  ctsumebgt05all+=ct_et;
	}

	if (ct_et>1.0) {
	  ncalotowerebgt1all++;
	  ctsumebgt1all+=ct_et;
	}

	if (ct_et>2.0) {
	  ncalotowerebgt2all++;
	  ctsumebgt2all+=ct_et;
	}

	if (ct_et>5.0) {
	  ncalotowerebgt5all++;
	  ctsumebgt5all+=ct_et;
	}

	if (ct_et>10.0) {
	  ncalotowerebgt10all++;
	  ctsumebgt10all+=ct_et;
	}



      }

      if (fabs(ct_eta)>1.48 && fabs(ct_eta)<3.0) {


	ncalotoweree++;

	if (ct_emf>0.5) {
	  ncalotowereegt05++;
	  ctsumeegt05+=ct_emf;
	}

	if (ct_emf>1.0) {
	  ncalotowereegt1++;
	  ctsumeegt1+=ct_emf;
	}

	if (ct_emf>2.0) {
	  ncalotowereegt2++;
	  ctsumeegt2+=ct_emf;
	}

	if (ct_emf>5.0) {
	  ncalotowereegt5++;
	  ctsumeegt5+=ct_emf;
	}

	if (ct_emf>10.0) {
	  ncalotowereegt10++;
	  ctsumeegt10+=ct_emf;
	}


	if (ct_had>0.5) {
	  ncalotowereegt05had++;
	  ctsumeegt05had+=ct_had;
	}


	if (ct_had>1.0) {
	  ncalotowereegt1had++;
	  ctsumeegt1had+=ct_had;
	}

	if (ct_had>2.0) {
	  ncalotowereegt2had++;
	  ctsumeegt2had+=ct_had;
	}

	if (ct_had>5.0) {
	  ncalotowereegt5had++;
	  ctsumeegt5had+=ct_had;
	}

	if (ct_had>10.0) {
	  ncalotowereegt10had++;
	  ctsumeegt10had+=ct_had;
	}


	if (ct_hof>0.5) {
	  ncalotowereegt05hof++;
	  ctsumeegt05hof+=ct_hof;
	}

	if (ct_hof>1.0) {
	  ncalotowereegt1hof++;
	  ctsumeegt1hof+=ct_hof;
	}

	if (ct_hof>2.0) {
	  ncalotowereegt2hof++;
	  ctsumeegt2hof+=ct_hof;
	}

	if (ct_hof>5.0) {
	  ncalotowereegt5hof++;
	  ctsumeegt5hof+=ct_hof;
	}

	if (ct_hof>10.0) {
	  ncalotowereegt10hof++;
	  ctsumeegt10hof+=ct_hof;
	}


	if (ct_et>0.5) {
	  ncalotowereegt05all++;
	  ctsumeegt05all+=ct_et;
	}

	if (ct_et>1.0) {
	  ncalotowereegt1all++;
	  ctsumeegt1all+=ct_et;
	}

	if (ct_et>2.0) {
	  ncalotowereegt2all++;
	  ctsumeegt2all+=ct_et;
	}

	if (ct_et>5.0) {
	  ncalotowereegt5all++;
	  ctsumeegt5all+=ct_et;
	}

	if (ct_et>10.0) {
	  ncalotowereegt10all++;
	  ctsumeegt10all+=ct_et;
	}



      }


      if (fabs(ct_eta)>3.0) {


	ncalotowerhf++;

	if (ct_emf>0.5) {
	  ncalotowerhfgt05++;
	  ctsumhfgt05+=ct_emf;
	}

	if (ct_emf>1.0) {
	  ncalotowerhfgt1++;
	  ctsumhfgt1+=ct_emf;
	}

	if (ct_emf>2.0) {
	  ncalotowerhfgt2++;
	  ctsumhfgt2+=ct_emf;
	}

	if (ct_emf>5.0) {
	  ncalotowerhfgt5++;
	  ctsumhfgt5+=ct_emf;
	}

	if (ct_emf>10.0) {
	  ncalotowerhfgt10++;
	  ctsumhfgt10+=ct_emf;
	}



	if (ct_had>0.5) {
	  ncalotowerhfgt05had++;
	  ctsumhfgt05had+=ct_had;
	}

	if (ct_had>1.0) {
	  ncalotowerhfgt1had++;
	  ctsumhfgt1had+=ct_had;
	}

	if (ct_had>2.0) {
	  ncalotowerhfgt2had++;
	  ctsumhfgt2had+=ct_had;
	}

	if (ct_had>5.0) {
	  ncalotowerhfgt5had++;
	  ctsumhfgt5had+=ct_had;
	}

	if (ct_had>10.0) {
	  ncalotowerhfgt10had++;
	  ctsumhfgt10had+=ct_had;
	}


	if (ct_hof>0.5) {
	  ncalotowerhfgt05hof++;
	  ctsumhfgt05hof+=ct_hof;
	}

	if (ct_hof>1.0) {
	  ncalotowerhfgt1hof++;
	  ctsumhfgt1hof+=ct_hof;
	}

	if (ct_hof>2.0) {
	  ncalotowerhfgt2hof++;
	  ctsumhfgt2hof+=ct_hof;
	}

	if (ct_hof>5.0) {
	  ncalotowerhfgt5hof++;
	  ctsumhfgt5hof+=ct_hof;
	}

	if (ct_hof>10.0) {
	  ncalotowerhfgt10hof++;
	  ctsumhfgt10hof+=ct_hof;
	}


	if (ct_et>0.5) {
	  ncalotowerhfgt05all++;
	  ctsumhfgt05all+=ct_et;
	}

	if (ct_et>1.0) {
	  ncalotowerhfgt1all++;
	  ctsumhfgt1all+=ct_et;
	}

	if (ct_et>2.0) {
	  ncalotowerhfgt2all++;
	  ctsumhfgt2all+=ct_et;
	}

	if (ct_et>5.0) {
	  ncalotowerhfgt5all++;
	  ctsumhfgt5all+=ct_et;
	}

	if (ct_et>10.0) {
	  ncalotowerhfgt10all++;
	  ctsumhfgt10all+=ct_et;
	}

      }




    }

  }

  edm::ESHandle<CaloGeometry> pG;
  iSetup.get<CaloGeometryRecord>().get(pG);
  const CaloGeometry* geo=pG.product();

  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);


  edm::ESHandle<EcalChannelStatus> chanstat;
  iSetup.get<EcalChannelStatusRcd>().get(chanstat);
  const EcalChannelStatus* cstat=chanstat.product();
 
  edm::ESHandle<EcalPedestals> ecalped;
  iSetup.get<EcalPedestalsRcd>().get(ecalped);
  const EcalPedestals* eped=ecalped.product();
 




  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  const EcalSeverityLevelAlgo *severity=sevlv.product();

  edm::Handle<reco::SuperClusterCollection> sc_eb;
  event.getByLabel("correctedHybridSuperClusters","",sc_eb);


  edm::Handle<reco::SuperClusterCollection> sc_ee;
  event.getByLabel("correctedMulti5x5SuperClustersWithPreshower","",sc_ee);



  // Supercluster loop:  EB


  ebscsumet_all=0;
  float ebscsumet_severity0=0;
  float ebscsumet_koot0=0;

  float ebscsumet_all_gt2=0;
  float ebscsumet_severity0_gt2=0;
  float ebscsumet_koot0_gt2=0;

  float ebscsumet_all_gt5=0;
  float ebscsumet_severity0_gt5=0;
  float ebscsumet_koot0_gt5=0;

  float ebscsumet_all_gt10=0;
  float ebscsumet_severity0_gt10=0;
  float ebscsumet_koot0_gt10=0;


  ebscsumet_all_eta15=0;
  ebscsumet_all_eta20=0;
  ebscsumet_all_eta25=0;


  float ebscsumet_all_eta15_pu10=0;
  float ebscsumet_all_eta20_pu10=0;
  float ebscsumet_all_eta25_pu10=0;

  float ebscsumet_all_eta15_pu20=0;
  float ebscsumet_all_eta20_pu20=0;
  float ebscsumet_all_eta25_pu20=0;

  float ebscsumet_all_eta15_pu30=0;
  float ebscsumet_all_eta20_pu30=0;
  float ebscsumet_all_eta25_pu30=0;

  float ebscsumet_all_eta15_pueq10=0;
  float ebscsumet_all_eta20_pueq10=0;
  float ebscsumet_all_eta25_pueq10=0;

  float ebscsumet_all_eta15_pueq20=0;
  float ebscsumet_all_eta20_pueq20=0;
  float ebscsumet_all_eta25_pueq20=0;



  eescsumet_all_eta15=0;
  eescsumet_all_eta20=0;
  eescsumet_all_eta25=0;


  float eescsumet_all_eta15_pu10=0;
  float eescsumet_all_eta20_pu10=0;
  float eescsumet_all_eta25_pu10=0;

  float eescsumet_all_eta15_pu20=0;
  float eescsumet_all_eta20_pu20=0;
  float eescsumet_all_eta25_pu20=0;

  float eescsumet_all_eta15_pu30=0;
  float eescsumet_all_eta20_pu30=0;
  float eescsumet_all_eta25_pu30=0;

  float eescsumet_all_eta15_pueq10=0;
  float eescsumet_all_eta20_pueq10=0;
  float eescsumet_all_eta25_pueq10=0;

  float eescsumet_all_eta15_pueq20=0;
  float eescsumet_all_eta20_pueq20=0;
  float eescsumet_all_eta25_pueq20=0;


  float ebscsumet_all_eta15_pueq01=0;
  float ebscsumet_all_eta20_pueq01=0;
  float ebscsumet_all_eta25_pueq01=0;
  float eescsumet_all_eta15_pueq01=0;
  float eescsumet_all_eta20_pueq01=0;
  float eescsumet_all_eta25_pueq01=0;

  float ebscsumet_all_eta15_pueq02=0;
  float ebscsumet_all_eta20_pueq02=0;
  float ebscsumet_all_eta25_pueq02=0;
  float eescsumet_all_eta15_pueq02=0;
  float eescsumet_all_eta20_pueq02=0;
  float eescsumet_all_eta25_pueq02=0;

  float ebscsumet_all_eta15_pueq03=0;
  float ebscsumet_all_eta20_pueq03=0;
  float ebscsumet_all_eta25_pueq03=0;
  float eescsumet_all_eta15_pueq03=0;
  float eescsumet_all_eta20_pueq03=0;
  float eescsumet_all_eta25_pueq03=0;

  float ebscsumet_all_eta15_pueq04=0;
  float ebscsumet_all_eta20_pueq04=0;
  float ebscsumet_all_eta25_pueq04=0;
  float eescsumet_all_eta15_pueq04=0;
  float eescsumet_all_eta20_pueq04=0;
  float eescsumet_all_eta25_pueq04=0;

  float ebscsumet_all_eta15_pueq05=0;
  float ebscsumet_all_eta20_pueq05=0;
  float ebscsumet_all_eta25_pueq05=0;
  float eescsumet_all_eta15_pueq05=0;
  float eescsumet_all_eta20_pueq05=0;
  float eescsumet_all_eta25_pueq05=0;

  float ebscsumet_all_eta15_pueq06=0;
  float ebscsumet_all_eta20_pueq06=0;
  float ebscsumet_all_eta25_pueq06=0;
  float eescsumet_all_eta15_pueq06=0;
  float eescsumet_all_eta20_pueq06=0;
  float eescsumet_all_eta25_pueq06=0;

  float ebscsumet_all_eta15_pueq07=0;
  float ebscsumet_all_eta20_pueq07=0;
  float ebscsumet_all_eta25_pueq07=0;
  float eescsumet_all_eta15_pueq07=0;
  float eescsumet_all_eta20_pueq07=0;
  float eescsumet_all_eta25_pueq07=0;

  float ebscsumet_all_eta15_pueq08=0;
  float ebscsumet_all_eta20_pueq08=0;
  float ebscsumet_all_eta25_pueq08=0;
  float eescsumet_all_eta15_pueq08=0;
  float eescsumet_all_eta20_pueq08=0;
  float eescsumet_all_eta25_pueq08=0;

  float ebscsumet_all_eta15_pueq09=0;
  float ebscsumet_all_eta20_pueq09=0;
  float ebscsumet_all_eta25_pueq09=0;
  float eescsumet_all_eta15_pueq09=0;
  float eescsumet_all_eta20_pueq09=0;
  float eescsumet_all_eta25_pueq09=0;

  ebnumsc_all=0;
  eenumsc_all=0;

  for (reco::SuperClusterCollection::const_iterator sCEB = sc_eb->begin(); sCEB != sc_eb->end(); ++sCEB) {

    DetId seedid=sCEB->seed()->seed();

    bool goodseed=false;

    int sevlev=0;
    float time=0;
    int koot=0;

    for (EcalRecHitCollection::const_iterator hitItr = EBhits->begin(); hitItr != EBhits->end(); ++hitItr) {

      EcalRecHit hit = (*hitItr);
      EBDetId det    = hit.id(); 
    
      if (det!=seedid) continue;

      sevlev=severity->severityLevel(det,*ebRecHits);
      time     = hit.time();

      koot=0;
      if (hit.checkFlag(EcalRecHit::kOutOfTime)) koot=1;

      goodseed=true;

    }

    if (goodseed) {

      Float_t scenergy=sCEB->energy();
      Float_t sc_eta=sCEB->eta();
      Float_t sc_phi=sCEB->phi();
      Float_t scet=scenergy/cosh(sc_eta);

      ebscsumet_all+=scet;
      if (sevlev==0) ebscsumet_severity0+=scet;
      if (koot==0) ebscsumet_koot0+=scet;
      ebnumsc_all++;
    


      scet_eb_all->Fill(scet);
      if (sevlev==0) scet_eb_severity0->Fill(scet);
      if (koot==0) scet_eb_koot0->Fill(scet);
      


      sceta_all->Fill(sc_eta);
 
      sceta_vs_bxtrain->Fill(bunchintrain,sc_eta);

      if (sevlev==0) sceta_severity0->Fill(sc_eta);
      if (koot==0) sceta_koot0->Fill(sc_eta);

 
      if (numgoodvtx==1) {
	sceta_all_pueq01->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq01->Fill(sc_eta);
      }

 
      if (numgoodvtx==2) {
	sceta_all_pueq02->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq02->Fill(sc_eta);
      }

 
      if (numgoodvtx==3) {
	sceta_all_pueq03->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq03->Fill(sc_eta);
      }

 
      if (numgoodvtx==4) {
	sceta_all_pueq04->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq04->Fill(sc_eta);
      }

 
      if (numgoodvtx==5) {
	sceta_all_pueq05->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq05->Fill(sc_eta);
      }

 
      if (numgoodvtx==6) {
	sceta_all_pueq06->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq06->Fill(sc_eta);
      }

 
      if (numgoodvtx==7) {
	sceta_all_pueq07->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq07->Fill(sc_eta);
      }

 
      if (numgoodvtx==8) {
	sceta_all_pueq08->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq08->Fill(sc_eta);
      }

 
      if (numgoodvtx==9) {
	sceta_all_pueq09->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq09->Fill(sc_eta);
      }




      if (fabs(sc_eta)<0.5) {

	scet_eb_all_eta15->Fill(scet);	
	ebscsumet_all_eta15+=scet;


	if (numgoodvtx<=10) {
	  scet_eb_all_eta15_pu10->Fill(scet);
	  ebscsumet_all_eta15_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_eb_all_eta15_pu20->Fill(scet);
	  ebscsumet_all_eta15_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_eb_all_eta15_pu30->Fill(scet);
	  ebscsumet_all_eta15_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_eb_all_eta15_pueq10->Fill(scet);
	  ebscsumet_all_eta15_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_eb_all_eta15_pueq20->Fill(scet);
	  ebscsumet_all_eta15_pueq20+=scet;
	}

	if (numgoodvtx==1) {
	  scet_eb_all_eta15_pueq01->Fill(scet);
	  ebscsumet_all_eta15_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_eb_all_eta15_pueq02->Fill(scet);
	  ebscsumet_all_eta15_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_eb_all_eta15_pueq03->Fill(scet);
	  ebscsumet_all_eta15_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_eb_all_eta15_pueq04->Fill(scet);
	  ebscsumet_all_eta15_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_eb_all_eta15_pueq05->Fill(scet);
	  ebscsumet_all_eta15_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_eb_all_eta15_pueq06->Fill(scet);
	  ebscsumet_all_eta15_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_eb_all_eta15_pueq07->Fill(scet);
	  ebscsumet_all_eta15_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_eb_all_eta15_pueq08->Fill(scet);
	  ebscsumet_all_eta15_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_eb_all_eta15_pueq09->Fill(scet);
	  ebscsumet_all_eta15_pueq09+=scet;
	}




      }


      if (fabs(sc_eta)>=0.5 && fabs(sc_eta)<1.0) {

	scet_eb_all_eta20->Fill(scet);
	ebscsumet_all_eta20+=scet;



	if (numgoodvtx<=10) {
	  scet_eb_all_eta20_pu10->Fill(scet);
	  ebscsumet_all_eta20_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_eb_all_eta20_pu20->Fill(scet);
	  ebscsumet_all_eta20_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_eb_all_eta20_pu30->Fill(scet);
	  ebscsumet_all_eta20_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_eb_all_eta20_pueq10->Fill(scet);
	  ebscsumet_all_eta20_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_eb_all_eta20_pueq20->Fill(scet);
	  ebscsumet_all_eta20_pueq20+=scet;
	}



	if (numgoodvtx==1) {
	  scet_eb_all_eta20_pueq01->Fill(scet);
	  ebscsumet_all_eta20_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_eb_all_eta20_pueq02->Fill(scet);
	  ebscsumet_all_eta20_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_eb_all_eta20_pueq03->Fill(scet);
	  ebscsumet_all_eta20_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_eb_all_eta20_pueq04->Fill(scet);
	  ebscsumet_all_eta20_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_eb_all_eta20_pueq05->Fill(scet);
	  ebscsumet_all_eta20_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_eb_all_eta20_pueq06->Fill(scet);
	  ebscsumet_all_eta20_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_eb_all_eta20_pueq07->Fill(scet);
	  ebscsumet_all_eta20_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_eb_all_eta20_pueq08->Fill(scet);
	  ebscsumet_all_eta20_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_eb_all_eta20_pueq09->Fill(scet);
	  ebscsumet_all_eta20_pueq09+=scet;
	}



	
      }


      if (fabs(sc_eta)>1.0) {

	scet_eb_all_eta25->Fill(scet);
	ebscsumet_all_eta25+=scet;


	if (numgoodvtx<=10) {
	  scet_eb_all_eta25_pu10->Fill(scet);
	  ebscsumet_all_eta25_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_eb_all_eta25_pu20->Fill(scet);
	  ebscsumet_all_eta25_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_eb_all_eta25_pu30->Fill(scet);
	  ebscsumet_all_eta25_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_eb_all_eta25_pueq10->Fill(scet);
	  ebscsumet_all_eta25_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_eb_all_eta25_pueq20->Fill(scet);
	  ebscsumet_all_eta25_pueq20+=scet;
	}


 	if (numgoodvtx==1) {
	  scet_eb_all_eta25_pueq01->Fill(scet);
	  ebscsumet_all_eta25_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_eb_all_eta25_pueq02->Fill(scet);
	  ebscsumet_all_eta25_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_eb_all_eta25_pueq03->Fill(scet);
	  ebscsumet_all_eta25_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_eb_all_eta25_pueq04->Fill(scet);
	  ebscsumet_all_eta25_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_eb_all_eta25_pueq05->Fill(scet);
	  ebscsumet_all_eta25_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_eb_all_eta25_pueq06->Fill(scet);
	  ebscsumet_all_eta25_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_eb_all_eta25_pueq07->Fill(scet);
	  ebscsumet_all_eta25_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_eb_all_eta25_pueq08->Fill(scet);
	  ebscsumet_all_eta25_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_eb_all_eta25_pueq09->Fill(scet);
	  ebscsumet_all_eta25_pueq09+=scet;
	}




     }





      if (scet>2.0) {

	sceta_all_gt2->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt2->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt2->Fill(sc_eta);

	ebscsumet_all_gt2+=scet;
	if (sevlev==0) ebscsumet_severity0_gt2+=scet;
	if (koot==0) ebscsumet_koot0_gt2+=scet;
   

      }

     if (scet>5.0) {

	sceta_all_gt5->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt5->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt5->Fill(sc_eta);

 	ebscsumet_all_gt5+=scet;
	if (sevlev==0) ebscsumet_severity0_gt5+=scet;
	if (koot==0) ebscsumet_koot0_gt5+=scet;
 
     }

     if (scet>10.0) {

	sceta_all_gt10->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt10->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt10->Fill(sc_eta);

	ebscsumet_all_gt10+=scet;
	if (sevlev==0) ebscsumet_severity0_gt10+=scet;
	if (koot==0) ebscsumet_koot0_gt10+=scet;
 
      }

     if (scet>50.0) {
       scocc_eb_gt50->Fill(sc_phi,sc_eta);
     }
     
    }
  }


  scsumet_eb_all->Fill(ebscsumet_all);
  scsumet_eb_severity0->Fill(ebscsumet_severity0);
  scsumet_eb_koot0->Fill(ebscsumet_koot0);

  scsumet_eb_all_eta15->Fill(ebscsumet_all_eta15);
  scsumet_eb_all_eta20->Fill(ebscsumet_all_eta20);
  scsumet_eb_all_eta25->Fill(ebscsumet_all_eta25);



  scsumet_eb_all_eta15_pu10->Fill(ebscsumet_all_eta15_pu10);
  scsumet_eb_all_eta20_pu10->Fill(ebscsumet_all_eta20_pu10);
  scsumet_eb_all_eta25_pu10->Fill(ebscsumet_all_eta25_pu10);

  scsumet_eb_all_eta15_pu20->Fill(ebscsumet_all_eta15_pu20);
  scsumet_eb_all_eta20_pu20->Fill(ebscsumet_all_eta20_pu20);
  scsumet_eb_all_eta25_pu20->Fill(ebscsumet_all_eta25_pu20);

  scsumet_eb_all_eta15_pu30->Fill(ebscsumet_all_eta15_pu30);
  scsumet_eb_all_eta20_pu30->Fill(ebscsumet_all_eta20_pu30);
  scsumet_eb_all_eta25_pu30->Fill(ebscsumet_all_eta25_pu30);

  scsumet_eb_all_eta15_pueq10->Fill(ebscsumet_all_eta15_pueq10);
  scsumet_eb_all_eta20_pueq10->Fill(ebscsumet_all_eta20_pueq10);
  scsumet_eb_all_eta25_pueq10->Fill(ebscsumet_all_eta25_pueq10);

  scsumet_eb_all_eta15_pueq20->Fill(ebscsumet_all_eta15_pueq20);
  scsumet_eb_all_eta20_pueq20->Fill(ebscsumet_all_eta20_pueq20);
  scsumet_eb_all_eta25_pueq20->Fill(ebscsumet_all_eta25_pueq20);


  scsumet_eb_all_eta15_pueq01->Fill(ebscsumet_all_eta15_pueq01);
  scsumet_eb_all_eta20_pueq01->Fill(ebscsumet_all_eta20_pueq01);
  scsumet_eb_all_eta25_pueq01->Fill(ebscsumet_all_eta25_pueq01);

  scsumet_eb_all_eta15_pueq02->Fill(ebscsumet_all_eta15_pueq02);
  scsumet_eb_all_eta20_pueq02->Fill(ebscsumet_all_eta20_pueq02);
  scsumet_eb_all_eta25_pueq02->Fill(ebscsumet_all_eta25_pueq02);

  scsumet_eb_all_eta15_pueq03->Fill(ebscsumet_all_eta15_pueq03);
  scsumet_eb_all_eta20_pueq03->Fill(ebscsumet_all_eta20_pueq03);
  scsumet_eb_all_eta25_pueq03->Fill(ebscsumet_all_eta25_pueq03);

  scsumet_eb_all_eta15_pueq04->Fill(ebscsumet_all_eta15_pueq04);
  scsumet_eb_all_eta20_pueq04->Fill(ebscsumet_all_eta20_pueq04);
  scsumet_eb_all_eta25_pueq04->Fill(ebscsumet_all_eta25_pueq04);

  scsumet_eb_all_eta15_pueq05->Fill(ebscsumet_all_eta15_pueq05);
  scsumet_eb_all_eta20_pueq05->Fill(ebscsumet_all_eta20_pueq05);
  scsumet_eb_all_eta25_pueq05->Fill(ebscsumet_all_eta25_pueq05);

  scsumet_eb_all_eta15_pueq06->Fill(ebscsumet_all_eta15_pueq06);
  scsumet_eb_all_eta20_pueq06->Fill(ebscsumet_all_eta20_pueq06);
  scsumet_eb_all_eta25_pueq06->Fill(ebscsumet_all_eta25_pueq06);

  scsumet_eb_all_eta15_pueq07->Fill(ebscsumet_all_eta15_pueq07);
  scsumet_eb_all_eta20_pueq07->Fill(ebscsumet_all_eta20_pueq07);
  scsumet_eb_all_eta25_pueq07->Fill(ebscsumet_all_eta25_pueq07);

  scsumet_eb_all_eta15_pueq08->Fill(ebscsumet_all_eta15_pueq08);
  scsumet_eb_all_eta20_pueq08->Fill(ebscsumet_all_eta20_pueq08);
  scsumet_eb_all_eta25_pueq08->Fill(ebscsumet_all_eta25_pueq08);

  scsumet_eb_all_eta15_pueq09->Fill(ebscsumet_all_eta15_pueq09);
  scsumet_eb_all_eta20_pueq09->Fill(ebscsumet_all_eta20_pueq09);
  scsumet_eb_all_eta25_pueq09->Fill(ebscsumet_all_eta25_pueq09);


  scsumet_eb_all_gt2->Fill(ebscsumet_all_gt2);
  scsumet_eb_severity0_gt2->Fill(ebscsumet_severity0_gt2);
  scsumet_eb_koot0_gt2->Fill(ebscsumet_koot0_gt2);

  scsumet_eb_all_gt5->Fill(ebscsumet_all_gt5);
  scsumet_eb_severity0_gt5->Fill(ebscsumet_severity0_gt5);
  scsumet_eb_koot0_gt5->Fill(ebscsumet_koot0_gt5);

  scsumet_eb_all_gt10->Fill(ebscsumet_all_gt10);
  scsumet_eb_severity0_gt10->Fill(ebscsumet_severity0_gt10);
  scsumet_eb_koot0_gt10->Fill(ebscsumet_koot0_gt10);


  // Supercluster loop:  EE


  eescsumet_all=0;
  float eescsumet_severity0=0;
  float eescsumet_koot0=0;

  float eescsumet_all_gt2=0;
  float eescsumet_severity0_gt2=0;
  float eescsumet_koot0_gt2=0;

  float eescsumet_all_gt5=0;
  float eescsumet_severity0_gt5=0;
  float eescsumet_koot0_gt5=0;

  float eescsumet_all_gt10=0;
  float eescsumet_severity0_gt10=0;
  float eescsumet_koot0_gt10=0;




  for (reco::SuperClusterCollection::const_iterator sCEE = sc_ee->begin(); sCEE != sc_ee->end(); ++sCEE) {

    // cout << "IN EE SC" << endl;
    DetId seedid=sCEE->seed()->seed();

    bool goodseed=false;

    int sevlev=0;
    float time=0;
    int koot=0;

    float seeden=0;


    for (EcalRecHitCollection::const_iterator hitItr = EEhits->begin(); hitItr != EEhits->end(); ++hitItr) {

      EcalRecHit hit = (*hitItr);
      EEDetId det    = hit.id(); 
  
      if (det!=seedid) continue;

      sevlev=severity->severityLevel(det,*eeRecHits);
      time     = hit.time();

      seeden = hit.energy();

      koot=0;

      if (hit.checkFlag(EcalRecHit::kOutOfTime)) koot=1;
      // cout << "sevlev=" << sevlev << endl;

      //      if (sevlev!=0) goodseed=false;
      
      goodseed=true;

    }

    if (goodseed) {

      Float_t scenergy=sCEE->energy();
      Float_t sc_eta=sCEE->eta();
      Float_t sc_phi=sCEE->phi();
      Float_t scet=scenergy/cosh(sc_eta);

      sceta_all->Fill(sc_eta);

      sceta_vs_bxtrain->Fill(bunchintrain,sc_eta);


      if (sevlev==0) sceta_severity0->Fill(sc_eta);
      if (koot==0) sceta_koot0->Fill(sc_eta);

      eenumsc_all++;

      if (numgoodvtx==1) {
	sceta_all_pueq01->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq01->Fill(sc_eta);
      }

 
      if (numgoodvtx==2) {
	sceta_all_pueq02->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq02->Fill(sc_eta);
      }

 
      if (numgoodvtx==3) {
	sceta_all_pueq03->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq03->Fill(sc_eta);
      }

 
      if (numgoodvtx==4) {
	sceta_all_pueq04->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq04->Fill(sc_eta);
      }

 
      if (numgoodvtx==5) {
	sceta_all_pueq05->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq05->Fill(sc_eta);
      }

 
      if (numgoodvtx==6) {
	sceta_all_pueq06->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq06->Fill(sc_eta);
      }

 
      if (numgoodvtx==7) {
	sceta_all_pueq07->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq07->Fill(sc_eta);
      }

 
      if (numgoodvtx==8) {
	sceta_all_pueq08->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq08->Fill(sc_eta);
      }

 
      if (numgoodvtx==9) {
	sceta_all_pueq09->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_pueq09->Fill(sc_eta);
      }





      scet_ee_all->Fill(scet);
      if (sevlev==0) scet_ee_severity0->Fill(scet);
      if (koot==0) scet_ee_koot0->Fill(scet);
      

	
      eescsumet_all+=scet;
      if (sevlev==0) eescsumet_severity0+=scet;
      if (koot==0) eescsumet_koot0+=scet;



      if (fabs(sc_eta)<2.0) {

	scet_ee_all_eta15->Fill(scet);	
	eescsumet_all_eta15+=scet;


	if (numgoodvtx<=10) {
	  scet_ee_all_eta15_pu10->Fill(scet);
	  eescsumet_all_eta15_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_ee_all_eta15_pu20->Fill(scet);
	  eescsumet_all_eta15_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_ee_all_eta15_pu30->Fill(scet);
	  eescsumet_all_eta15_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_ee_all_eta15_pueq10->Fill(scet);
	  eescsumet_all_eta15_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_ee_all_eta15_pueq20->Fill(scet);
	  eescsumet_all_eta15_pueq20+=scet;
	}



	if (numgoodvtx==1) {
	  scet_ee_all_eta15_pueq01->Fill(scet);
	  eescsumet_all_eta15_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_ee_all_eta15_pueq02->Fill(scet);
	  eescsumet_all_eta15_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_ee_all_eta15_pueq03->Fill(scet);
	  eescsumet_all_eta15_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_ee_all_eta15_pueq04->Fill(scet);
	  eescsumet_all_eta15_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_ee_all_eta15_pueq05->Fill(scet);
	  eescsumet_all_eta15_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_ee_all_eta15_pueq06->Fill(scet);
	  eescsumet_all_eta15_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_ee_all_eta15_pueq07->Fill(scet);
	  eescsumet_all_eta15_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_ee_all_eta15_pueq08->Fill(scet);
	  eescsumet_all_eta15_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_ee_all_eta15_pueq09->Fill(scet);
	  eescsumet_all_eta15_pueq09+=scet;
	}



      }


      if (fabs(sc_eta)>=2.0 && fabs(sc_eta)<2.5) {

	scet_ee_all_eta20->Fill(scet);
	eescsumet_all_eta20+=scet;



	if (numgoodvtx<=10) {
	  scet_ee_all_eta20_pu10->Fill(scet);
	  eescsumet_all_eta20_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_ee_all_eta20_pu20->Fill(scet);
	  eescsumet_all_eta20_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_ee_all_eta20_pu30->Fill(scet);
	  eescsumet_all_eta20_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_ee_all_eta20_pueq10->Fill(scet);
	  eescsumet_all_eta20_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_ee_all_eta20_pueq20->Fill(scet);
	  eescsumet_all_eta20_pueq20+=scet;
	}




	if (numgoodvtx==1) {
	  scet_ee_all_eta20_pueq01->Fill(scet);
	  eescsumet_all_eta20_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_ee_all_eta20_pueq02->Fill(scet);
	  eescsumet_all_eta20_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_ee_all_eta20_pueq03->Fill(scet);
	  eescsumet_all_eta20_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_ee_all_eta20_pueq04->Fill(scet);
	  eescsumet_all_eta20_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_ee_all_eta20_pueq05->Fill(scet);
	  eescsumet_all_eta20_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_ee_all_eta20_pueq06->Fill(scet);
	  eescsumet_all_eta20_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_ee_all_eta20_pueq07->Fill(scet);
	  eescsumet_all_eta20_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_ee_all_eta20_pueq08->Fill(scet);
	  eescsumet_all_eta20_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_ee_all_eta20_pueq09->Fill(scet);
	  eescsumet_all_eta20_pueq09+=scet;
	}



	
      }


      if (fabs(sc_eta)>2.5) {

	scet_ee_all_eta25->Fill(scet);
	eescsumet_all_eta25+=scet;


	if (numgoodvtx<=10) {
	  scet_ee_all_eta25_pu10->Fill(scet);
	  eescsumet_all_eta25_pu10+=scet;
	}

	if (numgoodvtx>10 && numgoodvtx<=20) {
	  scet_ee_all_eta25_pu20->Fill(scet);
	  eescsumet_all_eta25_pu20+=scet;
	}

	if (numgoodvtx>20) {
	  scet_ee_all_eta25_pu30->Fill(scet);
	  eescsumet_all_eta25_pu30+=scet;
	}


	if (numgoodvtx==10) {
	  scet_ee_all_eta25_pueq10->Fill(scet);
	  eescsumet_all_eta25_pueq10+=scet;
	}


	if (numgoodvtx==20) {
	  scet_ee_all_eta25_pueq20->Fill(scet);
	  eescsumet_all_eta25_pueq20+=scet;
	}



	if (numgoodvtx==1) {
	  scet_ee_all_eta25_pueq01->Fill(scet);
	  eescsumet_all_eta25_pueq01+=scet;
	}

	if (numgoodvtx==2) {
	  scet_ee_all_eta25_pueq02->Fill(scet);
	  eescsumet_all_eta25_pueq02+=scet;
	}

	if (numgoodvtx==3) {
	  scet_ee_all_eta25_pueq03->Fill(scet);
	  eescsumet_all_eta25_pueq03+=scet;
	}

	if (numgoodvtx==4) {
	  scet_ee_all_eta25_pueq04->Fill(scet);
	  eescsumet_all_eta25_pueq04+=scet;
	}

	if (numgoodvtx==5) {
	  scet_ee_all_eta25_pueq05->Fill(scet);
	  eescsumet_all_eta25_pueq05+=scet;
	}

	if (numgoodvtx==6) {
	  scet_ee_all_eta25_pueq06->Fill(scet);
	  eescsumet_all_eta25_pueq06+=scet;
	}

	if (numgoodvtx==7) {
	  scet_ee_all_eta25_pueq07->Fill(scet);
	  eescsumet_all_eta25_pueq07+=scet;
	}

	if (numgoodvtx==8) {
	  scet_ee_all_eta25_pueq08->Fill(scet);
	  eescsumet_all_eta25_pueq08+=scet;
	}

	if (numgoodvtx==9) {
	  scet_ee_all_eta25_pueq09->Fill(scet);
	  eescsumet_all_eta25_pueq09+=scet;
	}



      }
	

      if (scet>2.0) {

	sceta_all_gt2->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt2->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt2->Fill(sc_eta);
	
	eescsumet_all_gt2+=scet;
	if (sevlev==0) eescsumet_severity0_gt2+=scet;
	if (koot==0) eescsumet_koot0_gt2+=scet;
   
      }

      if (scet>5.0) {

	sceta_all_gt5->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt5->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt5->Fill(sc_eta);
	
	eescsumet_all_gt5+=scet;
	if (sevlev==0) eescsumet_severity0_gt5+=scet;
	if (koot==0) eescsumet_koot0_gt5+=scet;
 
      }
      
      if (scet>10.0) {

	sceta_all_gt10->Fill(sc_eta);
	if (sevlev==0) sceta_severity0_gt10->Fill(sc_eta);
	if (koot==0) sceta_koot0_gt10->Fill(sc_eta);
	
	eescsumet_all_gt10+=scet;
	if (sevlev==0) eescsumet_severity0_gt10+=scet;
	if (koot==0) eescsumet_koot0_gt10+=scet;
 
      }


     if (scet>50.0) {
       scocc_ee_gt50->Fill(sc_phi,sc_eta);
     }
 


      //      cout << "EE sc ene=" << scenergy << " sc et=" << scet << " sc_seeden=" << seeden << " sc eta=" << sc_eta << " sc_phi=" << sc_phi  << " time=" << time << " koot=" << koot << " sevlev=" << sevlev << endl;
	
    }


  }




  scsumet_ee_all->Fill(eescsumet_all);
  scsumet_ee_severity0->Fill(eescsumet_severity0);
  scsumet_ee_koot0->Fill(eescsumet_koot0);


  scsumet_ee_all_eta15->Fill(eescsumet_all_eta15);
  scsumet_ee_all_eta20->Fill(eescsumet_all_eta20);
  scsumet_ee_all_eta25->Fill(eescsumet_all_eta25);



  scsumet_ee_all_eta15_pu10->Fill(eescsumet_all_eta15_pu10);
  scsumet_ee_all_eta20_pu10->Fill(eescsumet_all_eta20_pu10);
  scsumet_ee_all_eta25_pu10->Fill(eescsumet_all_eta25_pu10);

  scsumet_ee_all_eta15_pu20->Fill(eescsumet_all_eta15_pu20);
  scsumet_ee_all_eta20_pu20->Fill(eescsumet_all_eta20_pu20);
  scsumet_ee_all_eta25_pu20->Fill(eescsumet_all_eta25_pu20);

  scsumet_ee_all_eta15_pu30->Fill(eescsumet_all_eta15_pu30);
  scsumet_ee_all_eta20_pu30->Fill(eescsumet_all_eta20_pu30);
  scsumet_ee_all_eta25_pu30->Fill(eescsumet_all_eta25_pu30);

  scsumet_ee_all_eta15_pueq10->Fill(eescsumet_all_eta15_pueq10);
  scsumet_ee_all_eta20_pueq10->Fill(eescsumet_all_eta20_pueq10);
  scsumet_ee_all_eta25_pueq10->Fill(eescsumet_all_eta25_pueq10);

  scsumet_ee_all_eta15_pueq20->Fill(eescsumet_all_eta15_pueq20);
  scsumet_ee_all_eta20_pueq20->Fill(eescsumet_all_eta20_pueq20);
  scsumet_ee_all_eta25_pueq20->Fill(eescsumet_all_eta25_pueq20);



  scsumet_ee_all_eta15_pueq01->Fill(eescsumet_all_eta15_pueq01);
  scsumet_ee_all_eta20_pueq01->Fill(eescsumet_all_eta20_pueq01);
  scsumet_ee_all_eta25_pueq01->Fill(eescsumet_all_eta25_pueq01);

  scsumet_ee_all_eta15_pueq02->Fill(eescsumet_all_eta15_pueq02);
  scsumet_ee_all_eta20_pueq02->Fill(eescsumet_all_eta20_pueq02);
  scsumet_ee_all_eta25_pueq02->Fill(eescsumet_all_eta25_pueq02);

  scsumet_ee_all_eta15_pueq03->Fill(eescsumet_all_eta15_pueq03);
  scsumet_ee_all_eta20_pueq03->Fill(eescsumet_all_eta20_pueq03);
  scsumet_ee_all_eta25_pueq03->Fill(eescsumet_all_eta25_pueq03);

  scsumet_ee_all_eta15_pueq04->Fill(eescsumet_all_eta15_pueq04);
  scsumet_ee_all_eta20_pueq04->Fill(eescsumet_all_eta20_pueq04);
  scsumet_ee_all_eta25_pueq04->Fill(eescsumet_all_eta25_pueq04);

  scsumet_ee_all_eta15_pueq05->Fill(eescsumet_all_eta15_pueq05);
  scsumet_ee_all_eta20_pueq05->Fill(eescsumet_all_eta20_pueq05);
  scsumet_ee_all_eta25_pueq05->Fill(eescsumet_all_eta25_pueq05);

  scsumet_ee_all_eta15_pueq06->Fill(eescsumet_all_eta15_pueq06);
  scsumet_ee_all_eta20_pueq06->Fill(eescsumet_all_eta20_pueq06);
  scsumet_ee_all_eta25_pueq06->Fill(eescsumet_all_eta25_pueq06);

  scsumet_ee_all_eta15_pueq07->Fill(eescsumet_all_eta15_pueq07);
  scsumet_ee_all_eta20_pueq07->Fill(eescsumet_all_eta20_pueq07);
  scsumet_ee_all_eta25_pueq07->Fill(eescsumet_all_eta25_pueq07);

  scsumet_ee_all_eta15_pueq08->Fill(eescsumet_all_eta15_pueq08);
  scsumet_ee_all_eta20_pueq08->Fill(eescsumet_all_eta20_pueq08);
  scsumet_ee_all_eta25_pueq08->Fill(eescsumet_all_eta25_pueq08);

  scsumet_ee_all_eta15_pueq09->Fill(eescsumet_all_eta15_pueq09);
  scsumet_ee_all_eta20_pueq09->Fill(eescsumet_all_eta20_pueq09);
  scsumet_ee_all_eta25_pueq09->Fill(eescsumet_all_eta25_pueq09);





  scsumet_ee_all_gt2->Fill(eescsumet_all_gt2);
  scsumet_ee_severity0_gt2->Fill(eescsumet_severity0_gt2);
  scsumet_ee_koot0_gt2->Fill(eescsumet_koot0_gt2);

  scsumet_ee_all_gt5->Fill(eescsumet_all_gt5);
  scsumet_ee_severity0_gt5->Fill(eescsumet_severity0_gt5);
  scsumet_ee_koot0_gt5->Fill(eescsumet_koot0_gt5);

  scsumet_ee_all_gt10->Fill(eescsumet_all_gt10);
  scsumet_ee_severity0_gt10->Fill(eescsumet_severity0_gt10);
  scsumet_ee_koot0_gt10->Fill(eescsumet_koot0_gt10);





  

  // trigger primitive loop

   edm::Handle<EcalTrigPrimDigiCollection> TP;
   event.getByLabel("simEcalTriggerPrimitiveDigis",TP);

   edm::Handle<EcalTrigPrimDigiCollection> offlineTPsH;
   event.getByLabel("simEcalTriggerPrimitiveDigis",offlineTPsH);
   EcalTrigPrimDigiCollection const &offlineTPs = *(offlineTPsH.product());

   edm::Handle<EcalTrigPrimDigiCollection> onlineTPsH;
   event.getByLabel(onlineTPTag_, onlineTPsH);
   EcalTrigPrimDigiCollection const &onlineTPs  = *(onlineTPsH.product());

   EcalTrigPrimDigiCollection::const_iterator onIt_d, offIt_d;
   for(onIt_d = onlineTPs.begin(); onIt_d != onlineTPs.end(); ++onIt_d){
     if(onIt_d->compressedEt() > 0) {
       //std::cout << "looping on online TPG" << std::endl;
       //std::cout << "ET=" << onIt_d->compressedEt() << std::endl;
       h_ttonline_et->Fill(onIt_d->compressedEt());
       //EcalTPGTowerStatusMapIterator towIt = towerMap.find(onIt_d->id().rawId());
       //std::cout << "TP reading " << towIt->second << std::endl;
     }//non zero TP
   }//loop online TPG

   for(offIt_d = offlineTPs.begin(); offIt_d != offlineTPs.end(); ++offIt_d){ 
     if(offIt_d->compressedEt() > 0) {
       h_ttoffline_et->Fill(offIt_d->compressedEt());
     }//non zero TP
   }//loop offline TPG




   numtp_eb=0, numtp_ee=0;

   numtp_samp2_eb=0, numtp_samp2_ee=0;

   numtp_geq1_eb=0, numtp_geq1_ee=0;

   tpsumet_eb=0, tpsumet_ee=0;

   tpsumet_samp2_eb=0, tpsumet_samp2_ee=0;

   tpsumet_geq1_eb=0, tpsumet_geq1_ee=0;

   tpsumet_mod1_eb=0,  tpsumet_mod2_eb=0,  tpsumet_mod3_eb=0,  tpsumet_mod4_eb=0;

   tpsumet_geq1_mod1_eb=0,  tpsumet_geq1_mod2_eb=0,  tpsumet_geq1_mod3_eb=0,  tpsumet_geq1_mod4_eb=0;

   tpsumet_ring1_ee=0, tpsumet_ring2_ee=0, tpsumet_ring3_ee=0;

   tpsumet_geq1_ring1_ee=0, tpsumet_geq1_ring2_ee=0, tpsumet_geq1_ring3_ee=0;




   for (EcalTrigPrimDigiCollection::const_iterator tp=TP->begin(); tp!=TP->end(); ++tp) {



       int tp_ieta = tp->id().ieta();
       int tp_iphi = tp->id().iphi();

       int tp_sfgvb=0;

       int etmax=0;
       int sampmax=0;

       int tp_compressed_et=0;





       for (Int_t j=0;j<5;j++) {
	  
	   EcalTriggerPrimitiveSample sam=tp->sample(j);
	   
	   int et_tmp=sam.compressedEt();

	   if (et_tmp > etmax) {
	     etmax=et_tmp;
	     sampmax=j;
	     tp_sfgvb=sam.sFGVB();
	   }

       }

       tp_compressed_et=etmax;







      

   
       //       correction for towers at ieta=27,28


         if (abs(tp_ieta)==27 || abs(tp_ieta)==28) tp_compressed_et*=2;


       float tp_et=tp_compressed_et*0.5;


       if (tp_compressed_et>0) tp_map->Fill(tp_iphi,tp_ieta,tp_et);
       if (tp_compressed_et>0 && sampmax==2) tp_map_samp2->Fill(tp_iphi,tp_ieta,tp_et);
       if (tp_compressed_et>0 && tp_sfgvb==1) tp_map2->Fill(tp_iphi,tp_ieta,tp_et);


       if (tp_compressed_et>0) {
	 numtp_vs_ieta->Fill(float(tp_ieta)-0.5+1*(tp_ieta<0));
	 if (tp_compressed_et>=2) numtp_geq1_vs_ieta->Fill(float(tp_ieta)-0.5+1*(tp_ieta<0));
	 numtp_etweighted_vs_ieta->Fill(float(tp_ieta)-0.5+1*(tp_ieta<0),tp_et);
       }

       if (tp_compressed_et>0 && sampmax==2) {
	 numtp_vs_ieta_samp2->Fill(float(tp_ieta)-0.5+1*(tp_ieta<0));
       }



       if (tp->id().subDet()==EcalBarrel) {
   
 
	 if (tp_compressed_et>0) {


	 
	   numtp_eb++;
	   tpsumet_eb+=tp_et;

	   if (sampmax==2) {
	     numtp_samp2_eb++;
	     tpsumet_samp2_eb+=tp_et;
	   }


	   if (tp_compressed_et>=2) {
	     numtp_geq1_eb++;
	     tpsumet_geq1_eb+=tp_et;
	   }

	   if (abs(tp_ieta)<=5) {
	     tpsumet_mod1_eb+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_mod1_eb+=tp_et;
	   }

	   if (abs(tp_ieta)>=6 && abs(tp_ieta)<=9) {
	     tpsumet_mod2_eb+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_mod2_eb+=tp_et;
	   }

	   if (abs(tp_ieta)>=10 && abs(tp_ieta)<=13) {
	     tpsumet_mod3_eb+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_mod3_eb+=tp_et;
	   }

	   if (abs(tp_ieta)>=14) {
	     tpsumet_mod4_eb+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_mod4_eb+=tp_et;
	   }
	 
	   //   cout << " Barrel TP:  ieta=" << tp_ieta << " iphi=" << tp_iphi << " compressedEt=" << tp_compressed_et << " sfgvb=" << tp_sfgvb << endl;
//  	   cout << "Barrel:  mod1=" << tpsumet_mod1_eb 
//  		<< " mod2=" << tpsumet_mod2_eb
//  		<< " mod3=" << tpsumet_mod3_eb
//  		<< " mod4=" << tpsumet_mod4_eb
// 		<< " mod 1-4" << tpsumet_mod1_eb+tpsumet_mod2_eb+tpsumet_mod3_eb+tpsumet_mod4_eb
// 		<< " tot=" << tpsumet_eb << endl;
	 }

       }

       if (tp->id().subDet()==EcalEndcap) {
    
	 if (tp_compressed_et>0) {
	

 

	   numtp_ee++;
	   tpsumet_ee+=tp_et;

	   if (sampmax==2) {
	     numtp_samp2_ee++;
	     tpsumet_samp2_ee+=tp_et;
	     sumtrigprim_ee_bx0+=tp_et;
	   }


	   if (sampmax==0) {
	     sumtrigprim_ee_bxm2+=tp_et;
	   }

	   if (sampmax==1) {
	     sumtrigprim_ee_bxm1+=tp_et;
	   }

	   if (sampmax==3) {
	     sumtrigprim_ee_bxp1+=tp_et;
	   }

	   if (sampmax==4) {
	     sumtrigprim_ee_bxp2+=tp_et;
	   }




	   if (tp_compressed_et>=2) {
	     numtp_geq1_ee++;
	     tpsumet_geq1_ee+=tp_et;
	   }


	   if (abs(tp_ieta)<=21) {
	     tpsumet_ring1_ee+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_ring1_ee+=tp_et;
	   }

	   if (abs(tp_ieta)>=22 && abs(tp_ieta)<=25) {
	     tpsumet_ring2_ee+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_ring2_ee+=tp_et;
	   }

	   if (abs(tp_ieta)>=26) {
	     tpsumet_ring3_ee+=tp_et;
	     if (tp_compressed_et>=2) tpsumet_geq1_ring3_ee+=tp_et;
	   }


	  

	 }

       }




   }

  

   // RecHit loop: EB

  EcalRecHit maxebhit;
  EBDetId  maxebid(0);
  EcalRecHit maxeehit;
  EcalRecHit maxebhit2;
  EcalRecHit maxeehit2;


  float ensum1=0;
  float ensum2=0;
  float ensum3=0;
  float ensum4=0;

  float errsum1=0;
  float errsum2=0;

  int numgt1=0;



  rechitsumet_eb_all=0;
  rechitsumet_eb_01=0;
  rechitsumet_eb_05=0;

  rechitsumet_eb_0105=0;
  rechitsumet_eb_0530=0;


  rechitsumet_ee_all=0;
  rechitsumet_ee_01=0;
  rechitsumet_ee_05=0;

  rechitsumet_ee_0105=0;
  rechitsumet_ee_0530=0;



  for (EcalRecHitCollection::const_iterator hitItr = EBhits->begin(); hitItr != EBhits->end(); ++hitItr) {

      EcalRecHit hit = (*hitItr);
      EBDetId det    = hit.id(); 
      const DetId det2     = hit.id(); 
      float ampli    = hit.energy();
      float time     = hit.time()-toffset;
      int   ebflag   = hit.recoFlag();
      int sm         = det.ism();

      int ebhashedid = det.hashedIndex();

      const EcalPedestals::Item *aped=0;
      aped=&eped->barrel(ebhashedid);
      float pedg12=aped->mean_x12;
      //     cout << "pedestal g12=" << pedg12 << endl;

      float_t chi2   = hit.chi2();
      float_t chi2oot= hit.outOfTimeChi2();
      //  float lasercalib_eb=laser->getLaserCorrection(EBDetId(det), event.time());

      int chanstatus=-1;
 
      EcalChannelStatus::const_iterator chIt = cstat->find( det );
      if ( chIt != cstat->end() ) {
	chanstatus = chIt->getStatusCode() & 0x1F;
      } 
	

      GlobalPoint poseb=geo->getPosition(hit.detid());


      float eta_eb=poseb.eta();
      float phi_eb=poseb.phi();
      float pf=1.0/cosh(eta_eb);
      float eteb=ampli*pf;
      int ieta=det.ieta();      
      int iphi=det.iphi();


          

      float eb_r4t=0;
      float eb_r9t=0;

      float eb_r4a=0;


      // r4 for all hits
      float e4x1a=0;

   


      if (ampli!=0) eb_r4a=e4x1a/ampli;
      eb_r4a=1-eb_r4a;

        


      // calculate e2/e9, e2/e25 for all hits with E>5.0

	Float_t emax=-999;
	Int_t n1max=-111;
	Int_t n2max=-111;

	if (pTopology.isValid() && pG.isValid() && eteb>3.0) { 


	  if (eb_r4a>0.95) r4count++;
	  
	  const CaloTopology *topology=pTopology.product();
	  
	  const reco::BasicCluster *cluster=0;
	  /*
	    float e3x3tmp=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),-1,1,-1,1);
	    
	    float e5x5tmp=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),-2,2,-2,2);
	  */
	  
	  
	  
	  
	}





	Int_t ee_kGood=0;
	Int_t ee_kPoorReco=0;
	Int_t ee_kOutOfTime=0;
	Int_t ee_kFaultyHardware=0;
	Int_t ee_kNoisy=0;
	Int_t ee_kPoorCalib=0;
	Int_t ee_kSaturated=0;
	Int_t ee_kLeadingEdgeRecovered=0;
	Int_t ee_kNeighboursRecovered=0;
	Int_t ee_kTowerRecovered=0;
	Int_t ee_kDead=0;
	Int_t ee_kKilled=0;
	Int_t ee_kTPSaturated=0;
	Int_t ee_kL1SpikeFlag=0;
	Int_t ee_kWeird=0;
	Int_t ee_kDiWeird=0;
	Int_t ee_kUnknown=0;
	
	
	
	if (hit.checkFlag(EcalRecHit::kGood)) ee_kGood=1;
	if (hit.checkFlag(EcalRecHit::kPoorReco)) ee_kPoorReco=1;
	if (hit.checkFlag(EcalRecHit::kOutOfTime)) ee_kOutOfTime=1;
	if (hit.checkFlag(EcalRecHit::kFaultyHardware)) ee_kFaultyHardware=1;
	if (hit.checkFlag(EcalRecHit::kNoisy)) ee_kNoisy=1;
	if (hit.checkFlag(EcalRecHit::kPoorCalib)) ee_kPoorCalib=1;
	if (hit.checkFlag(EcalRecHit::kSaturated)) ee_kSaturated=1;
	if (hit.checkFlag(EcalRecHit::kLeadingEdgeRecovered)) ee_kLeadingEdgeRecovered=1;
	if (hit.checkFlag(EcalRecHit::kNeighboursRecovered)) ee_kNeighboursRecovered=1;
	if (hit.checkFlag(EcalRecHit::kTowerRecovered)) ee_kTowerRecovered=1;
	if (hit.checkFlag(EcalRecHit::kDead)) ee_kDead=1;
	if (hit.checkFlag(EcalRecHit::kKilled)) ee_kKilled=1;
	if (hit.checkFlag(EcalRecHit::kTPSaturated)) ee_kTPSaturated=1;
	if (hit.checkFlag(EcalRecHit::kL1SpikeFlag)) ee_kL1SpikeFlag=1;
	if (hit.checkFlag(EcalRecHit::kWeird)) ee_kWeird=1;
	if (hit.checkFlag(EcalRecHit::kDiWeird)) ee_kDiWeird=1;
	if (hit.checkFlag(EcalRecHit::kUnknown)) ee_kUnknown=1;
	



	// check for hits close to boundaries
	
	//if (abs(ieta)==85) ee_kGood=0;
	//if (EcalTools::isNextToDeadFromNeighbours(det,*cstat,12)) ee_kGood=0;
	

	// mask out noisy channels in run D

	if ((ieta==11 && iphi==68) || (ieta==68 && iphi==74) || (ieta==-83 && iphi==189) || (ieta==-75 && iphi==199)) ee_kWeird=1;

	
	
	if (!(ee_kWeird || ee_kDiWeird)) {




	  EcalTrigTowerDetId ebttid=det.tower();

	  int towerieta=ebttid.ieta();
	  int toweriphi=ebttid.iphi();
	  int tcc=ebttid.iDCC();
	  int tt=ebttid.iTT();
	  rechit_map->Fill(toweriphi,towerieta,eteb);



 

	  if (eteb>0.1) ebtime_vs_bxtrain_01->Fill(bunchintrain+0.5,time);
	  if (eteb>0.5) ebtime_vs_bxtrain_05->Fill(bunchintrain+0.5,time);

	  
	  if (ampli>1.0) { 
	    ebsum_gt1+=ampli;
	    ebhits1GeV++;
	  }
	  
	  if (ampli>2.0) {
	    ebsum_gt2+=ampli;
	    ebhits2GeV++;
	  }
	  
	  if (ampli>4.0) {
	    ebsum_gt4+=ampli;
	    ebhits4GeV++;
	  }
	  
	  
	  if (eteb>1.0) {
	    ebsum_gt1et+=eteb;
	    ebhits1GeVet++;
	  }
	  
	  
	  if (eteb>2.0) {
	    ebsum_gt2et+=eteb;
	    ebhits2GeVet++;
	  }
	  
	  if (eteb>4.0) {
	    ebsum_gt4et+=eteb;
	    ebhits4GeVet++;
	  }
	  
	
	// rechit et sums

	rechitsumet_eb_all+=eteb;
	if (eteb>0.1) rechitsumet_eb_01+=eteb;
	if (eteb>0.5) rechitsumet_eb_05+=eteb;


	if (eteb>0.1 && eteb<=0.5) rechitsumet_eb_0105+=eteb;
	if (eteb>0.5 && eteb<=3.0) rechitsumet_eb_0530+=eteb;



	if (eteb>0.1) ebnumrechits_01++;
	if (eteb>0.5) ebnumrechits_05++;


	if (eteb>0.1 && eteb<=0.5) ebnumrechits_0105++;
	if (eteb>0.5 && eteb<=3.0) ebnumrechits_0530++;

 	if (eteb>0.1) rechiteta_vs_bxtrain_01->Fill(bunchintrain,eta_eb);
 	if (eteb>0.5) rechiteta_vs_bxtrain_05->Fill(bunchintrain,eta_eb);


	// *** DIGI LOOP - EB ***

	if (eteb>0.1) {



	  float pedoff=0;

	  if (isMC_) pedoff=0;
	  if (!isMC_) pedoff=200-pedg12;

	  
	  //cout << "EB DIGI printout:\n\n";
	 
	  //  cout << "rechit energy=" << ampli << " rechit et=" << eteb << endl; 
	  //cout << "SAMPLE   ADC   GAIN\n";
	  
	  Int_t ebdigiadc=0;
	  Int_t ebdigigain=0;
	  char mytxteb[80];
	  
	  // find digi corresponding to rechit

	  EBDigiCollection::const_iterator digiItrEB= EBdigis->begin();
	  while(digiItrEB != EBdigis->end() && digiItrEB->id() != hitItr->id())
	    {
	      digiItrEB++;
	    }
	  
	  if (digiItrEB != EBdigis->end()) {
	    
	    EBDataFrame df(*digiItrEB);
	    for(int i=0; i<10;++i)
	      {
		ebdigiadc=df.sample(i).adc();
		ebdigiadc+=pedoff;
		ebdigigain=df.sample(i).gainId();
		//		sprintf(mytxteb,"  %02d    %04d     %d",i+1,ebdigiadc,ebdigigain);
		//		cout << mytxteb << endl;
		

		eb_digi_01->Fill(i+0.5,ebdigiadc);
		
		if (eteb>0.5) eb_digi_05->Fill(i+0.5,ebdigiadc);


		if (eteb>0.1 && eteb<=0.5) { 
		  eb_digi_0105->Fill(i+0.5,float(ebdigiadc));
		  eb_digi_0105_vs_time->Fill(time,i+0.5,float(ebdigiadc));
		  eb_digi_0105_vs_bxtrain->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));

		  eb_digi_0105_vs_time_norm->Fill(time,i+0.5,1.0);
		  eb_digi_0105_vs_bxtrain_norm->Fill(bunchintrain+0.5,i+0.5,1.0);


		  if (abs(eb_eta)<0.5) {

		    eb_digi_0105_vs_time_eta15->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0105_vs_bxtrain_eta15->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0105_vs_time_norm_eta15->Fill(time,i+0.5,1.0);
		    eb_digi_0105_vs_bxtrain_norm_eta15->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(eb_eta)>=0.5 && abs(eb_eta)<1.0) {

		    eb_digi_0105_vs_time_eta20->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0105_vs_bxtrain_eta20->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0105_vs_time_norm_eta20->Fill(time,i+0.5,1.0);
		    eb_digi_0105_vs_bxtrain_norm_eta20->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(eb_eta)>=1.0) {

		    eb_digi_0105_vs_time_eta25->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0105_vs_bxtrain_eta25->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0105_vs_time_norm_eta25->Fill(time,i+0.5,1.0);
		    eb_digi_0105_vs_bxtrain_norm_eta25->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }



		}


		if (eteb>0.5 && eteb<=3.0) {

		  eb_digi_0530->Fill(i+0.5,ebdigiadc);
		  eb_digi_0530_vs_time->Fill(time,i+0.5,float(ebdigiadc));
		  eb_digi_0530_vs_bxtrain->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));

		  eb_digi_0530_vs_time_norm->Fill(time,i+0.5,1.0);
		  eb_digi_0530_vs_bxtrain_norm->Fill(bunchintrain+0.5,i+0.5,1.0);


		  if (abs(eb_eta)<0.5) {

		    eb_digi_0530_vs_time_eta15->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0530_vs_bxtrain_eta15->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0530_vs_time_norm_eta15->Fill(time,i+0.5,1.0);
		    eb_digi_0530_vs_bxtrain_norm_eta15->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(eb_eta)>=0.5 && abs(eb_eta)<1.0) {

		    eb_digi_0530_vs_time_eta20->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0530_vs_bxtrain_eta20->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0530_vs_time_norm_eta20->Fill(time,i+0.5,1.0);
		    eb_digi_0530_vs_bxtrain_norm_eta20->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(eb_eta)>=1.0) {

		    eb_digi_0530_vs_time_eta25->Fill(time,i+0.5,float(ebdigiadc));
		    eb_digi_0530_vs_bxtrain_eta25->Fill(bunchintrain+0.5,i+0.5,float(ebdigiadc));
		    
		    eb_digi_0530_vs_time_norm_eta25->Fill(time,i+0.5,1.0);
		    eb_digi_0530_vs_bxtrain_norm_eta25->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }




	
		}


	      }
	    
	  }

	}
	


	}
  
	  
	  
	  
	  
   
        eb_timing_0->Fill(time);
        eb_r4_0->Fill(eb_r4a);
        eb_timing_vs_r4_0->Fill(eb_r4a,time);
        if (eb_r4a<0.95) eb_timing_r4_0->Fill(time);

    

        if (ampli>0.2) {
           eb_timing_200->Fill(time);
           eb_r4_200->Fill(eb_r4a);
           eb_timing_vs_r4_200->Fill(eb_r4a,time);
           if (eb_r4a<0.95) eb_timing_r4_200->Fill(time);

        }


        if (ampli>0.4) {
           eb_timing_400->Fill(time);
           eb_r4_400->Fill(eb_r4a);
           eb_timing_vs_r4_400->Fill(eb_r4a,time);
           if (eb_r4a<0.95) eb_timing_r4_400->Fill(time);

        }


        if (ampli>0.6) {
           eb_timing_600->Fill(time);
           eb_r4_600->Fill(eb_r4a);
           eb_timing_vs_r4_600->Fill(eb_r4a,time);
           if (eb_r4a<0.95) eb_timing_r4_600->Fill(time);
        }


        if (ampli>0.8) {
           eb_timing_800->Fill(time);
           eb_r4_800->Fill(eb_r4a);
           eb_timing_vs_r4_800->Fill(eb_r4a,time);
           if (eb_r4a<0.95) eb_timing_r4_800->Fill(time);
        }


        if (eteb>3.00) {

	  float sigmat=28.51/(ampli/0.04)+sqrt(2)*0.2565;

           eb_timing_1000->Fill(time);
           eb_r4_1000->Fill(eb_r4a);
           eb_timing_vs_r4_1000->Fill(eb_r4a,time);
           if (eb_r4a<0.95) eb_timing_r4_1000->Fill(time);

           numgt1++;
           ensum1+=time*ampli;
           ensum2+=ampli;
	   ensum3+=pow(ampli*sigmat,2);
	   ensum4+=pow(ampli,2);

           errsum1+=time/pow(sigmat,2);
           errsum2+=1.0/pow(sigmat,2);



	     int koutoftime=0;
	      if (hit.recoFlag()==EcalRecHit::kOutOfTime) koutoftime=1;
	      //	      cout << "koutoftime=" << koutoftime << endl;

	   
	   int sevlev=severity->severityLevel(det2,*ebRecHits);

    
	   int koot=0;
	   int kweird=0;
	   
	   if (hit.checkFlag(EcalRecHit::kOutOfTime)) koot=1;
	   if (hit.checkFlag(EcalRecHit::kWeird)) kweird=1;
	   
	   // cout << "recHit kOutOfTime=" << koot << endl;
	   //  cout << "recHit kWeird=" << kweird << endl;
	   
	   
	   int kweird_cleaning=0;
	   
	   if (cleaningAlgo_->checkTopology(det2,*ebRecHits)==EcalRecHit::kWeird) kweird_cleaning=1;
	   
	   
	   //  cout << "cleaning algo kWeird=" << kweird_cleaning << endl;
	   
	   
	   
	   // cout << "swisscross1=" << severity.spikeFromNeighbours(det2,*ebRecHits,5.0,EcalSeverityLevelAlgo::kSwissCross) << endl;
	   // cout << "swisscross2=" << severity.spikeFromNeighbours(det2,*ebRecHits,5.0,EcalSeverityLevelAlgo::kSwissCrossBordersIncluded) << endl;
	   
	   //             int sevlev=severity.severityLevel(det2,*ebRecHits,*cstat,5.0,EcalSeverityLevelAlgo::kSwissCrossBordersIncluded);
	   
	   // cout << "SeverityLevel=" << sevlev << endl;
	   
        }
	

	if (ampli>2.0) {
	  eb_timing_2000->Fill(time);
	  eb_r4_2000->Fill(eb_r4a);
	  eb_timing_vs_r4_2000->Fill(eb_r4a,time);
	  if (eb_r4a<0.95) eb_timing_r4_2000->Fill(time);
        }
	
	
	if (ampli>3.0) {
	  eb_timing_3000->Fill(time);
	  eb_r4_3000->Fill(eb_r4a);
	  eb_timing_vs_r4_3000->Fill(eb_r4a,time);
	  if (eb_r4a<0.95) eb_timing_r4_3000->Fill(time);
        }
	
	
	if (ampli>5.0) {
	  eb_timing_5000->Fill(time);
	  eb_r4_5000->Fill(eb_r4a);
	  eb_timing_vs_r4_5000->Fill(eb_r4a,time);
	  if (eb_r4a<0.95) eb_timing_r4_5000->Fill(time);
        }
	
	
	
	
	if (ampli>0.0) {
	  
	  if (pTopology.isValid() && pG.isValid()) { 
	    
	    const CaloTopology *topology=pTopology.product();
	    
	    const reco::BasicCluster *cluster=0;
	    
	    
	    float e4x1t=0;
	    /*
	      e4x1t+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),-1,-1,0,0);
	      e4x1t+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),1,1,0,0);
	      e4x1t+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),0,0,-1,-1);
	      e4x1t+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,hit.id(),0,0,+1,+1);
	    */
	    
	    if (ampli!=0) eb_r4t=e4x1t/ampli;
	    
	    eb_r4t=1-eb_r4t;
	    
	    //         float e3x3t=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),-1,1,-1,1);
	    
	    
	    //         if (e3x3t>0) eb_r9t=ampli/e3x3t;
	    
	    
	  }
	}
	
	//     if (ampli>ebmax && ((ampli<3) || (ampli>3 && abs(time)<4 && eb_r9t<0.90))) {

	if (ampli>ebmax && ee_kGood) {
	  ebmax=ampli;
	  ebmaxet=eteb;
	  eb_ieta=ieta;      
	  eb_iphi=iphi;      
	  eb_eta=eta_eb;
	  eb_phi=phi_eb;
	  ebtime=time;
	  ebflags=ebflag;
	  ebchi2=chi2;
	  ebchi2oot=chi2oot;
	  maxebhit=*hitItr;
	  maxebid=det;
	  
	  //       swisscross_algo=swisscross;
	  //  e2e9_algo=e2overe9tmp;
	  
	}
	
	
	
	
	Float_t ieta_tmp=float(ieta)+0.5;
	if (ieta>0) ieta_tmp=ieta_tmp-1;
	Float_t iphi_tmp=float(iphi)-0.5;
	


	
	if (!(ee_kWeird || ee_kDiWeird)) {
	
	

        eb_rechitenergy_->Fill(ampli);



	rechiteta_all->Fill(eta_eb);
	rechiteta_etweight->Fill(eta_eb,eteb);

	if (eteb>1.0) {
	  rechiteta_gt1et->Fill(eta_eb);
	  rechiteta_etweight_gt1et->Fill(eta_eb,eteb);

	  if (numgoodvtx<10) 	  rechiteta_gt1et_pu10->Fill(eta_eb);
	  if (numgoodvtx>=10 && numgoodvtx<20) 	  rechiteta_gt1et_pu20->Fill(eta_eb);
	  if (numgoodvtx>20) 	  rechiteta_gt1et_pu30->Fill(eta_eb);

	}


        if (fabs(eta_eb)<0.2) eb_rechitenergy_02->Fill(ampli);
	if (fabs(eta_eb)>=0.2 && fabs(eta_eb)<0.4) eb_rechitenergy_04->Fill(ampli);
	if (fabs(eta_eb)>=0.4 && fabs(eta_eb)<0.6) eb_rechitenergy_06->Fill(ampli);
	if (fabs(eta_eb)>=0.6 && fabs(eta_eb)<0.8) eb_rechitenergy_08->Fill(ampli);
	if (fabs(eta_eb)>=0.7 && fabs(eta_eb)<1.0) eb_rechitenergy_10->Fill(ampli);
	if (fabs(eta_eb)>=1.0 && fabs(eta_eb)<1.2) eb_rechitenergy_12->Fill(ampli);
	if (fabs(eta_eb)>=1.2 && fabs(eta_eb)<1.4) eb_rechitenergy_14->Fill(ampli);
	if (fabs(eta_eb)>=1.4) eb_rechitenergy_148->Fill(ampli);

      

      
        eb_rechitet_->Fill(eteb);


        if (fabs(eta_eb)<0.2) eb_rechitet_02->Fill(eteb);
	if (fabs(eta_eb)>=0.2 && fabs(eta_eb)<0.4) eb_rechitet_04->Fill(eteb);
	if (fabs(eta_eb)>=0.4 && fabs(eta_eb)<0.6) eb_rechitet_06->Fill(eteb);
	if (fabs(eta_eb)>=0.6 && fabs(eta_eb)<0.8) eb_rechitet_08->Fill(eteb);
	if (fabs(eta_eb)>=0.7 && fabs(eta_eb)<1.0) eb_rechitet_10->Fill(eteb);
	if (fabs(eta_eb)>=1.0 && fabs(eta_eb)<1.2) eb_rechitet_12->Fill(eteb);
	if (fabs(eta_eb)>=1.2 && fabs(eta_eb)<1.4) eb_rechitet_14->Fill(eteb);
	if (fabs(eta_eb)>=1.4) eb_rechitet_148->Fill(eteb);

 


        if (fabs(eta_eb)<0.5) eb_rechitetvspu_05->Fill(float(numgoodvtx)-0.5,eteb);
	if (fabs(eta_eb)>=0.5 && fabs(eta_eb)<1.0) eb_rechitetvspu_10->Fill(float(numgoodvtx)-0.5,eteb);
	if (fabs(eta_eb)>=1.0 && fabs(eta_eb)<1.5) eb_rechitetvspu_15->Fill(float(numgoodvtx)-0.5,eteb);
       



        ebocc->Fill(iphi_tmp,ieta_tmp,1.);
        eboccen->Fill(iphi_tmp,ieta_tmp,ampli);
        eboccet->Fill(iphi_tmp,ieta_tmp,eteb);


	
        eb_rechiten_vs_eta->Fill(eta_eb,ampli);
        eb_rechitet_vs_eta->Fill(eta_eb,eteb);
	

        if (eteb>1.0) {

           eboccgt1et->Fill(iphi_tmp,ieta_tmp,1.);
           eboccetgt1et->Fill(iphi_tmp,ieta_tmp,eteb);

        }

	if (ampli>1.0) {

           eboccgt1->Fill(iphi_tmp,ieta_tmp,1.);
           eboccengt1->Fill(iphi_tmp,ieta_tmp,ampli);

        }




	

      ebhits++;
      if (ampli>1.0) ebhits1GeV++;
	}
	
  }


  // end of eb loop



  if (numgt1>0) {

    tmean_en=ensum1/ensum2;
    terr_en=sqrt(ensum3/ensum4);

    tmean_sig=errsum1/errsum2;
    terr_sig=sqrt(1/errsum2);

  }




  // EE rechits


	   float_t badsc1_et=0;
	   Int_t badsc1_hits=0;

	   float_t badsc2_et=0;
	   Int_t badsc2_hits=0;


  for (EcalRecHitCollection::const_iterator hitItr = EEhits->begin(); hitItr != EEhits->end(); ++hitItr) {

      EcalRecHit hit = (*hitItr);
      EEDetId det = hit.id(); 
            
      int dee=0;

      //   float lasercalib_ee=laser->getLaserCorrection(EEDetId(det), event.time());

      float ampli = hit.energy();
      float time     = hit.time()-toffset;
      int   ebflag   = hit.recoFlag();
      float_t chi2   = hit.chi2();
      float_t chi2oot= hit.outOfTimeChi2();


      int eehashedid = det.hashedIndex();

      const EcalPedestals::Item *aped=0;
      aped=&eped->endcap(eehashedid);
      float pedg12=aped->mean_x12;
      //cout << "EE pedestal g12=" << pedg12 << endl;



      GlobalPoint posee=geo->getPosition(hit.detid());
      float eta_ee=posee.eta();
      float phi_ee=posee.phi();
      float pf=1.0/cosh(eta_ee);
      float etee=ampli*pf;
      int ix=det.ix();
      int iy=det.iy();
      int side=det.zside();

      int iz=0;
      if (side==1) iz=1;
      if (side==-1) iz=-1;


      if (ix<51 && side==1) dee=1;
      if (ix>50 && side==1) dee=2;
      if (ix<51 && side==-1) dee=4;
      if (ix>50 && side==-1) dee=3;

      float ee_r4a=0;
      float ee4x1a=0;

      if (pTopology.isValid() && pG.isValid()) { 

        const CaloTopology *topology=pTopology.product();
  
        const reco::BasicCluster *cluster=0;
	/*
        ee4x1a+=EcalClusterTools::matrixEnergy(*cluster,eeRecHits,topology,hit.id(),-1,-1,0,0);
        ee4x1a+=EcalClusterTools::matrixEnergy(*cluster,eeRecHits,topology,hit.id(),1,1,0,0);
        ee4x1a+=EcalClusterTools::matrixEnergy(*cluster,eeRecHits,topology,hit.id(),0,0,-1,-1);
        ee4x1a+=EcalClusterTools::matrixEnergy(*cluster,eeRecHits,topology,hit.id(),0,0,+1,+1);
	*/
      }


      if (ampli!=0) ee_r4a=ee4x1a/ampli;
      ee_r4a=1-ee_r4a;


        if (etee>3.0) {

	  float sigmat=36.08/(ampli/0.06)+sqrt(2)*0.18;




        }



      Int_t ee_kGood=0;
           Int_t ee_kPoorReco=0;
           Int_t ee_kOutOfTime=0;
           Int_t ee_kFaultyHardware=0;
           Int_t ee_kNoisy=0;
           Int_t ee_kPoorCalib=0;
           Int_t ee_kSaturated=0;
           Int_t ee_kLeadingEdgeRecovered=0;
           Int_t ee_kNeighboursRecovered=0;
           Int_t ee_kTowerRecovered=0;
           Int_t ee_kDead=0;
           Int_t ee_kKilled=0;
           Int_t ee_kTPSaturated=0;
           Int_t ee_kL1SpikeFlag=0;
           Int_t ee_kWeird=0;
           Int_t ee_kDiWeird=0;
           Int_t ee_kUnknown=0;



	   if (hit.checkFlag(EcalRecHit::kGood)) ee_kGood=1;
	   if (hit.checkFlag(EcalRecHit::kPoorReco)) ee_kPoorReco=1;
	   if (hit.checkFlag(EcalRecHit::kOutOfTime)) ee_kOutOfTime=1;
	   if (hit.checkFlag(EcalRecHit::kFaultyHardware)) ee_kFaultyHardware=1;
	   if (hit.checkFlag(EcalRecHit::kNoisy)) ee_kNoisy=1;
	   if (hit.checkFlag(EcalRecHit::kPoorCalib)) ee_kPoorCalib=1;
	   if (hit.checkFlag(EcalRecHit::kSaturated)) ee_kSaturated=1;
	   if (hit.checkFlag(EcalRecHit::kLeadingEdgeRecovered)) ee_kLeadingEdgeRecovered=1;
	   if (hit.checkFlag(EcalRecHit::kNeighboursRecovered)) ee_kNeighboursRecovered=1;
	   if (hit.checkFlag(EcalRecHit::kTowerRecovered)) ee_kTowerRecovered=1;
	   if (hit.checkFlag(EcalRecHit::kDead)) ee_kDead=1;
	   if (hit.checkFlag(EcalRecHit::kKilled)) ee_kKilled=1;
	   if (hit.checkFlag(EcalRecHit::kTPSaturated)) ee_kTPSaturated=1;
	   if (hit.checkFlag(EcalRecHit::kL1SpikeFlag)) ee_kL1SpikeFlag=1;
	   if (hit.checkFlag(EcalRecHit::kWeird)) ee_kWeird=1;
	   if (hit.checkFlag(EcalRecHit::kDiWeird)) ee_kDiWeird=1;
	   if (hit.checkFlag(EcalRecHit::kUnknown)) ee_kUnknown=1;



	   //   if (EcalTools::isNextToDeadFromNeighbours(det,*cstat,12)) ee_kGood=0;
	   // if (EEDetId::isNextToRingBoundary(det)) ee_kGood=0;




	// mask out noisy channels in run D

	   if (ix==58 && iy==94 && iz==-1)  ee_kWeird=1;


	   if (!(ee_kWeird)) {





	     int towerieta=ietamapTT->GetBinContent(ix,iy);
	     int toweriphi=iphimapTT->GetBinContent(ix,iy);

	     if (iz<0) towerieta*=-1;

	     rechit_map->Fill(toweriphi,towerieta,etee);
	     rechit_map_time->Fill(toweriphi,towerieta,time);

	     rechit_map_ixiy->SetBinContent(ix+100*(iz>0),iy,etee);
	     rechit_map_time_ixiy->SetBinContent(ix+100*(iz>0),iy,time);

	     //  cout << " EE Rechit: energy=" << ampli << " ix,iy,iz=" << ix << "," << iy << "," << iz << " tower ieta,iphi=" 
	     //	  << towerieta << "," << toweriphi <<endl;  




	     if (etee>0.1) eetime_vs_bxtrain_01->Fill(bunchintrain+0.5,time);
	     if (etee>0.5) eetime_vs_bxtrain_05->Fill(bunchintrain+0.5,time);


	     if (ampli>1.0) { 
	       eesum_gt1+=ampli;
	       eehits1GeV++;
	     }

	     if (ampli>2.0) {
	       eesum_gt2+=ampli;
	       eehits2GeV++;
	     }

	     if (ampli>4.0) {
	       eesum_gt4+=ampli;
	       eehits4GeV++;
	     }


	     if (etee>1.0) {
	       eesum_gt1et+=etee;
	       eehits1GeVet++;
	     }


	     if (etee>2.0) {
	       eesum_gt2et+=etee;
	       eehits2GeVet++;
	     }

	     if (etee>4.0) {
	       eesum_gt4et+=etee;
	       eehits4GeVet++;
	     }
	    


	     // rechit et sums

	     rechitsumet_ee_all+=etee;
	     if (etee>0.1) rechitsumet_ee_01+=etee;
	     if (etee>0.5) rechitsumet_ee_05+=etee;

	     if (etee>0.1 && etee<=0.5) rechitsumet_ee_0105+=etee;
	     if (etee>0.5 && etee<=3.0) rechitsumet_ee_0530+=etee;


	     if (etee>0.1) eenumrechits_01++;
	     if (etee>0.5) eenumrechits_05++;
	     

	     if (etee>0.1 && etee<=0.5) eenumrechits_0105++;
	     if (etee>0.5 && etee<=3.0) eenumrechits_0530++;

	     if (etee>0.1) rechiteta_vs_bxtrain_01->Fill(bunchintrain,eta_ee);
	     if (etee>0.5) rechiteta_vs_bxtrain_05->Fill(bunchintrain,eta_ee);
	   


	// ** DIGI LOOP - EE **

	if (etee>0.1) {

	  float pedoff=0;

	  if (isMC_) pedoff=0;
	  if (!isMC_) pedoff=200-pedg12;

	  
	  //cout << "EE DIGI printout:\n\n";
	 
	  //  cout << "rechit energy=" << ampli << " rechit et=" << etee << endl; 
	  //cout << "SAMPLE   ADC   GAIN\n";
	  
	  Int_t eedigiadc=0;
	  Int_t eedigigain=0;
	  char mytxtee[80];
	  
	  
	  EEDigiCollection::const_iterator digiItrEE= EEdigis->begin();
	  while(digiItrEE != EEdigis->end() && digiItrEE->id() != hitItr->id())
	    {
	      digiItrEE++;
	    }
	  
	  if (digiItrEE != EEdigis->end()) {
	    
	    EEDataFrame df(*digiItrEE);
	    for(int i=0; i<10;++i)
	      {
		eedigiadc=df.sample(i).adc();
		eedigiadc+=pedoff;
		eedigigain=df.sample(i).gainId();
		//		sprintf(mytxtee,"  %02d    %04d     %d",i+1,eedigiadc,eedigigain);
		//		cout << mytxtee << endl;
		

		ee_digi_01->Fill(i+0.5,eedigiadc);
		
		if (etee>0.5) ee_digi_05->Fill(i+0.5,eedigiadc);


		if (etee>0.1 && etee<=0.5) {

		  ee_digi_0105->Fill(i+0.5,eedigiadc);
		  ee_digi_0105_vs_time->Fill(time,i+0.5,float(eedigiadc));
		  ee_digi_0105_vs_bxtrain->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));

		  ee_digi_0105_vs_time_norm->Fill(time,i+0.5,1.0);
		  ee_digi_0105_vs_bxtrain_norm->Fill(bunchintrain+0.5,i+0.5,1.0);


		  if (abs(ee_eta)<2.0) {

		    ee_digi_0105_vs_time_eta15->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0105_vs_bxtrain_eta15->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0105_vs_time_norm_eta15->Fill(time,i+0.5,1.0);
		    ee_digi_0105_vs_bxtrain_norm_eta15->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(ee_eta)>=2.0 && abs(ee_eta)<2.5) {

		    ee_digi_0105_vs_time_eta20->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0105_vs_bxtrain_eta20->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0105_vs_time_norm_eta20->Fill(time,i+0.5,1.0);
		    ee_digi_0105_vs_bxtrain_norm_eta20->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(ee_eta)>=2.5) {

		    ee_digi_0105_vs_time_eta25->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0105_vs_bxtrain_eta25->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0105_vs_time_norm_eta25->Fill(time,i+0.5,1.0);
		    ee_digi_0105_vs_bxtrain_norm_eta25->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }






		}

		if (etee>0.5 && etee<=3.0) {

		  ee_digi_0530->Fill(i+0.5,eedigiadc);
		  ee_digi_0530_vs_time->Fill(time,i+0.5,float(eedigiadc));
		  ee_digi_0530_vs_bxtrain->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));

		  ee_digi_0530_vs_time_norm->Fill(time,i+0.5,1.0);
		  ee_digi_0530_vs_bxtrain_norm->Fill(bunchintrain+0.5,i+0.5,1.0);



		  if (abs(ee_eta)<2.0) {

		    ee_digi_0530_vs_time_eta15->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0530_vs_bxtrain_eta15->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0530_vs_time_norm_eta15->Fill(time,i+0.5,1.0);
		    ee_digi_0530_vs_bxtrain_norm_eta15->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(ee_eta)>=2.0 && abs(ee_eta)<2.5) {

		    ee_digi_0530_vs_time_eta20->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0530_vs_bxtrain_eta20->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0530_vs_time_norm_eta20->Fill(time,i+0.5,1.0);
		    ee_digi_0530_vs_bxtrain_norm_eta20->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }


		  if (abs(ee_eta)>=2.5) {

		    ee_digi_0530_vs_time_eta25->Fill(time,i+0.5,float(eedigiadc));
		    ee_digi_0530_vs_bxtrain_eta25->Fill(bunchintrain+0.5,i+0.5,float(eedigiadc));
		    
		    ee_digi_0530_vs_time_norm_eta25->Fill(time,i+0.5,1.0);
		    ee_digi_0530_vs_bxtrain_norm_eta25->Fill(bunchintrain+0.5,i+0.5,1.0);

		  }




		}


	      }
	    
	  }

	}
	





	   }

 

	   if (ampli>eemax && ee_kGood) {
	     eemax=ampli;
	     eemaxet=etee;
	     eeix=ix;      
	     eeiy=iy;
	     eeiz=side;
	     ee_eta=eta_ee;
	     ee_phi=phi_ee;
	     eetime=time;
	     eeflags=ebflag;
	     maxeehit=*hitItr;
	     
	   }
	   
	   if (etee>eemaxet2 && ee_kGood) {
	     eemax2=ampli;
	     eemaxet2=etee;
	     eeix2=ix;      
	     eeiy2=iy;
	     eeiz2=side;
	     ee_eta2=eta_ee;
	     ee_phi2=phi_ee;
	     eetime2=time;
	     eeflags2=ebflag;
	     maxeehit2=*hitItr;
	     
	   }



      Float_t ix_tmp=ix-0.5;
      if (side==-1) ix_tmp=ix_tmp+100;
      Float_t iy_tmp=iy-0.5;


      goodevent=true;

      if (!(ee_kWeird))  {

	//   cout << "EErechit" << endl;





	

        ee_rechitenergy_->Fill(ampli);



	rechiteta_all->Fill(eta_ee);
	rechiteta_etweight->Fill(eta_ee,etee);

	if (etee>1.0) {
	  rechiteta_gt1et->Fill(eta_ee);
	  rechiteta_etweight_gt1et->Fill(eta_ee,etee);

	  if (numgoodvtx<10) 	  rechiteta_gt1et_pu10->Fill(eta_ee);
	  if (numgoodvtx>=10 && numgoodvtx<20) 	  rechiteta_gt1et_pu20->Fill(eta_ee);
	  if (numgoodvtx>20) 	  rechiteta_gt1et_pu30->Fill(eta_ee);

	}





    
        if (fabs(eta_ee)<1.6)    ee_rechitenergy_16->Fill(ampli);
        if (fabs(eta_ee)>=1.6 && fabs(eta_ee)<1.8)    ee_rechitenergy_18->Fill(ampli);
        if (fabs(eta_ee)>=1.8 && fabs(eta_ee)<2.0)    ee_rechitenergy_20->Fill(ampli);
        if (fabs(eta_ee)>=2.0 && fabs(eta_ee)<2.2)    ee_rechitenergy_22->Fill(ampli);
        if (fabs(eta_ee)>=2.2 && fabs(eta_ee)<2.4)    ee_rechitenergy_24->Fill(ampli);
        if (fabs(eta_ee)>=2.4 && fabs(eta_ee)<2.6)    ee_rechitenergy_26->Fill(ampli);
        if (fabs(eta_ee)>=2.6 && fabs(eta_ee)<2.8)    ee_rechitenergy_28->Fill(ampli);
        if (fabs(eta_ee)>=2.8)    ee_rechitenergy_30->Fill(ampli);

	ee_rechitet_->Fill(etee);

	if (fabs(eta_ee)<1.6)    ee_rechitet_16->Fill(etee);
        if (fabs(eta_ee)>=1.6 && fabs(eta_ee)<1.8)    ee_rechitet_18->Fill(etee);
        if (fabs(eta_ee)>=1.8 && fabs(eta_ee)<2.0)    ee_rechitet_20->Fill(etee);
        if (fabs(eta_ee)>=2.0 && fabs(eta_ee)<2.2)    ee_rechitet_22->Fill(etee);
        if (fabs(eta_ee)>=2.2 && fabs(eta_ee)<2.4)    ee_rechitet_24->Fill(etee);
        if (fabs(eta_ee)>=2.4 && fabs(eta_ee)<2.6)    ee_rechitet_26->Fill(etee);
        if (fabs(eta_ee)>=2.6 && fabs(eta_ee)<2.8)    ee_rechitet_28->Fill(etee);
        if (fabs(eta_ee)>=2.8)    ee_rechitet_30->Fill(etee);


	if (fabs(eta_ee)<2.0) ee_rechitetvspu_20->Fill(float(numgoodvtx)-0.5,etee);
	if (fabs(eta_ee)>=2.0 && fabs(eta_ee)<2.5) ee_rechitetvspu_25->Fill(float(numgoodvtx)-0.5,etee);
	if (fabs(eta_ee)>=2.5) ee_rechitetvspu_30->Fill(float(numgoodvtx)-0.5,etee);
 

	//     ee_rechiten_vs_sm->Fill(dee-0.5,ampli);
	// ee_rechitet_vs_sm->Fill(dee-0.5,etee);

        if (eeRecHits->size()>500)  ee_rechitenergy_notypeb_->Fill(ampli);

	if (side==1) {

           eep_rechiten_vs_eta->Fill(fabs(eta_ee),ampli);
           eep_rechiten_vs_phi->Fill(phi_ee,ampli);

           eep_rechitet_vs_eta->Fill(fabs(eta_ee),etee);
           eep_rechitet_vs_phi->Fill(phi_ee,etee);

        }

        if (side==-1) {

           eem_rechiten_vs_eta->Fill(fabs(eta_ee),ampli);
           eem_rechiten_vs_phi->Fill(phi_ee,ampli);
           eem_rechitet_vs_eta->Fill(fabs(eta_ee),etee);
           eem_rechitet_vs_phi->Fill(phi_ee,etee);
 
        }


         eeocc->Fill(ix_tmp,iy_tmp,1.);
         eeoccen->Fill(ix_tmp,iy_tmp,ampli);
         eeoccet->Fill(ix_tmp,iy_tmp,etee);

         if (etee>1.0) {

            eeoccgt1et->Fill(ix_tmp,iy_tmp,1.);
            eeoccetgt1et->Fill(ix_tmp,iy_tmp,etee);
           
         }

        if (ampli>1.0) {

            eeoccgt1->Fill(ix_tmp,iy_tmp,1.);
            eeoccengt1->Fill(ix_tmp,iy_tmp,ampli);
           
         }


 	   int sevlev=severity->severityLevel(det,*eeRecHits);


	   bool badsc1=false;
	   bool badsc2=false;

	   if (ix>=21 && ix<=25 && iy>=21 && iy<=25 && iz==-1) badsc1=true;
	   if (ix>=46 && ix<=50 && iy>=96 && iy<=100 && iz==+1) badsc2=true;





	   if (badsc1) {

	     badsc1_et+=etee;
	     if (ampli>1000.0 && ee_kGood==0) badsc1_hits++;

	   }

	   if (badsc2) {

	     badsc2_et+=etee;
	     if (ampli>1000.0 && ee_kGood==0) badsc2_hits++;

	   }






   

      eehits++;
      if (side==1) eephits++;
      if (side==-1) eemhits++;

      
	   }

  }
  

  // end of ee loop

  

   // tp and calotower matching algo

    numcalotower=0;
    numtrigprim=0;
    numcalotower_eb=0;
    numtrigprim_eb=0;
    numcalotower_ee=0;
    numtrigprim_ee=0;
    numrechittower=0;
    numrechittower_eb=0;
    numrechittower_ee=0;

    numcalotower_notmatched1=0;
    numtrigprim_notmatched1=0;
    numcalotower_notmatched1_eb=0;
    numtrigprim_notmatched1_eb=0;
    numcalotower_notmatched1_ee=0;
    numtrigprim_notmatched1_ee=0;
 

    numcalotower_notmatched2=0;
    numtrigprim_notmatched2=0;
    numcalotower_notmatched2_eb=0;
    numtrigprim_notmatched2_eb=0;
    numcalotower_notmatched2_ee=0;
    numtrigprim_notmatched2_ee=0;
 
    numcalotower_matched=0;
    numtrigprim_matched=0;
    numcalotower_matched_eb=0;
    numtrigprim_matched_eb=0;
    numcalotower_matched_ee=0;
    numtrigprim_matched_ee=0;
 
    numcalotower_matched2=0;
    numtrigprim_matched2=0;
    numcalotower_matched_eb2=0;
    numtrigprim_matched_eb2=0;
    numcalotower_matched_ee2=0;
    numtrigprim_matched_ee2=0;
 
    numrechittower_matched_rh=0;
    numtrigprim_matched_rh=0;
    numrechittower_matched_eb_rh=0;
    numtrigprim_matched_eb_rh=0;
    numrechittower_matched_ee_rh=0;
    numtrigprim_matched_ee_rh=0;
 
    numrechittower_matched_rh2=0;
    numtrigprim_matched_rh2=0;
    numrechittower_matched_eb_rh2=0;
    numtrigprim_matched_eb_rh2=0;
    numrechittower_matched_ee_rh2=0;
    numtrigprim_matched_ee_rh2=0;
 


    sumcalotower=0;
    sumtrigprim=0;
    sumcalotower_eb=0;
    sumtrigprim_eb=0;
    sumcalotower_ee=0;
    sumtrigprim_ee=0;


    sumtrigprim_samp2_eb=0;
    sumtrigprim_samp2_ee=0;

    sumrechittower=0;
    sumrechittower_eb=0;
    sumrechittower_ee=0;


    sumcalotower_notmatched1=0;
    sumtrigprim_notmatched1=0;
    sumcalotower_notmatched1_eb=0;
    sumtrigprim_notmatched1_eb=0;
    sumcalotower_notmatched1_ee=0;
    sumtrigprim_notmatched1_ee=0;
 

    sumcalotower_notmatched2=0;
    sumtrigprim_notmatched2=0;
    sumcalotower_notmatched2_eb=0;
    sumtrigprim_notmatched2_eb=0;
    sumcalotower_notmatched2_ee=0;
    sumtrigprim_notmatched2_ee=0;
 
    sumcalotower_matched=0;
    sumtrigprim_matched=0;
    sumcalotower_matched_eb=0;
    sumtrigprim_matched_eb=0;
    sumcalotower_matched_ee=0;
    sumtrigprim_matched_ee=0;


    sumcalotower_matched2=0;
    sumtrigprim_matched2=0;
    sumcalotower_matched_eb2=0;
    sumtrigprim_matched_eb2=0;
    sumcalotower_matched_ee2=0;
    sumtrigprim_matched_ee2=0;
 
    sumrechittower_matched_rh=0;
    sumtrigprim_matched_rh=0;
    sumrechittower_matched_eb_rh=0;
    sumtrigprim_matched_eb_rh=0;
    sumrechittower_matched_ee_rh=0;
    sumtrigprim_matched_ee_rh=0;


    sumrechittower_matched_rh2=0;
    sumtrigprim_matched_rh2=0;
    sumrechittower_matched_eb_rh2=0;
    sumtrigprim_matched_eb_rh2=0;
    sumrechittower_matched_ee_rh2=0;
    sumtrigprim_matched_ee_rh2=0;



    numtrigprim_ee24=0;
    sumtrigprim_ee24=0;
    numrechittower_ee24=0;
    sumrechittower_ee24=0;

    numtrigprim_matched_ee24_rh=0;
    sumtrigprim_matched_ee24_rh=0;

    numrechittower_matched_ee24_rh=0;
    sumrechittower_matched_ee24_rh=0;




    for (Int_t i=0;i<72;i++) {
      for (Int_t j=0;j<57;j++) {

       int iphi_tmp=i+1;
       int ieta_tmp=-28+j;
       if (ieta_tmp==0) continue;


       float tp_et=tp_map->GetBinContent(i+1,j+1);
       float tp_et2=tp_map2->GetBinContent(i+1,j+1);
       float tp_et_samp2=tp_map_samp2->GetBinContent(i+1,j+1);
       float ct_et=calotower_map->GetBinContent(i+1,j+1);
       float rh_et=rechit_map->GetBinContent(i+1,j+1);


       if (rh_et>1.0) {

	 // cout << "ieta,iphi=" << ieta_tmp << " , " << iphi_tmp << " RH ET=" << rh_et << endl;  


	 rechitmap_cumul->Fill(iphi_tmp,ieta_tmp);
	 numrechittower++;
	 sumrechittower+=rh_et;

	 if (fabs(ieta_tmp)<=17) {
	   numrechittower_eb++;
	   sumrechittower_eb+=rh_et;
	 }


	 if (fabs(ieta_tmp)>17) {
	   numrechittower_ee++;
	   sumrechittower_ee+=rh_et;

	   rh_et_all_ee->Fill(rh_et);

	   if (fabs(ieta_tmp)<25) {
	     numrechittower_ee24++;
	     sumrechittower_ee24+=rh_et;
	     rh_et_all_ee24->Fill(rh_et);

	   }

	 }

       }


       if (tp_et>0.01) {

	 // cout << "ieta,iphi=" << ieta_tmp << " , " << iphi_tmp << " TP ET=" << tp_et << endl;  


	 tpmap_cumul->Fill(iphi_tmp,ieta_tmp);

	 numtrigprim++;
	 sumtrigprim+=tp_et;

	 tp_et_all->Fill(tp_et);
	 tp_et_vs_ieta_all->Fill(ieta_tmp,tp_et);

	 if (fabs(ieta_tmp)<=17) {
	   numtrigprim_eb++;
	   sumtrigprim_eb+=tp_et;

	   tp_et_all_eb->Fill(tp_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numtrigprim_ee++;
	   sumtrigprim_ee+=tp_et;
	   tp_et_all_ee->Fill(tp_et);

	   if (fabs(ieta_tmp)<25) {
	     numtrigprim_ee24++;
	     sumtrigprim_ee24+=tp_et;
	     tp_et_all_ee24->Fill(tp_et);
	   }


	   if (fabs(ieta_tmp)>=18 && fabs(ieta_tmp)<=20)  sumtrigprim_ee_eta1820+=tp_et;
	   if (fabs(ieta_tmp)>=21 && fabs(ieta_tmp)<=22)  sumtrigprim_ee_eta2122+=tp_et;
	   if (fabs(ieta_tmp)>=23 && fabs(ieta_tmp)<=24)  sumtrigprim_ee_eta2324+=tp_et;
	   if (fabs(ieta_tmp)>=25 && fabs(ieta_tmp)<=26)  sumtrigprim_ee_eta2526+=tp_et;
	   if (fabs(ieta_tmp)>=27 && fabs(ieta_tmp)<=28)  sumtrigprim_ee_eta2728+=tp_et;
	   


	 }

       }



      if (tp_et_samp2>0.01) {

	 // cout << "ieta,iphi=" << ieta_tmp << " , " << iphi_tmp << " TP ET=" << tp_et << endl;  


	 if (fabs(ieta_tmp)<=17) {

	   sumtrigprim_samp2_eb+=tp_et_samp2;

	 }

	 if (fabs(ieta_tmp)>17) {

	   sumtrigprim_samp2_ee+=tp_et_samp2;

	 }

      }


      if (ct_et>0.01) {

	 calotower_map_cumul->Fill(iphi_tmp,ieta_tmp);

	 numcalotower++;
	 sumcalotower+=ct_et;

	 ct_et_all->Fill(ct_et);
	 ct_et_vs_ieta_all->Fill(ieta_tmp,ct_et);


	 if (fabs(ieta_tmp)<=17) {
	   numcalotower_eb++;
	   sumcalotower_eb+=ct_et;
	   ct_et_all_eb->Fill(ct_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numcalotower_ee++;
	   sumcalotower_ee+=ct_et;
	   ct_et_all_ee->Fill(ct_et);

	 }

       }




      // matched calotower and TP
      if (ct_et>0.01 && tp_et>0.01) {

	 tower_map_cumul_matched->Fill(iphi_tmp,ieta_tmp);

	 numcalotower_matched++;
	 sumcalotower_matched+=ct_et;


	 numtrigprim_matched++;
	 sumtrigprim_matched+=tp_et;



	 tp_et_matched->Fill(tp_et);
	 tp_et_vs_ieta_matched->Fill(ieta_tmp,tp_et);

	 ct_et_matched->Fill(ct_et);
	 ct_et_vs_ieta_matched->Fill(ieta_tmp,ct_et);


	 ctminustp_et_matched->Fill(ct_et-tp_et);
	 ctminustp_et_vs_ieta_matched->Fill(ieta_tmp,ct_et-tp_et);
	 tp_et_vs_ct_et_matched->Fill(ct_et,tp_et);



	 if (fabs(ieta_tmp)<=17) {
	   numcalotower_matched_eb++;
	   sumcalotower_matched_eb+=ct_et;
	   numtrigprim_matched_eb++;
	   sumtrigprim_matched_eb+=tp_et;
	   tp_et_matched_eb->Fill(tp_et);
	   ct_et_matched_eb->Fill(ct_et);
	   ctminustp_et_matched_eb->Fill(ct_et-tp_et);
	   tp_et_vs_ct_et_matched_eb->Fill(ct_et,tp_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numcalotower_matched_ee++;
	   sumcalotower_matched_ee+=ct_et;
	   numtrigprim_matched_ee++;
	   sumtrigprim_matched_ee+=tp_et;
	   tp_et_matched_ee->Fill(tp_et);
	   ct_et_matched_ee->Fill(ct_et);
	   ctminustp_et_matched_ee->Fill(ct_et-tp_et);
	   tp_et_vs_ct_et_matched_ee->Fill(ct_et,tp_et);

	 }

       }



      // matched calotower and spike-cleaned TP
      if (ct_et>0.01 && tp_et2>0.01) {

	 tower_map_cumul_matched2->Fill(iphi_tmp,ieta_tmp);

	 numcalotower_matched2++;
	 sumcalotower_matched2+=ct_et;


	 numtrigprim_matched2++;
	 sumtrigprim_matched2+=tp_et2;



	 tp_et_matched2->Fill(tp_et2);
	 tp_et_vs_ieta_matched2->Fill(ieta_tmp,tp_et2);

	 ct_et_matched2->Fill(ct_et);
	 ct_et_vs_ieta_matched2->Fill(ieta_tmp,ct_et);


	 ctminustp_et_matched2->Fill(ct_et-tp_et2);
	 ctminustp_et_vs_ieta_matched2->Fill(ieta_tmp,ct_et-tp_et2);
	 tp_et_vs_ct_et_matched2->Fill(ct_et,tp_et2);



	 if (fabs(ieta_tmp)<=17) {

	   numcalotower_matched_eb2++;
	   sumcalotower_matched_eb2+=ct_et;
	   numtrigprim_matched_eb2++;
	   sumtrigprim_matched_eb+=tp_et2;
	   tp_et_matched_eb2->Fill(tp_et2);
	   ct_et_matched_eb2->Fill(ct_et);
	   ctminustp_et_matched_eb2->Fill(ct_et-tp_et2);
	   tp_et_vs_ct_et_matched_eb2->Fill(ct_et,tp_et2);

	 }

	 if (fabs(ieta_tmp)>17) {
       
	   numcalotower_matched_ee2++;
	   sumcalotower_matched_ee2+=ct_et;
	   numtrigprim_matched_ee2++;
	   sumtrigprim_matched_ee2+=tp_et2;
	   tp_et_matched_ee2->Fill(tp_et2);
	   ct_et_matched_ee2->Fill(ct_et);
	   ctminustp_et_matched_ee2->Fill(ct_et-tp_et2);
	   tp_et_vs_ct_et_matched_ee2->Fill(ct_et,tp_et2);

	 }

       }




      // matched rechit tower and TP
      if (rh_et>1.0 && tp_et>0.01) {

	if (rh_et>3 && tp_et<1) {
	  //  cout << "ieta,iphi=" << ieta_tmp << "," << iphi_tmp << " rh_et=" << rh_et << " tp_et=" << tp_et << endl;

	}

	 tower_map_cumul_matched_rh->Fill(iphi_tmp,ieta_tmp);

	 numrechittower_matched_rh++;
	 sumrechittower_matched_rh+=rh_et;


	 numtrigprim_matched_rh++;
	 sumtrigprim_matched_rh+=tp_et;



	 tp_et_matched_rh->Fill(tp_et);
	 tp_et_vs_ieta_matched_rh->Fill(ieta_tmp,tp_et);

	 rh_et_matched_rh->Fill(rh_et);
	 rh_et_vs_ieta_matched_rh->Fill(ieta_tmp,rh_et);


	 rhminustp_et_matched_rh->Fill(rh_et-tp_et);
	 rhminustp_et_vs_ieta_matched_rh->Fill(ieta_tmp,rh_et-tp_et);
	 tp_et_vs_rh_et_matched_rh->Fill(rh_et,tp_et);



	 if (fabs(ieta_tmp)<=17) {
	   numrechittower_matched_eb_rh++;
	   sumrechittower_matched_eb_rh+=rh_et;
	   numtrigprim_matched_eb_rh++;
	   sumtrigprim_matched_eb_rh+=tp_et;
	   tp_et_matched_eb_rh->Fill(tp_et);
	   rh_et_matched_eb_rh->Fill(rh_et);
	   rhminustp_et_matched_eb_rh->Fill(rh_et-tp_et);
	   tp_et_vs_rh_et_matched_eb_rh->Fill(rh_et,tp_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numrechittower_matched_ee_rh++;
	   sumrechittower_matched_ee_rh+=rh_et;
	   numtrigprim_matched_ee_rh++;
	   sumtrigprim_matched_ee_rh+=tp_et;
	   tp_et_matched_ee_rh->Fill(tp_et);
	   rh_et_matched_ee_rh->Fill(rh_et);
	   rhminustp_et_matched_ee_rh->Fill(rh_et-tp_et);
	   tp_et_vs_rh_et_matched_ee_rh->Fill(rh_et,tp_et);


	   if (fabs(ieta_tmp)<25) {
	     numtrigprim_matched_ee24_rh++;
	     sumtrigprim_matched_ee24_rh+=tp_et;
	     tp_et_matched_ee24_rh->Fill(tp_et);
	     numrechittower_matched_ee24_rh++;
	     sumrechittower_matched_ee24_rh+=rh_et;
	     rh_et_matched_ee24_rh->Fill(rh_et);
	   }


	   if(fabs(ieta_tmp)>=18 && fabs(ieta_tmp)<=20) {
	     tp_et_vs_rh_et_matched_ee_ieta1820_rh->Fill(rh_et,tp_et);
	     sumtrigprim_ee_matched_eta1820+=tp_et;
	     sumrechittower_ee_matched_eta1820+=rh_et;
	   }

	   if(fabs(ieta_tmp)>=21 && fabs(ieta_tmp)<=22) {
	     tp_et_vs_rh_et_matched_ee_ieta2122_rh->Fill(rh_et,tp_et);
	     sumtrigprim_ee_matched_eta2122+=tp_et;
	     sumrechittower_ee_matched_eta2122+=rh_et;
	   }

	   if(fabs(ieta_tmp)>=23 && fabs(ieta_tmp)<=24) {
	     tp_et_vs_rh_et_matched_ee_ieta2324_rh->Fill(rh_et,tp_et);
	     sumtrigprim_ee_matched_eta2324+=tp_et;
	     sumrechittower_ee_matched_eta2324+=rh_et;
	   }

	   if(fabs(ieta_tmp)>=25 && fabs(ieta_tmp)<=26) {
	     tp_et_vs_rh_et_matched_ee_ieta2526_rh->Fill(rh_et,tp_et);
	     sumtrigprim_ee_matched_eta2526+=tp_et;
	     sumrechittower_ee_matched_eta2526+=rh_et;
	   }

	   if(fabs(ieta_tmp)>=27 && fabs(ieta_tmp)<=28) {
	     tp_et_vs_rh_et_matched_ee_ieta2728_rh->Fill(rh_et,tp_et);
	     sumtrigprim_ee_matched_eta2728+=tp_et;
	     sumrechittower_ee_matched_eta2728+=rh_et;
	   }


	 if (fabs(ieta_tmp)>17) {


// 	   cout << "Matched tower:  ieta,iphi = " << ieta_tmp << "," << iphi_tmp << 
// 	    " TP E_t = " << tp_et << " RH E_t = " << rh_et << endl;


	   ietamatched->Fill(ieta_tmp);
	   ietatotal->Fill(ieta_tmp);



// 	   cout << " Rechits: " << endl;

	   Int_t iztmp=1;
	   if (ieta_tmp<0) iztmp=-1;


	   for (Int_t ix=0;ix<100;ix++) {
	     for (Int_t iy=0;iy<100;iy++) {


	       int towerieta=ietamapTT->GetBinContent(ix+1,iy+1);
	       int toweriphi=iphimapTT->GetBinContent(ix+1,iy+1);


	       if (iztmp<0) towerieta*=-1;
	       
	       if (towerieta==ieta_tmp && toweriphi==iphi_tmp) {
	       

		 float rechit_et=rechit_map_ixiy->GetBinContent(ix+1+100*(iztmp>0),iy+1);
		 float rechit_time=rechit_map_time_ixiy->GetBinContent(ix+1+100*(iztmp>0),iy+1);

		 if (rechit_et!=0) {

		   matched_timevset->Fill(rechit_time,rechit_et);
		   matched_timevsieta->Fill(rechit_time,fabs(ieta_tmp));
		   matched_timevsieta_etweighted->Fill(rechit_time,fabs(ieta_tmp),rechit_et);
		   matched_timevsetvsieta->Fill(rechit_time,rechit_et,fabs(ieta_tmp));

		 }

	// 	 cout << "   ix=" << ix+1 << " iy=" << iy+1 << " iz=" << iztmp
// 		      << " et=" << rechit_et << " time=" << rechit_time << endl;
		 
	       }
	       
	     }
	   }

	 }

	 }

       }

      // unmatched rechit tower and TP
      if (rh_et>1.0 && tp_et<0.01) {


	 if (fabs(ieta_tmp)>17) {

// 	   cout << "Unmatched tower:  ieta,iphi = " << ieta_tmp << "," << iphi_tmp << 
// 	    " TP E_t = " << tp_et << " RH E_t = " << rh_et << endl;

// 	   cout << " Rechits: " << endl;

	   Int_t iztmp=1;
	   if (ieta_tmp<0) iztmp=-1;

	   ietatotal->Fill(ieta_tmp);

	   for (Int_t ix=0;ix<100;ix++) {
	     for (Int_t iy=0;iy<100;iy++) {


	       int towerieta=ietamapTT->GetBinContent(ix+1,iy+1);
	       int toweriphi=iphimapTT->GetBinContent(ix+1,iy+1);


	       if (iztmp<0) towerieta*=-1;
	       
	       if (towerieta==ieta_tmp && toweriphi==iphi_tmp) {
	       

		 float rechit_et=rechit_map_ixiy->GetBinContent(ix+1+100*(iztmp>0),iy+1);
		 float rechit_time=rechit_map_time_ixiy->GetBinContent(ix+1+100*(iztmp>0),iy+1);

		 if (rechit_et!=0) {

		   unmatched_timevset->Fill(rechit_time,rechit_et);
		   unmatched_timevsieta->Fill(rechit_time,fabs(ieta_tmp));
		   unmatched_timevsieta_etweighted->Fill(rechit_time,fabs(ieta_tmp),rechit_et);
		   unmatched_timevsetvsieta->Fill(rechit_time,rechit_et,fabs(ieta_tmp));

		 }

// 		 cout << "   ix=" << ix+1 << " iy=" << iy+1 << " iz=" << iztmp
// 		      << " et=" << rechit_et << " time=" << rechit_time << endl;
		 
	       }
	       
	     }
	   }
	 }
      }


      // matched rechit tower and spike-cleaned TP
     if (rh_et>1.0 && tp_et2>0.01) {

	 tower_map_cumul_matched_rh2->Fill(iphi_tmp,ieta_tmp);

	 numrechittower_matched_rh2++;
	 sumrechittower_matched_rh2+=rh_et;


	 numtrigprim_matched_rh2++;
	 sumtrigprim_matched_rh2+=tp_et2;



	 tp_et_matched_rh2->Fill(tp_et2);
	 tp_et_vs_ieta_matched_rh2->Fill(ieta_tmp,tp_et2);

	 rh_et_matched_rh2->Fill(rh_et);
	 rh_et_vs_ieta_matched_rh2->Fill(ieta_tmp,rh_et);


	 rhminustp_et_matched_rh2->Fill(rh_et-tp_et2);
	 rhminustp_et_vs_ieta_matched_rh2->Fill(ieta_tmp,rh_et-tp_et2);
	 tp_et_vs_rh_et_matched_rh2->Fill(rh_et,tp_et2);



	 if (fabs(ieta_tmp)<=17) {
	   numrechittower_matched_eb_rh2++;
	   sumrechittower_matched_eb_rh2+=rh_et;
	   numtrigprim_matched_eb_rh2++;
	   sumtrigprim_matched_eb_rh2+=tp_et2;
	   tp_et_matched_eb_rh2->Fill(tp_et2);
	   rh_et_matched_eb_rh2->Fill(rh_et);
	   rhminustp_et_matched_eb_rh2->Fill(rh_et-tp_et2);
	   tp_et_vs_rh_et_matched_eb_rh2->Fill(rh_et,tp_et2);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numrechittower_matched_ee_rh2++;
	   sumrechittower_matched_ee_rh2+=rh_et;
	   numtrigprim_matched_ee_rh2++;
	   sumtrigprim_matched_ee_rh2+=tp_et2;
	   tp_et_matched_ee_rh2->Fill(tp_et2);
	   rh_et_matched_ee_rh2->Fill(rh_et);
	   rhminustp_et_matched_ee_rh2->Fill(rh_et-tp_et2);
	   tp_et_vs_rh_et_matched_ee_rh2->Fill(rh_et,tp_et2);

	 }

       }






     // unmatched calotower and TP
     if (ct_et>0.01 && tp_et<0.01) {

       tower_map_cumul_unmatched1->Fill(iphi_tmp,ieta_tmp);

	 numcalotower_notmatched1++;
	 sumcalotower_notmatched1+=ct_et;


	 numtrigprim_notmatched1++;
	 sumtrigprim_notmatched1+=tp_et;


	 ct_et_unmatched->Fill(ct_et);
	 ct_et_vs_ieta_unmatched->Fill(ieta_tmp,ct_et);



	 if (fabs(ieta_tmp)<=17) {
	   numcalotower_notmatched1_eb++;
	   sumcalotower_notmatched1_eb+=ct_et;
	   numtrigprim_notmatched1_eb++;
	   sumtrigprim_notmatched1_eb+=tp_et;
	   ct_et_unmatched_eb->Fill(ct_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numcalotower_notmatched1_ee++;
	   sumcalotower_notmatched1_ee+=ct_et;
	   numtrigprim_notmatched1_ee++;
	   sumtrigprim_notmatched1_ee+=tp_et;
	   ct_et_unmatched_ee->Fill(ct_et);

	 }

       }


     // unmatched calotower and TP
     if (ct_et<0.01 && tp_et>0.01) {


       tower_map_cumul_unmatched2->Fill(iphi_tmp,ieta_tmp);

	 numcalotower_notmatched2++;
	 sumcalotower_notmatched2+=ct_et;


	 numtrigprim_notmatched2++;
	 sumtrigprim_notmatched2+=tp_et;


	 tp_et_unmatched->Fill(tp_et);
	 tp_et_vs_ieta_unmatched->Fill(ieta_tmp,tp_et);




	 if (fabs(ieta_tmp)<=17) {
	   numcalotower_notmatched2_eb++;
	   sumcalotower_notmatched2_eb+=ct_et;
	   numtrigprim_notmatched2_eb++;
	   sumtrigprim_notmatched2_eb+=tp_et;
	   tp_et_unmatched_eb->Fill(tp_et);

	 }

	 if (fabs(ieta_tmp)>17) {
	   numcalotower_notmatched2_ee++;
	   sumcalotower_notmatched2_ee+=ct_et;
	   numtrigprim_notmatched2_ee++;
	   sumtrigprim_notmatched2_ee+=tp_et;
	   tp_et_unmatched_ee->Fill(tp_et);

	 }

       }




     }
   }


//    cout << "MATCHING SUMMARY: " << endl;

//    cout << " TOTAL:  numcalotower=" << numcalotower << " numcalotower_eb=" << numcalotower_eb << " numcalotower_ee=" << numcalotower_ee << " sumct=" << sumcalotower << " sumct_eb=" << sumcalotower_eb << " sumct_ee=" << sumcalotower_ee << endl;


//    cout << " TOTAL:  numtrigprim=" << numtrigprim << " numtrigprim_eb=" << numtrigprim_eb << " numtrigprim_ee=" << numtrigprim_ee << " sumtp=" << sumtrigprim << " sumtp_eb=" << sumtrigprim_eb << " sumtp_ee=" << sumtrigprim_ee << endl;





//    cout << " MATCHED:  numcalotower=" << numcalotower_matched << " numcalotower_eb=" << numcalotower_matched_eb << " numcalotower_ee=" << numcalotower_matched_ee << " sumct=" << sumcalotower_matched << " sumct_eb=" << sumcalotower_matched_eb << " sumct_ee=" << sumcalotower_matched_ee << endl;


//    cout << " MATCHED:  numtrigprim=" << numtrigprim_matched << " numtrigprim_eb=" << numtrigprim_matched_eb << " numtrigprim_ee=" << numtrigprim_matched_ee << " sumtp=" << sumtrigprim_matched << " sumtp_eb=" << sumtrigprim_matched_eb << " sumtp_ee=" << sumtrigprim_matched_ee << endl;



//    cout << " NOTMATCHED1:  numcalotower=" << numcalotower_notmatched1 << " numcalotower_eb=" << numcalotower_notmatched1_eb << " numcalotower_ee=" << numcalotower_notmatched1_ee << " sumct=" << sumcalotower_notmatched1 << " sumct_eb=" << sumcalotower_notmatched1_eb << " sumct_ee=" << sumcalotower_notmatched1_ee << endl;


//    cout << " NOTMATCHED1:  numtrigprim=" << numtrigprim_notmatched1 << " numtrigprim_eb=" << numtrigprim_notmatched1_eb << " numtrigprim_ee=" << numtrigprim_notmatched1_ee << " sumtp=" << sumtrigprim_notmatched1 << " sumtp_eb=" << sumtrigprim_notmatched1_eb << " sumtp_ee=" << sumtrigprim_notmatched1_ee << endl;


//    cout << " NOTMATCHED2:  numcalotower=" << numcalotower_notmatched2 << " numcalotower_eb=" << numcalotower_notmatched2_eb << " numcalotower_ee=" << numcalotower_notmatched2_ee << " sumct=" << sumcalotower_notmatched2 << " sumct_eb=" << sumcalotower_notmatched2_eb << " sumct_ee=" << sumcalotower_notmatched2_ee << endl;


//    cout << " NOTMATCHED2:  numtrigprim=" << numtrigprim_notmatched2 << " numtrigprim_eb=" << numtrigprim_notmatched2_eb << " numtrigprim_ee=" << numtrigprim_notmatched2_ee << " sumtp=" << sumtrigprim_notmatched2 << " sumtp_eb=" << sumtrigprim_notmatched2_eb << " sumtp_ee=" << sumtrigprim_notmatched2_ee << endl;







   // end of tp and calotower matching algo







  // Rechit flags:

  if (maxebhit.recoFlag()==EcalRecHit::kGood) ebflag_kgood=1;
  if (maxebhit.recoFlag()==EcalRecHit::kPoorReco) ebflag_kpoorreco=1;
  if (maxebhit.recoFlag()==EcalRecHit::kOutOfTime) ebflag_koutoftime=1;
  // if (maxebhit.recoFlag()==EcalRecHit::kFake) ebflag_kfake=1;



  // some legacy EE code - not important    

 if (pTopology.isValid() && pG.isValid()) { 

    const CaloTopology *topology=pTopology.product();




   

    const reco::BasicCluster *cluster=0;
    float e3x3=0; //EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),-1,1,-1,1);


    float e5x5=0; //EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),-2,2,-2,2);

    eb_e9=0; //e3x3;
    eb_e25=0; //e5x5;

   
    

    float e4x1=0;

    //    e4x1+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),-1,-1,0,0);

    //    e4x1+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),1,1,0,0);

    //    e4x1+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),0,0,-1,-1);

    //    e4x1+=EcalClusterTools::matrixEnergy(*cluster,ebRecHits,topology,maxebhit.id(),0,0,+1,+1);




   
    if (e3x3>0) eb_r9=ebmax/e3x3;

    if (ebmax!=0) eb_r4=e4x1/ebmax;




    const reco::BasicCluster *cluster2=0;
     float e3x32=0; //EcalClusterTools::matrixEnergy(*cluster2,ebRecHits,topology,maxebhit2.id(),-1,1,-1,1);

    float e4x12=0;
 
    //    e4x12+=EcalClusterTools::matrixEnergy(*cluster2,ebRecHits,topology,maxebhit2.id(),-1,-1,0,0);

    //e4x12+=EcalClusterTools::matrixEnergy(*cluster2,ebRecHits,topology,maxebhit2.id(),1,1,0,0);

    //e4x12+=EcalClusterTools::matrixEnergy(*cluster2,ebRecHits,topology,maxebhit2.id(),0,0,-1,-1);

    // e4x12+=EcalClusterTools::matrixEnergy(*cluster2,ebRecHits,topology,maxebhit2.id(),0,0,+1,+1);


   
    if (e3x32>0) eb_r92=ebmax2/e3x32;

    if (ebmax2!=0) eb_r42=e4x12/ebmax2;




    const reco::BasicCluster *clusteree=0;
    float e3x3ee=0; //EcalClusterTools::matrixEnergy(*clusteree,eeRecHits,topology,maxeehit.id(),-1,1,-1,1);


   
    if (e3x3ee>0) ee_r9=eemax/e3x3ee;


    const reco::BasicCluster *clusteree2=0;
    float e3x3ee2=0; //EcalClusterTools::matrixEnergy(*clusteree2,eeRecHits,topology,maxeehit2.id(),-1,1,-1,1);


   
    if (e3x3ee2>0) ee_r92=eemax2/e3x3ee2;





  
 }




  dataTree_->Fill();

  //   }

  /*


 // cout << "begin jets" << endl;


  edm::Handle<PFJetCollection> jets;
  PFJetCollection::const_iterator i_jet;
  event.getByLabel (jets_,jets);
  std::vector<reco::PFCandidatePtr> PFJetPart;
  
  histo_event_->Fill(event.id().event());





  //  const JetCorrector* corrector = JetCorrector::getJetCorrector(jetcorr_,iSetup);
    

  int njet(0);
  if (jets->size()>0)
    {
      for (i_jet = jets->begin();i_jet != jets->end(); i_jet++)
       {    
              	
	   PFJetPart  = i_jet->getPFConstituents();
	  
	   //       scale=corrector->correction(i_jet->p4());
           scale=1.0;
           ncr_ = PFJetPart.size() ;

	   energy_pf_     = 0;
	   energyc_pf_    = 0;
	   energyn_pf_    = 0;
	   energyg_pf_    = 0;
	    
           energy_ecal    = 0;
           energy_hcal    = 0;

	   	   		  
	  for(Int_t j=0;j<PFJetPart.size();j++){
	    
	    PFCandidate::ParticleType type = PFJetPart[j]->particleId();
            int tpid  = PFJetPart[j]->translateTypeToPdgId(type);
	       
	    bool chargedp   = false;
	    bool neutralp   = false;
	    bool gammap     = false;
	    bool electronsp = false;
	   	      	    
	    chargedp   = TMath::Abs(tpid)==211 || TMath::Abs(tpid) == 321 || TMath::Abs(tpid)  == 2212 || TMath::Abs(tpid) == 3122|| TMath::Abs(tpid) == 3312|| TMath::Abs(tpid) == 3112 || TMath::Abs(tpid) == 3222;
	    neutralp   = TMath::Abs(tpid)==130 || TMath::Abs(tpid) == 2112 || TMath::Abs(tpid) == 310  || TMath::Abs(tpid) == 3322;
	    gammap     = TMath::Abs(tpid)==22;
	    electronsp = TMath::Abs(tpid)==11;
	    		    	   	  		   energy_pf_  = energy_pf_+ PFJetPart[j]->p4().energy(); 	 	  
	    if(chargedp)           energyc_pf_ = energyc_pf_+PFJetPart[j]->p4().energy();
	    if(neutralp)           energyn_pf_ = energyn_pf_+PFJetPart[j]->p4().energy();
	    if(gammap||electronsp) energyg_pf_ = energyg_pf_+PFJetPart[j]->p4().energy();
             //	    	   	    
            energy_ecal = energy_ecal + PFJetPart[j]->ecalEnergy();
            energy_hcal = energy_hcal + PFJetPart[j]->hcalEnergy();


	  }
	  	
          ptJet_      = i_jet->pt()*scale;
          etaJet_     = i_jet->eta();
          phiJet_     = i_jet->phi();
          chfJet_     = i_jet->chargedHadronEnergyFraction();
          nhfJet_     = i_jet->neutralHadronEnergyFraction();
          cemfJet_    = i_jet->chargedEmEnergyFraction();
          nemfJet_    = i_jet->neutralEmEnergyFraction();
          cmultiJet_  = i_jet->chargedMultiplicity();
          nmultiJet_  = i_jet->neutralMultiplicity();
         
	
          nrjets_     = jets->size();
         
	  dataTree_->Fill();
          njet++;
             
          PFJetPart.clear();
         
       }  
  
    }
 
  */
  
  //  cout << "end jets" << endl;
      
}
//////////////////////////////////////////////////////////////////////////////////////////
PFDataTreeProducer::~PFDataTreeProducer() 
{
  delete file_;
  delete dataTree_;
}
//}
