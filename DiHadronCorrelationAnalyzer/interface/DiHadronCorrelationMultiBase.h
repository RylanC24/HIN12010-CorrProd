#ifndef DiHadronCorrelationMultiBase_
#define DiHadronCorrelationMultiBase_

#include <exception>
#include <vector>
#include <string>
#include <memory>
#include <iostream>
#include <math.h>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/BasicClusterFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/HcalRecHit/interface/HBHERecHit.h"
#include "DataFormats/HcalRecHit/interface/HFRecHit.h"
#include "DataFormats/HcalRecHit/interface/HORecHit.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "RecoHI/HiCentralityAlgos/interface/CentralityProvider.h"
//#include "DataFormats/HeavyIonEvent/interface/CentralityProvider.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidate.h"
#include "DataFormats/Candidate/interface/VertexCompositeCandidateFwd.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "RecoJets/JetAlgorithms/interface/JetAlgoHelper.h"
#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/CutParameters.h"
#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/DiHadronCorrelationEvent.h"

class TBranch;
class TClonesArray;
class TH1D;
class TH2D;
class TH3D;
class TList;
class TTree;
class TNtuple;
class TVector3;
class TLorentzVector;
class DiHadronCorrelationEvent;

using namespace std;
using namespace reco;
using namespace edm;

#define PI 3.141592

class DiHadronCorrelationMultiBase : public edm::EDAnalyzer {
   
 protected:
   
   // particle types
   enum ParticleType {
     kGenParticle=0,
     kTrack=1,
     kTrackP=2,
     kTrackM=3,
     kJet=4,
     kMuon=5,
     kPhoton=6,
     kEcalSC=7,
     kEcalBC=8,
     kPFHadron=9,
     kPFPhoton=10,
     kPFHadronPhoton=11,
     kPFPhotonTrack=12,
     kPFElectron=13,
     kPFMuon=14,
     kPFNeutral=15,
     kPFHadronHF=16,
     kPFEgammaHF=17,
     kPionZero=18,
     kPionZeroBkgUp=19,
     kPionZeroBkgDown=20,
     kCaloTower=21,
     kEcalRecHit=22,
     kHcalRecHit=23,
     kKshort=24,
     kLambda=25
   };

   ParticleType  trgID;
   ParticleType  assID;

   //  vectors for tracks, hits etc
   DiHadronCorrelationEvent* eventcorr;
   vector<DiHadronCorrelationEvent> eventcorrArray;

   CentralityProvider * cent;
   edm::Service<TFileService> theOutputs;

   // histograms
   TH2D*  hEffWeight;
   TH1D*  hTrgWeight;
   TH1D*  hNVtx;
   TH2D*  hNVtxVsNMult;
   TH1D*  hZVtx;
   TH2D*  hXYVtx;
   TH1D*  hMultRawAll;
   TH1D*  hMultCorrAll;
   TH1D*  hPtAll_trg;
   TH2D*  hdNdetadptAll_trg;
   TH2D*  hdNdetadphiAll_trg;
   TH1D*  hPtAll_ass;
   TH2D*  hdNdetadptAll_ass;
   TH2D*  hdNdetadphiAll_ass;
   TH1D*  hPtCorrAll_trg;
   TH2D*  hdNdetadptCorrAll_trg;
   TH2D*  hdNdetadphiCorrAll_trg;
   TH1D*  hPtCorrAll_ass;
   TH2D*  hdNdetadptCorrAll_ass;
   TH2D*  hdNdetadphiCorrAll_ass;
   TH2D*  hEventEngineer;
   TH1D*  hCentrality;
   TH1D*  hHFTowerSum; 
   TH2D*  hHFvsNpixel;
   TH2D*  hHFvsZDC;
   TH1D*  hNpart;
   TH2D*  hV0InvMassVsPt;
   TH2D*  hV0InvMassVsPt_masspipicut;
   TH2D*  hV0AP;
   TH2D*  hV0AP_masspipicut;
   TH1D*  hV0MassPiPi;
   TH1D*  hV0MassEE;
   TH2D*  hInvMassVsPt_Signal;

   TH2D*  hdNdetadphi_trg[MAXPTTRGBINS];
   TH2D*  hdNdetadphi_ass[MAXPTASSBINS];
   TH2D*  hdNdetadphiCorr_trg[MAXPTTRGBINS];
   TH2D*  hdNdetadphiCorr_ass[MAXPTASSBINS];
   TH1D*  hMult_trg[MAXPTTRGBINS];
   TH1D*  hMult_ass[MAXPTASSBINS];
   TH1D*  hMultCorr_trg[MAXPTTRGBINS];
   TH1D*  hMultCorr_ass[MAXPTASSBINS];
   TH1D*  hpT_Signal_trg[MAXPTTRGBINS];
   TH1D*  hpT_Signal_ass[MAXPTASSBINS];
   TH1D*  hpTCorr_Signal_trg[MAXPTTRGBINS];
   TH1D*  hpTCorr_Signal_ass[MAXPTASSBINS];
   TH1D*  hpT_Background_trg[MAXPTTRGBINS];
   TH1D*  hpT_Background_ass[MAXPTASSBINS];
   TH1D*  hpTCorr_Background_trg[MAXPTTRGBINS];
   TH1D*  hpTCorr_Background_ass[MAXPTASSBINS];

   // parameters
   CutParameters cutPara;
   int  checksign;
   int  eventClass;
   int  nEvent; 
   int  nMult;
   int  nMultAll_trg;
   int  nMultAll_ass;
   double nMultCorr;
   double nMultAllCorr_trg;
   double nMultAllCorr_ass;
   int  nVertices;
   unsigned int  maxofflinetracks;
   int nCentBins;
   int  NEtaBins;
   int  NPhiBins;      
   double  xVtx;
   double  yVtx;
   double  zVtx;
   double  xVtxError;
   double  yVtxError;
   double  zVtxError;
   int  hiCentrality;
   unsigned int nMult_trg[MAXPTTRGBINS];
   unsigned int nMult_ass[MAXPTASSBINS];
   double nMultCorr_trg[MAXPTTRGBINS];
   double nMultCorr_ass[MAXPTASSBINS];
   double ptMean_trg[MAXPTTRGBINS];
   double ptMean_ass[MAXPTASSBINS];
   double ptMean2_trg[MAXPTTRGBINS];
   double ptMean2_ass[MAXPTASSBINS];
  
   virtual void analyze(const edm::Event&, const edm::EventSetup&);
   //virtual void endJob() ;
   virtual void beginRun(const edm::Run&, const edm::EventSetup&);
   virtual void endRun(const edm::Run&, const edm::EventSetup&);

   virtual void GetMult(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   virtual void LoopParticles(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool istrg, int pdgid, bool isstable, bool ischarge);    
   virtual void LoopTracks(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool istrg, int icharge=999);
   virtual void LoopCaloTower(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool istrg);
   virtual void LoopV0Candidates(const edm::Event& iEvent, const edm::EventSetup& iSetup, bool istrg, TString candtype);
   virtual void GetVertices(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   virtual void NormalizeHists() {};
   virtual void AssignTrgPtBins(double pt, double eta, double phi, double mass, double charge, double effweight);
   virtual void AssignAssPtBins(double pt, double eta, double phi, double mass, double charge, double effweight);
   virtual double GetEventEngineer(const edm::Event& iEvent, const edm::EventSetup& iSetup, int nn);
   virtual int  GetCentralityBin(const edm::Event& iEvent, const edm::EventSetup& iSetup);
   virtual double GetDeltaEta(double eta_trg, double eta_ass);
   virtual double GetDeltaPhi(double phi_trg, double phi_ass);
   virtual bool   GetPttrgBin(double pt, double eta, int itrg);
   virtual bool   GetPtassBin(double pt, double eta, int jass);
   virtual double GetTrgWeight(double nmult);
   virtual double GetEffWeight(double eta, double pt, double zvtx, int centbin);
   virtual ParticleType GetParticleID(TString particleid);
   
 public:

   explicit DiHadronCorrelationMultiBase(const edm::ParameterSet&);
   ~DiHadronCorrelationMultiBase() {};
};
#endif  // DiHadronCorrelationMultiBase_
