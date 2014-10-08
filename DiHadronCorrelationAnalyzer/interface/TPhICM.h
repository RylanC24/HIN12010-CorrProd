#ifndef TPhICM_
#define TPhICM_

#ifndef ROOT_Rtypes
#include "Rtypes.h"
#endif

#include <fstream>
#include <vector>

//////////////////////////////////////////////////////////////////////////
//
// $Id:
// Author: B. Alver alver@mit.edu
// 
// Isotropic Cluster Model Class
//
//////////////////////////////////////////////////////////////////////////

class TF1;
class TList;
class TString;
class TFile;
class TTree;
class TH1F;
class TH1D;
class TH2D;
class TH3D;
class TClonesArray;

typedef struct { 
   Double_t x;
   Double_t y;
   Double_t z;
   Double_t px; 
   Double_t py; 
   Double_t pz; 
   Double_t m; 
   Int_t n;
   Int_t type;
} Cluster_t;

class TPhICM {
 protected:
   Int_t fNevents;
   Int_t fMultiplicity;
   Int_t fInitMultiplicity;

   Double_t fNuclR;
   Double_t fB;
   Double_t fEccA;
   Double_t fEccB;
   Double_t fAbsPar;
   Double_t fBetaRadial; 

   Double_t fClusProb[6];
   Int_t fClusNchild[6];
   Double_t fClusM[6];
   Double_t fClusEtaBroad[6];
   TF1* fClusMf[6];

   TF1 *fClusPtDist[6];
   TF1 *fClusPzDist[6];

   Bool_t fErrorflag;
   Bool_t fConserveMom;

   std::vector<Double_t> *fPhiVec;
   std::vector<Double_t> *fEtaVec;
   std::vector<Double_t> *fPtVec;
   std::vector<Bool_t> *fChargeVec;
   TTree  *fOutTree;
   TFile  *fOutTreef;

   int EventKind;
   TClonesArray* MonteCarlo;

   TF1 *fEtaShift;
   Double_t fEtaBroaden;

   void RunEvent(Int_t evtno);
   void WriteTreeEvent();
  
   Double_t ProbPathLength(Double_t PathLength); 
   Int_t    GenClusType();
   Int_t    GenClusNChild(Int_t clustype);
   Double_t  GenClusMass(Int_t clustype);
   void     GenInitEcc();
   void     GenClusPos(Double_t& clusX, Double_t& clusY);

   void FixModelParameters(Int_t modelnum);
   void ConserveMomentumInEtaRange(std::vector<Cluster_t> &clusterList,
				   Double_t etamin=-6,Double_t etamax=6);
 public:
   TPhICM(Int_t nevents, Int_t seed, Int_t multiplicity=100, Int_t modelnum=0);
   void    Run();
   void    SetInitParameters(Double_t r=13, Double_t b=0, Double_t absp=3) {fNuclR=r; fB=b; fAbsPar=absp;}
   void    SetBetaRadialFlow(Double_t beta=0.0) {fBetaRadial=beta;}
   void    SetEtaBroaden(Double_t etabroaden=0) {fEtaBroaden=etabroaden;}

   void    SetWriteTree(TString outfname);

   virtual ~TPhICM() {}
};

#endif
