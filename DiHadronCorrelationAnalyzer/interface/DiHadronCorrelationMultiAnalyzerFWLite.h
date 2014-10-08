#ifndef DiHadronCorrelationMultiAnalyzerFWLite_
#define DiHadronCorrelationMultiAnalyzerFWLite_

#ifndef DiHadronCorrelationMultiBaseFWLite_
#include "DiHadronCorrelationMultiBaseFWLite.h"
#endif 

#include <vector>
#include <string>

class TH1D;
class TH2D;
class TH3D;
class TTree;
class TList;

class DiHadronCorrelationMultiAnalyzerFWLite : public DiHadronCorrelationMultiBaseFWLite {
   
 protected:

   // histograms
   TH1D*  hDeltaZvtx;
   TH2D* hSignal[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackground[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hCorrelation[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hSignalCosn[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackgroundCosn[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hSignalSinn[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackgroundSinn[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hSignalRhon[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackgroundRhon[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal2PEPCorrelator[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground2PEPCorrelator[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hSignal_eta1eta2[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackground_eta1eta2[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hCorrelation_eta1eta2[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hSignal_phi1phi2[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hBackground_phi1phi2[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hCorrelation_phi1phi2[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_eta[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_eta[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_eta[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_phi[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_phi[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_phi[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_phi_short[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_phi_short[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_phi_short[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_phi_long[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_phi_long[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_phi_long[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_eta_near[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_eta_near[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_eta_near[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignal_eta_away[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackground_eta_away[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hCorrelation_eta_away[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hMultSignalPair[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hMultBackgroundPair[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hSignalPhiTotal[MAXPTTRGBINS][MAXPTASSBINS];
   TH1D* hBackgroundPhiTotal[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hMultSignalPairVsPhiTotal[MAXPTTRGBINS][MAXPTASSBINS];
   TH2D* hMultBackgroundPairVsPhiTotal[MAXPTTRGBINS][MAXPTASSBINS];

   // parameters
   double signalTrgEffWeight;
   double bkgTrgEffWeight;
   double bkgAssEffWeight;
   int    bkgFactor;

   virtual void Process(); 
   virtual void MakeHists();
   virtual void NormalizeHists();
   virtual void DeleteHists();
   virtual void FillHistsSignal(const DiHadronCorrelationEvent& eventcorr);
   virtual void FillHistsBackground(const DiHadronCorrelationEvent& eventcorr_trg, const DiHadronCorrelationEvent& eventcorr_ass);
   virtual TList* GetOutputs();
   
 public:
   DiHadronCorrelationMultiAnalyzerFWLite(fwlite::ChainEvent&);
   virtual ~DiHadronCorrelationMultiAnalyzerFWLite();

   void SetBkgFactor(int bkgfactor) {bkgFactor = bkgfactor;}
};
#endif  // DiHadronCorrelationMultiAnalyzerFWLite_
