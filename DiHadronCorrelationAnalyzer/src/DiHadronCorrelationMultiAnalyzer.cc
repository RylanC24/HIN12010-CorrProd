#include "../interface/DiHadronCorrelationMultiAnalyzer.h"
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
#include <TFile.h>
#include <TList.h>
#include <TIterator.h>
#include <TLorentzVector.h>
#include <TClonesArray.h>
#include <TString.h>
#include <TVector3.h>
#include <vector> 
#include <iostream>
#include "Math/Vector3D.h"

using namespace std;

DiHadronCorrelationMultiAnalyzer::DiHadronCorrelationMultiAnalyzer(const edm::ParameterSet& iConfig) :
  DiHadronCorrelationMultiBase(iConfig),
  signalTrgEffWeight(0),
  bkgTrgEffWeight(0),
  bkgAssEffWeight(0)
{
  cutPara.IsSymmetrize=1;
  cutPara.IsHarmonics = iConfig.getParameter<bool>("IsHarmonics");
  cutPara.IsHarmonicsEta1Eta2=0;
  bkgFactor = 10; 
}

DiHadronCorrelationMultiAnalyzer::~DiHadronCorrelationMultiAnalyzer()
{}

void DiHadronCorrelationMultiAnalyzer::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  double etabinwidth = (cutPara.etatrgmax-cutPara.etaassmin-cutPara.etatrgmin+cutPara.etaassmax)/NEtaBins;
  double phibinwidth = 2*PI/NPhiBins;

  hDeltaZvtx = theOutputs->make<TH1D>("deltazvtx",";#Delta z_{vtx}",200,-1.0,1.0);

  for(int itrg=0;itrg<(int)(cutPara.pttrgmin.size());itrg++)
  {
    for(int jass=0;jass<(int)(cutPara.ptassmin.size());jass++)
    {
      hSignal[itrg][jass] = theOutputs->make<TH2D>(Form("signal_trg%d_ass%d",itrg,jass),";#Delta#eta;#Delta#phi",
                                     NEtaBins+1,cutPara.etatrgmin-cutPara.etaassmax-etabinwidth/2.,cutPara.etatrgmax-cutPara.etaassmin+etabinwidth/2.,                    
                                 NPhiBins-1,-(PI-phibinwidth)/2.0,(PI*3.0-phibinwidth)/2.0);

      hBackground[itrg][jass] = theOutputs->make<TH2D>(Form("background_trg%d_ass%d",itrg,jass),";#Delta#eta;#Delta#phi",
                                     NEtaBins+1,cutPara.etatrgmin-cutPara.etaassmax-etabinwidth/2.,cutPara.etatrgmax-cutPara.etaassmin+etabinwidth/2.,
                                     NPhiBins-1,-(PI-phibinwidth)/2.0,(PI*3.0-phibinwidth)/2.0);

      hCorrelation[itrg][jass] = theOutputs->make<TH2D>(Form("correlation_trg%d_ass%d",itrg,jass),";#Delta#eta;#Delta#phi",
                                     NEtaBins+1,cutPara.etatrgmin-cutPara.etaassmax-etabinwidth/2.,cutPara.etatrgmax-cutPara.etaassmin+etabinwidth/2.,
                                     NPhiBins-1,-(PI-phibinwidth)/2.0,(PI*3.0-phibinwidth)/2.0);

      hSignal_pt1pt2 = theOutputs->make<TH2D>("signal_pt1pt2",";p_{T,1};p_{T,2}", 50, 0, 5.0, 50, 0, 5.0);
      hBackground_pt1pt2 = theOutputs->make<TH2D>("background_pt1pt2",";p_{T,1};p_{T,2}", 50, 0, 5.0, 50, 0, 5.0);
      hCorrelation_pt1pt2 = theOutputs->make<TH2D>("correlation_pt1pt2",";p_{T,1};p_{T,2}", 50, 0, 5.0, 50, 0, 5.0);

      hSignal_eta1eta2[itrg][jass] = theOutputs->make<TH2D>(Form("signal_eta1eta2_trg%d_ass%d",itrg,jass),";#eta_{1};#eta_{2}",
                                     NEtaBins*2+1,cutPara.etatrgmin-etabinwidth/2.,cutPara.etatrgmax+etabinwidth/2.,
                                     NEtaBins*2+1,cutPara.etaassmin-etabinwidth/2.,cutPara.etaassmax+etabinwidth/2.);
      hBackground_eta1eta2[itrg][jass] = theOutputs->make<TH2D>(Form("background_eta1eta2_trg%d_ass%d",itrg,jass),";#eta_{1};#eta_{2}",
                                     NEtaBins*2+1,cutPara.etatrgmin-etabinwidth/2.,cutPara.etatrgmax+etabinwidth/2.,
                                     NEtaBins*2+1,cutPara.etaassmin-etabinwidth/2.,cutPara.etaassmax+etabinwidth/2.);
      hCorrelation_eta1eta2[itrg][jass] = theOutputs->make<TH2D>(Form("correlation_eta1eta2_trg%d_ass%d",itrg,jass),";#eta_{1};#eta_{2}",
                                     NEtaBins*2+1,cutPara.etatrgmin-etabinwidth/2.,cutPara.etatrgmax+etabinwidth/2.,
                                     NEtaBins*2+1,cutPara.etaassmin-etabinwidth/2.,cutPara.etaassmax+etabinwidth/2.);

      hSignal_phi1phi2[itrg][jass] = theOutputs->make<TH2D>(Form("signal_phi1phi2_trg%d_ass%d",itrg,jass),";#phi_{1};#phi_{2}",
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.,
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.);
      hBackground_phi1phi2[itrg][jass] = theOutputs->make<TH2D>(Form("background_phi1phi2_trg%d_ass%d",itrg,jass),";#phi_{1};#phi_{2}",
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.,
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.);
      hCorrelation_phi1phi2[itrg][jass] = theOutputs->make<TH2D>(Form("correlation_phi1phi2_trg%d_ass%d",itrg,jass),";#phi_{1};#phi_{2}",
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.,
                                     NPhiBins*2+1,-PI-phibinwidth/2.,PI+phibinwidth/2.);
    }
  }
  DiHadronCorrelationMultiBase::beginRun(iRun, iSetup);
}

void DiHadronCorrelationMultiAnalyzer::endRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{
  DiHadronCorrelationMultiBase::endRun(iRun, iSetup);

  if(!cutPara.IsCorr) return;
  
  cout<< "Start sorting the events!" << endl;
  std::sort(eventcorrArray.begin(),eventcorrArray.end());
  cout<< "Finish sorting the events!" << endl;
  
  cout<< "Start running correlation analysis!" << endl;
  for(unsigned int i=0;i<eventcorrArray.size();i++)
  {
    if( i % 100 == 0 ) cout << "Processing " << i << "th event" << endl;
    FillHistsSignal(eventcorrArray[i]);

    unsigned int mixstart = i+1;
    unsigned int mixend = i+1+bkgFactor;

    if(mixend>eventcorrArray.size()) mixend=eventcorrArray.size();
    for(unsigned int j=mixstart;j<mixend;j++)
    {
      if(eventcorrArray[i].centbin != eventcorrArray[j].centbin) break;
      double deltazvtx = eventcorrArray[i].zvtx-eventcorrArray[j].zvtx;
      hDeltaZvtx->Fill(deltazvtx);

      FillHistsBackground(eventcorrArray[i],eventcorrArray[j]);
    }
  }
  cout<< "Finish running correlation analysis!" << endl;

  NormalizeHists();
  cout<< "Finish normalizing the histograms!" << endl;
}

void DiHadronCorrelationMultiAnalyzer::NormalizeHists()
{
  for(int itrg=0;itrg<(int)(cutPara.pttrgmin.size());itrg++)
  {
    for(int jass=0;jass<(int)(cutPara.ptassmin.size());jass++)
    {
      if(hSignal[itrg][jass]->Integral()==0) continue;
      if(hBackground[itrg][jass]->Integral()==0) continue;

      double  etabinwidth = hSignal[itrg][jass]->GetXaxis()->GetBinWidth(1);
      double  phibinwidth = hSignal[itrg][jass]->GetYaxis()->GetBinWidth(1);
 
      hSignal[itrg][jass]->Scale(1.0/etabinwidth/phibinwidth);
      hSignal_eta1eta2[itrg][jass]->Scale(1.0/etabinwidth/etabinwidth);
      hSignal_phi1phi2[itrg][jass]->Scale(1.0/phibinwidth/phibinwidth);

      hBackground[itrg][jass]->Scale(1.0/etabinwidth/phibinwidth);
      hBackground_eta1eta2[itrg][jass]->Scale(1.0/etabinwidth/etabinwidth);
      hBackground_phi1phi2[itrg][jass]->Scale(1.0/phibinwidth/phibinwidth);

      hCorrelation[itrg][jass]->Add(hSignal[itrg][jass]);
      hCorrelation[itrg][jass]->Divide(hBackground[itrg][jass]);
      hCorrelation[itrg][jass]->Scale(hBackground[itrg][jass]->GetBinContent(hBackground[itrg][jass]->FindBin(0,0)));  

     if(hBackground_eta1eta2[itrg][jass]) hCorrelation_eta1eta2[itrg][jass]->Divide(hBackground_eta1eta2[itrg][jass]);
     if(hBackground_phi1phi2[itrg][jass]) hCorrelation_phi1phi2[itrg][jass]->Divide(hBackground_phi1phi2[itrg][jass]);
    }
  }
  if(hBackground_pt1pt2)
  {
    hCorrelation_pt1pt2->Add(hSignal_pt1pt2);
    hCorrelation_pt1pt2->Divide(hBackground_pt1pt2);
  }
}

//--------------- Calculate signal distributions ----------------------
void DiHadronCorrelationMultiAnalyzer::FillHistsSignal(const DiHadronCorrelationEvent& eventcorr)
{
  for(unsigned int itrg=0;itrg<cutPara.pttrgmin.size();itrg++)
    for(unsigned int jass=0;jass<cutPara.ptassmin.size();jass++)
    {
      unsigned int ntrgsize = eventcorr.pVect_trg[itrg].size();
      unsigned int nasssize = eventcorr.pVect_ass[jass].size();
      double nMultCorr_trg = eventcorr.nMultCorrVect_trg[itrg];
//      double nMultCorr_ass = eventcorr.nMultCorrVect_ass[jass];

      for(unsigned int ntrg=0;ntrg<ntrgsize;ntrg++)
      {
        TLorentzVector pvector_trg = (eventcorr.pVect_trg[itrg])[ntrg];   
        double effweight_trg = (eventcorr.effVect_trg[itrg])[ntrg];
        double chg_trg = (eventcorr.chgVect_trg[itrg])[ntrg];
        double eta_trg = pvector_trg.Eta();
        double phi_trg = pvector_trg.Phi();
        double pt_trg = pvector_trg.Pt();

        for(unsigned int nass=0;nass<nasssize;nass++)
        {
          TLorentzVector pvector_ass = (eventcorr.pVect_ass[jass])[nass];   
          double effweight_ass = (eventcorr.effVect_ass[jass])[nass];
          double chg_ass = (eventcorr.chgVect_ass[jass])[nass];
          double eta_ass = pvector_ass.Eta();
          double phi_ass = pvector_ass.Phi();
          double pt_ass = pvector_ass.Pt();

          // check charge sign
          if( (checksign == 0) && (chg_trg != chg_ass)) continue;
          if( (checksign == 1) && (chg_trg == chg_ass)) continue;

          double deltaPhi=GetDeltaPhi(phi_trg,phi_ass);
          double deltaEta=GetDeltaEta(eta_trg,eta_ass);

          if(deltaEta==0.0 && deltaPhi==0.0 && pt_trg==pt_ass) continue; // two particles are identical
//          if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue; // two particles are close 
//          if(fabs(deltaEta)<0.05 && fabs(deltaPhi)<0.05) continue; // two particles are close 

          // total weight
          double effweight = effweight_trg * effweight_ass;
//          if(cutPara.IsPtWeightAss) effweight = effweight / (pt_ass-ptMean2_ass[jass]/ptMean_ass[jass]) / (pt_trg-ptMean2_trg[itrg]/ptMean_trg[itrg]);
//          if(cutPara.IsPtWeightAss) effweight = effweight / pt_ass;
//          if(cutPara.IsPtWeightTrg) effweight = effweight / pt_trg;

//          nMultCorr_trg = 1; // turn off normalization temperorily
          // Fill dihadron correlation functions
          if(!cutPara.IsSymmetrize)
          {
            hSignal[itrg][jass]->Fill(deltaEta,deltaPhi,1.0/effweight/nMultCorr_trg);
            hSignal_eta1eta2[itrg][jass]->Fill(eta_trg,eta_ass,1.0/effweight/nMultCorr_trg);
            hSignal_phi1phi2[itrg][jass]->Fill(phi_trg,phi_ass,1.0/effweight/nMultCorr_trg);
          }
          else
          {
            hSignal[itrg][jass]->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal[itrg][jass]->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal[itrg][jass]->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal[itrg][jass]->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal[itrg][jass]->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal[itrg][jass]->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/effweight/nMultCorr_trg);
            hSignal_eta1eta2[itrg][jass]->Fill(eta_trg,eta_ass,1.0/effweight/nMultCorr_trg);
            hSignal_phi1phi2[itrg][jass]->Fill(phi_trg,phi_ass,1.0/effweight/nMultCorr_trg);
            hSignal_eta1eta2[itrg][jass]->Fill(eta_ass,eta_trg,1.0/effweight/nMultCorr_trg);
            hSignal_phi1phi2[itrg][jass]->Fill(phi_ass,phi_trg,1.0/effweight/nMultCorr_trg);
          }

          if(fabs(deltaPhi)<PI/8. && fabs(deltaEta)>2) hSignal_pt1pt2->Fill(pt_trg,pt_ass,1.0/effweight/nMultCorr_trg);
        }
      }
    } 
}

void DiHadronCorrelationMultiAnalyzer::FillHistsBackground(const DiHadronCorrelationEvent& eventcorr_trg, const DiHadronCorrelationEvent& eventcorr_ass)
{
  for(unsigned int itrg=0;itrg<cutPara.pttrgmin.size();itrg++)
    for(unsigned int jass=0;jass<cutPara.ptassmin.size();jass++)
    {
      unsigned int ntrgsize = eventcorr_trg.pVect_trg[itrg].size();
      unsigned int nasssize = eventcorr_ass.pVect_ass[jass].size();
      double nMultCorr_trg = eventcorr_trg.nMultCorrVect_trg[itrg];
//      double nMultCorr_ass = eventcorr_ass.nMultCorrVect_ass[jass];

      for(unsigned int ntrg=0;ntrg<ntrgsize;ntrg++)
      {
        TLorentzVector pvector_trg = (eventcorr_trg.pVect_trg[itrg])[ntrg];   
        double effweight_trg = (eventcorr_trg.effVect_trg[itrg])[ntrg];
        double chg_trg = (eventcorr_trg.chgVect_trg[itrg])[ntrg];
        double eta_trg = pvector_trg.Eta();
        double phi_trg = pvector_trg.Phi();
        double pt_trg = pvector_trg.Pt();

        for(unsigned int nass=0;nass<nasssize;nass++)
        {
          TLorentzVector pvector_ass = (eventcorr_ass.pVect_ass[jass])[nass];   
          double effweight_ass = (eventcorr_ass.effVect_ass[jass])[nass];
          double chg_ass = (eventcorr_ass.chgVect_ass[jass])[nass];
          double eta_ass = pvector_ass.Eta();
          double phi_ass = pvector_ass.Phi();
          double pt_ass = pvector_ass.Pt();

          // check charge sign
          if( (checksign == 0) && (chg_trg != chg_ass)) continue;
          if( (checksign == 1) && (chg_trg == chg_ass)) continue;

          double deltaPhi=GetDeltaPhi(phi_trg,phi_ass);
          double deltaEta=GetDeltaEta(eta_trg,eta_ass);

          if(deltaEta==0.0 && deltaPhi==0.0 && pt_trg==pt_ass) continue; // two particles are identical
//          if(fabs(deltaEta)<0.028 && fabs(deltaPhi)<0.02) continue; // two particles are close 

//          nMultCorr_trg = 1; // turn off normalization temperorily

          // total weight
          double effweight = effweight_trg * effweight_ass * nMultCorr_trg;
//          if(cutPara.IsPtWeightAss) effweight = effweight / (pt_ass-ptMean2_ass[jass]/ptMean_ass[jass]) / (pt_trg-ptMean2_trg[itrg]/ptMean_trg[itrg]);
//          if(cutPara.IsPtWeightAss) effweight = effweight / pt_ass;
//          if(cutPara.IsPtWeightTrg) effweight = effweight / pt_trg;

          // Fill dihadron correlation functions
          if(!cutPara.IsSymmetrize)
          {
            hBackground[itrg][jass]->Fill(deltaEta,deltaPhi,1.0/effweight);
            hBackground_eta1eta2[itrg][jass]->Fill(eta_trg,eta_ass,1.0/effweight);
            hBackground_phi1phi2[itrg][jass]->Fill(phi_trg,phi_ass,1.0/effweight);
          }
          else
          {
            hBackground[itrg][jass]->Fill(fabs(deltaEta),fabs(deltaPhi),1.0/4.0/effweight);
            hBackground[itrg][jass]->Fill(-fabs(deltaEta),fabs(deltaPhi),1.0/4.0/effweight);
            hBackground[itrg][jass]->Fill(fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/effweight);
            hBackground[itrg][jass]->Fill(-fabs(deltaEta),-fabs(deltaPhi),1.0/4.0/effweight);
            hBackground[itrg][jass]->Fill(fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/effweight);
            hBackground[itrg][jass]->Fill(-fabs(deltaEta),2*PI-fabs(deltaPhi),1.0/4.0/effweight);
            hBackground_eta1eta2[itrg][jass]->Fill(eta_trg,eta_ass,1.0/effweight);
            hBackground_phi1phi2[itrg][jass]->Fill(phi_trg,phi_ass,1.0/effweight);
            hBackground_eta1eta2[itrg][jass]->Fill(eta_ass,eta_trg,1.0/effweight);
            hBackground_phi1phi2[itrg][jass]->Fill(phi_ass,phi_trg,1.0/effweight);
          }

          if(fabs(deltaPhi)<PI/8. && fabs(deltaEta)>2) hBackground_pt1pt2->Fill(pt_trg,pt_ass,1.0/effweight);
        }
      }

    } 
}
