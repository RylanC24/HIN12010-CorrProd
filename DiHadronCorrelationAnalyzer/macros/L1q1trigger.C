#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TGraph.h>
#include <TMath.h>
#include <TString.h>
#include <stdio.h> 
#include <TLegendEntry.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>

#include <map>

using namespace std;
//#include "l1ExtraTree.h"
//#include "l1Tree.h"

long makeKey(long run, long event){
  return (10000000*run + event);
}
//long makeKey(long run, long lumi, long event){
//  return (10000000000*run + 10000000*lumi + event);
//}


void L1q1trigger(){
  const TString l1_input = "~/scratch1/DiHadronCorrelations/L1UpgradeAnalyzer.root";

  TFile *lFile = TFile::Open(l1_input);
  Int_t l1Up_evt, l1Up_run, l1Up_et, l1Up_q2;
  TTree *l1UpTree = (TTree*)lFile->Get("L1UpgradeAnalyzer/L1UpgradeTree");

  l1UpTree->SetBranchStatus("*",0);
  l1UpTree->SetBranchStatus("event",1);
  l1UpTree->SetBranchStatus("run",1);
  l1UpTree->SetBranchStatus("centrality_hwPt",1);
  l1UpTree->SetBranchStatus("v2_hwPt",1);

  l1UpTree->SetBranchAddress("event",&l1Up_evt);
  l1UpTree->SetBranchAddress("run",&l1Up_run);
  l1UpTree->SetBranchAddress("centrality_hwPt",&l1Up_et);
  l1UpTree->SetBranchAddress("v2_hwPt",&l1Up_q2);

  const TString forest_input = "~/scratch1/DiHadronCorrelations/outputs_312/HIData_Minbias_2760GeV/ebyeflow_promptreco_pt009.root";
  TFile *fFile = TFile::Open(forest_input);
  TTree *fTree = (TTree*)fFile->Get("ebyeflow_ana_HI/q2ntuple");

  Float_t f_evt, f_run, f_lumi, hiBin;
  float hiHF,q2hf;

  fTree->SetBranchStatus("*",0);
  fTree->SetBranchStatus("evtnum",1);
  fTree->SetBranchStatus("run",1);
  fTree->SetBranchStatus("lumi",1);
  fTree->SetBranchStatus("cent",1);
  fTree->SetBranchStatus("q2hf",1);
  fTree->SetBranchStatus("hf",1);

  fTree->SetBranchAddress("evtnum",&f_evt);
  fTree->SetBranchAddress("run",&f_run);
  fTree->SetBranchAddress("lumi",&f_lumi);
  fTree->SetBranchAddress("cent",&hiBin);
  fTree->SetBranchAddress("hf",&hiHF);
  fTree->SetBranchAddress("q2hf",&q2hf);

  map<long, int> kmap;
  map<long, int> kmapcal;
  
  cout << "Filling the map..." << endl;
  int l1up_entries = l1UpTree->GetEntries();
  //int l_entries = lEvtTree->GetEntries();
  for(long j = 0; j < l1up_entries; ++j){
    if(j % 10000 == 0) printf("%ld / %d\n",j,l1up_entries);
    l1UpTree->GetEntry(j);
    long key = makeKey(l1Up_run, l1Up_evt);

    pair<long,int> p(key,j);
    kmap.insert(p);
    kmapcal.insert(p);
  }
  cout << "map filled!!!" << endl;

//q2 histos
  TH2D* q2CorrHist = new TH2D("q2CorrHist",";offline q_{2};online q_{2}",1000,0.,1.0,1000,0.0,1.0);
  TH2D* HFCorrHist = new TH2D("HFCorrHist",";offline HF E_{T} sum (GeV); online HF E_{T} sum (GeV)",100,0.,8000.0,100,0.0,8000.0);
  TH2D* q2HFCorrHist = new TH2D("q2HFCorrHist",";online HF*q^{2}_{2} (GeV); online HF E_{T} sum (GeV)",300,0.,30000.0,100,0.0,8000.0);

  Float_t q2On   = 0.0; 

  int entries = fTree->GetEntries();
  for(long int j = 1; j < entries; ++j){
    if(j % 10000 == 0) printf("%ld / %d \n",j,entries);
    fTree->GetEntry(j);
    long keycal = makeKey(f_run, f_evt);
    
    map<long,int>::const_iterator gotcal = kmapcal.find(keycal);

    if(gotcal == kmapcal.end()){
      continue;      
    }
    else {
      l1UpTree->GetEntry(gotcal->second);
      kmapcal.erase(keycal);

      if(l1Up_et == 0) continue; 
      //q2 part
                
      q2On = sqrt(l1Up_q2)/l1Up_et;

//      if(l1Up_run<182060 && l1Up_run>182050 & l1Up_et<2839 && l1Up_et>1198)
    //  if(l1Up_et<2839 && l1Up_et>1198)
      {
        q2CorrHist->Fill(q2hf,q2On);
        HFCorrHist->Fill(hiHF,l1Up_et);
        q2HFCorrHist->Fill(l1Up_q2,l1Up_et);
      }
      q2On   = 0.0; 
    }  
  }  

  TCanvas* c = new TCanvas("c","c",800,400);
  c->Divide(2,1);
  c->cd(1);
  c->GetPad(1)->SetLogz();
  q2CorrHist->Rebin2D(2,2);
  q2CorrHist->Draw("colz");
  TProfile* hprof = q2CorrHist->ProfileX("prof",-1,-1,"e");
  hprof->Draw("PEsame");
  c->cd(2);
  c->GetPad(2)->SetLogz();
  HFCorrHist->Draw("colz");

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  q2HFCorrHist->Draw("colz");
}
