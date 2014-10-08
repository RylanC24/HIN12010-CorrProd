#include <TStopwatch.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "setupDiHadronCorrelation.h"
#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/DiHadronCorrelationMultiAnalyzerFWLite.h"
#endif

void RunDiHadronCorrelationMultiAnalysisFWLite(
           int ffrom, int fto,
           double zvtxmin, double zvtxmax,
           int nmin, int nmax,
           int centmin, int centmax)
{
  setup(ffrom, fto, zvtxmin, zvtxmax, nmin, nmax, centmin, centmax);

  TStopwatch *watch = new TStopwatch();

  // Run the analysis 
  TList* files = GetFileList(gFilesetname.Data());
  if(!files) {
    cerr << "Could not retrieve fileset: " << gFilesetname << endl;
    return;
  }

  vector<string> fileNames = GetDataVector(files,gFfrom,gFto);
//  fileNames.push_back("dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/edwenger/trkAnaSkim/MB_C10-PromptReco_v8_pptrkana_skim/v4/trkAnaSkimAOD_114_1.root");
  fwlite::ChainEvent event(fileNames);  

  //Initialize the analyzer      
  DiHadronCorrelationMultiAnalyzerFWLite* corr = new DiHadronCorrelationMultiAnalyzerFWLite(event);
  if(gTag.Contains("Isotropic")) corr->SetEventClass(2);
  if(gTag.Contains("Hard")) corr->SetEventClass(1);
  if(gTag.Contains("Soft")) corr->SetEventClass(0);
  corr->SetCheckSign(gCheckSign);
  corr->SetCutParameters(gCut);
  corr->SetNEtaPhiBins(gNEtaBins,gNPhiBins);
  corr->SetBkgFactor(gBkgFactor);
  if(gEPIndex>=0) corr->SetEventPlaneIndex(gEPIndex);

  TString strtrgid = "Track";
  if(gSystem->Getenv("TRGID")) strtrgid = gSystem->Getenv("TRGID");
  DiHadronCorrelationMultiBaseFWLite::ParticleType gTrgID = corr->GetParticleID(strtrgid);
  TString strassid = "Track";
  if(gSystem->Getenv("ASSID")) strassid = gSystem->Getenv("ASSID");
  DiHadronCorrelationMultiBaseFWLite::ParticleType gAssID = corr->GetParticleID(strassid);
  corr->SetTrgID(gTrgID);
  corr->SetAssID(gAssID);

  if(gCentfilename.Contains("root"))
  {
    cout<<"Running correlation analysis in Heavy Ion!"<<endl;
    TFile* fcentfile = new TFile(gCentfilename.Data());
    if(fcentfile->IsOpen()) corr->SetCentrality(fcentfile,gCenttablename,gNCentBins,gCentRunnum);
    else cout<<"Centrality table cannot be opened!"<<endl;
  }
  else cout<<"Running correlation analysis in pp!"<<endl;

  if(gEffhistname.Contains("root"))
  {
    TFile* feffhist = new TFile(gEffhistname.Data());    
    if(feffhist->IsOpen())
    {
      TH3D* htotaleff = 0;
      if(!gCut.IsHI) htotaleff = (TH3D*)feffhist->Get("rTotalEff3D");
      else
      {
        int histcentbin = -1;
        if((gCut.centmin>=0 && gCut.centmax<=2) || gCut.centmin>=50) histcentbin = 0;
        if(gCut.centmin>=2 && gCut.centmax<=4) histcentbin = 1;
        if(gCut.centmin>=4 && gCut.centmax<=12) histcentbin = 2;
        if(gCut.centmin>=12 && gCut.centmax<=20) histcentbin = 3;
        if(gCut.centmin>=20 && gCut.centmax<=40) histcentbin = 4;
        if(gCut.centmin==-1 && gCut.centmax==-1) histcentbin = 4;

        htotaleff = (TH3D*)feffhist->Get(Form("Tot_%d",histcentbin));
      }
      corr->LoadEffWeight(htotaleff);
      cout<<"Tracking efficiency weighting histogram is loaded!"<<endl;
    }
    else cout<<"Tracking efficiency weighting file cannot be opened!"<<endl;
  }
  else cout<<"No Tracking efficiency weighting histogram is found, or it's running MC!"<<endl;

  if(gTrghistname.Contains("root"))
  {
    TFile* ftrghist = new TFile(gTrghistname.Data());
    if(ftrghist->IsOpen())
    {
      corr->LoadTrgWeight((TH1D*)ftrghist->Get("trgEff"));
      cout<<"Triggering efficiency weighting histogram is loaded!"<<endl;
    }
    else cout<<"Triggering efficiency weighting file cannot be opened!"<<endl;
  }
  else cout<<"No Triggering efficiency weighting histogram is found, or it's running MC!"<<endl;

  if(gPileupdistfunchistname.Contains("root"))
  {
    TFile* fpileupdistfunchist = new TFile(gPileupdistfunchistname.Data());
    if(fpileupdistfunchist->IsOpen())
    {
      corr->LoadPileUpDistFunc((TH1D*)fpileupdistfunchist->Get("distfunc"));
      cout<<"Pileup distfunc histogram is loaded!"<<endl;
    }
    else cout<<"Pileup distfunc file cannot be opened!"<<endl;
  }

  corr->Process();

  // open output file
//  TString outputfilename = Form("%s/%s/unmerged/%s_%s_ffrom%d_fto%d_vtxmin%.1f_vtxmax%.1f_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",gOutputDir.Data(),gEventtype.Data(),gEventtype.Data(),gTag.Data(),gFfrom,gFto,gCut.zvtxmin,gCut.zvtxmax,gCut.nmin,gCut.nmax,gCut.etatrgmin,gCut.etatrgmax,gCut.etaassmin,gCut.etaassmax,gCut.centmin,gCut.centmax);
  TString outputfilename = Form("/tmp/%s_%s_ffrom%d_fto%d_vtxmin%.1f_vtxmax%.1f_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",gEventtype.Data(),gTag.Data(),gFfrom,gFto,gCut.zvtxmin,gCut.zvtxmax,gCut.nmin,gCut.nmax,gCut.etatrgmin,gCut.etatrgmax,gCut.etaassmin,gCut.etaassmax,gCut.centmin,gCut.centmax);
  TFile* outf = new TFile(outputfilename.Data(),"recreate");

  // Save the outputs
  TList* output = corr->GetOutputs();
  if(outf && output) {
    outf->cd();
    output->Write();
    cout << "Output file " <<  outf->GetName()
    << " written: " << endl;
    outf->ls();
    outf->Close();
    if(outf) delete outf;
  }  
  delete output;

  watch->Stop();
  watch->Print();
}
