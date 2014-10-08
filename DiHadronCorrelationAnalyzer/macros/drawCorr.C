#include "/net/hisrv0001/home/davidlw/useful_macros/CF.C"
#include "/net/hisrv0001/home/davidlw/useful_macros/GetMultiJetCorrFunc.C"

void drawCorr_2D(TString eventtype1="PPData_Minbias_7TeV", TString tag1="TRIGGER85_All_INCLEFF1TRG1v6", int nmin=110, int nmax=1000, int centmin=-1, int centmax=-1, int ptmin=2, int ptmax=1, double zmin=0.95, double zmax=1.05, double etatrgmin=-2.4, double etatrgmax=2.4, double etaassmin=-2.4, double etaassmax=2.4)
{
  TString filename = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",eventtype1.Data(),eventtype1.Data(),tag1.Data(),nmin,nmax,etatrgmin,etatrgmax,etaassmin,etaassmax,centmin,centmax);
  TH2D* h = (TH2D*)Get2DCFv3(filename.Data(),ptmin,ptmax);

  h->SetAxisRange(h->GetMinimum()*zmin,h->GetMaximum()*zmax,"Z");

  TCanvas* c = new TCanvas("c","",550,500);
  c->SetBottomMargin(0.08);
  c->SetTheta(60.839);
  c->SetPhi(38.0172);
  h->SetAxisRange(-4.5,4.5,"X");
  h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*1.4);
  h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*1.4);
  h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*0.85);
  h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*0.85);
  h->Draw("surf1 fb");

  TLatex* latex = new TLatex();
  latex->SetNDC();
  latex->SetTextSize(0.8*latex->GetTextSize());
//  latex->DrawLatex(0.045,0.95,Form("CMS pPb  #sqrt{s_{NN}} = 5.02 TeV, %d #leq N < %d",nmin,nmax));
  latex->DrawLatex(0.045,0.95,Form("CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV, %d #leq N < %d",nmin,nmax));
  latex->DrawLatex(0.045,0.86,"1 < p_{T}^{trig} < 3 GeV/c");
  latex->DrawLatex(0.045,0.8,"1 < p_{T}^{assoc} < 3 GeV/c");

//  SaveCanvas(c,"pPb/corr",Form("corr2D_pPbNew_pt%d-%d_nmin%d_nmax%d",ptmin,ptmax,nmin,nmax));
  SaveCanvas(c,"pPb/corr",Form("corr2D_PbPb_pt%d-%d_nmin%d_nmax%d",ptmin,ptmax,nmin,nmax));

  TH2D* h_signal = (TH2D*)Get2DSignal(filename.Data(),ptmin,ptmax);
  TH2D* h_background = (TH2D*)Get2DBackground(filename.Data(),ptmin,ptmax);
  
  TCanvas* c_signal = new TCanvas("c_signal","c_signal",550,500);
  h_signal->Draw("surf1");
  TCanvas* c_background = new TCanvas("c_background","c_background",550,500);
  h_background->Draw("surf1");

  TFile* fout = new TFile(Form("results/results_2D_N%d%d_eff.root",nmin,nmax),"recreate");
//  TFile* fout = new TFile(Form("results/results_2D_HF170_noeff.root",nmin,nmax),"recreate");
  h->Write();
  h_signal->Write();
  h_background->Write();
  fout->Close();
}

void drawCorr_dphi(TString eventtype1="PPData_Minbias_7TeV", TString tag1="TRIGGER85_All_INCLEFF1TRG1v6", int nmin=110, int nmax=1000, int centmin=-1, int centmax=-1, int ptmin=2, int ptmax=1, double detamin=2.0, double detamax=4.0, double ymin=0.95, double ymax=1.05, double etatrgmin=-2.4, double etatrgmax=2.4, double etaassmin=-2.4, double etaassmax=2.4)
{ 
  TString filename = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",eventtype1.Data(),eventtype1.Data(),tag1.Data(),nmin,nmax,etatrgmin,etatrgmax,etaassmin,etaassmax,centmin,centmax);
  TH1D* h = (TH1D*)Get1DCFdPhiv3(filename.Data(),ptmin,ptmax,detamin,detamax);

/*  
  h->GetYaxis()->SetTitleSize(h->GetYaxis()->GetTitleSize()*0.6);
  h->GetYaxis()->SetLabelSize(h->GetYaxis()->GetLabelSize()*0.6);
  h->GetXaxis()->SetTitleSize(h->GetXaxis()->GetTitleSize()*0.6);
  h->GetXaxis()->SetLabelSize(h->GetXaxis()->GetLabelSize()*0.6);
  h->GetYaxis()->SetTitleOffset(h->GetYaxis()->GetTitleOffset()*1.4);
  h->GetXaxis()->SetTitleOffset(h->GetXaxis()->GetTitleOffset()*1.4);
  h->SetMarkerColor(kRed);
*/
  h->SetAxisRange(h->GetMinimum()*ymin,h->GetMaximum()*ymax,"Y");

  TCanvas* c = new TCanvas("c","",550,500);
  h->Draw("PE");
}

void compareCorr_dphi(double ymin=0.95, double ymax=1.05)
{
//  TString eventtype1="PythiaReco_Minbias_7TeV";
  TString eventtype1="PPData_Minbias_7TeV";
//  TString tag1="TuneZ2_MinBias_42x_INCLHFPHFPPtWeight00Inf";
  TString tag1="TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032";
//  TString tag1="MINBIAS_All_Apr21ReReco_INCLHFPHFMPtWeight035";
//  TString tag1="MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight0310";
  int nmin1=150;
  int nmax1=1000;
  int centmin1=-1;
  int centmax1=-1;
  double etatrgmin1=-5.0;
  double etatrgmax1=-3.0;
  double etaassmin1=3.0;
  double etaassmax1=5.0;

//  TString eventtype2="PythiaReco_Minbias_7TeV";
  TString eventtype2="PPData_Minbias_7TeV";
//  TString tag2="TuneZ2_MinBias_42x_INCLHFPHFPPtWeight00Inf";
//  TString tag2="MINBIAS_All_Apr21ReReco_INCLHFPHFMPtWeight03Inf";
  TString tag2="TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032";
//  TString tag2="MINBIAS_All_Apr21ReReco_INCLHFPHFM033";
  int nmin2=150;
  int nmax2=1000;
  int centmin2=-1; 
  int centmax2=-1;
  double etatrgmin2=-5.0;
  double etatrgmax2=-3.0;
  double etaassmin2=3.0;
  double etaassmax2=5.0;

  TString filename1 = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",eventtype1.Data(),eventtype1.Data(),tag1.Data(),nmin1,nmax1,etatrgmin1,etatrgmax1,etaassmin1,etaassmax1,centmin1,centmax1);
  TString filename2 = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",eventtype2.Data(),eventtype2.Data(),tag2.Data(),nmin2,nmax2,etatrgmin2,etatrgmax2,etaassmin2,etaassmax2,centmin2,centmax2);

  int ptmin1=0;
  int ptmax1=0;
  double detamin1=-15;
  double detamax1=15;

  int ptmin2=0;
  int ptmax2=0;
  double detamin2=-15;
  double detamax2=15;

  TH1D* h1 = (TH1D*)Get1DCFdPhiv3(filename1.Data(),ptmin1,ptmax1,detamin1,detamax1);
  TH1D* h1_zyam=GetJetCorrFunc1D_ZYAM_pp(h1);
  h1_zyam->SetAxisRange(h1_zyam->GetMinimum()*ymin,h1_zyam->GetMaximum()*ymax,"Y");

  TH1D* h2 = (TH1D*)Get1DCFdPhiv3(filename2.Data(),ptmin2,ptmax2,detamin2,detamax2);
  TH1D* h2_zyam=GetJetCorrFunc1D_ZYAM_pp(h2);
  h2_zyam->SetAxisRange(h2_zyam->GetMinimum()*ymin,h2_zyam->GetMaximum()*ymax,"Y");

//  TF1* fit1 = new TF1("fit1","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x)+2*[5]*cos(5*x)+2*[6]*cos(6*x)+2*[7]*cos(7*x)+2*[8]*cos(8*x)+2*[9]*cos(9*x)+2*[10]*cos(10*x)+2*[11]*cos(11*x)+2*[12]*cos(12*x))",-PI/2,3*PI/2);
//  TF1* fit2 = new TF1("fit2","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x)+2*[4]*cos(4*x)+2*[5]*cos(5*x)+2*[6]*cos(6*x)+2*[7]*cos(7*x)+2*[8]*cos(8*x)+2*[9]*cos(9*x)+2*[10]*cos(10*x)+2*[11]*cos(11*x)+2*[12]*cos(12*x))",-PI/2,3*PI/2);
  TF1* fit1 = new TF1("fit1","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x))",-PI/2,3*PI/2);
  TF1* fit2 = new TF1("fit2","[0]*(1+2*[1]*cos(x)+2*[2]*cos(2*x)+2*[3]*cos(3*x))",-PI/2,3*PI/2);
  

  h1->Fit("fit1","RNO");
  h2->Fit("fit2","RNO");

//cout<<"v2="<<sqrt(fit2->GetParameter(2)-3*fit1->GetParameter(2)*fit1->GetParameter(0)/fit2->GetParameter(0))<<endl;
//  cout<<"v2="<<sqrt(fit2->GetParameter(2))<<endl;

  TCanvas* c1 = new TCanvas("c1","",800,400);
  c1->Divide(2,1);
  c1->cd(1);
  h1->Draw("PE");
  fit1->Draw("Lsame");
  c1->cd(2);
  h2->Draw("PE");
  fit2->Draw("Lsame");

  TCanvas* c = new TCanvas("c","",550,500);
  h2_zyam->Draw("PE");
  h1_zyam->SetMarkerStyle(24);
  h1_zyam->Scale(1.);
  h1_zyam->Draw("PESAME");

}
