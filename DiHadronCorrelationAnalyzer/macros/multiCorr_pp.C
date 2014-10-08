#include <fstream>
#include <iomanip>

void makeMultiCorr()
{
  const int ntrgptbins=5;
  const int nassptbins=5;

  TString dataname2 = "HIData_Minbias_2760GeV";
  TString datatag2 = "PPRereco_INCLEFF1v3HydjetEff";
  int nmin2=220;
  int nmax2=260;

  TString dataname = "PAData_Minbias_5TeV";
  TString datatag = "HM_PromptRecoNew_All_INCLEFF1v3";
  int nmin=220;
  int nmax=260;
/*
  TString dataname3 = "PPData_Minbias_7TeV";
  TString datatag3 = "TRIGGER_All_Apr21ReReco_INCLEFF1v3";
  int nmin3=150;
  int nmax3=1000;
*/
  TString dataname3 = "PAData_Minbias_5TeV";
  TString datatag3 = "MB_PromptRecoNew_All_INCLEFF1v3";
  int nmin3=0;
  int nmax3=20;

  TString dataname4 = "HIData_Minbias_2760GeV";
  TString datatag4 = "PPRereco_INCLEFF1v3HydjetEff";
  int nmin4=0;
  int nmax4=20;

  TH1D ***hdata_long = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname.Data(),datatag.Data(),nmin,nmax,-1,-1,-1.0,1.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D ***hdata_short = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname.Data(),datatag.Data(),nmin,nmax,-1,-1,2.0,4.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D* hdata_long_zyam[ntrgptbins][nassptbins];
  TH1D* hdata_short_zyam[ntrgptbins][nassptbins];

  TH1D ***hdata2_long = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname2.Data(),datatag2.Data(),nmin2,nmax2,-1,-1,-1.0,1.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D ***hdata2_short = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname2.Data(),datatag2.Data(),nmin2,nmax2,-1,-1,2.0,4.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D* hdata2_long_zyam[ntrgptbins][nassptbins];
  TH1D* hdata2_short_zyam[ntrgptbins][nassptbins];

  TH1D ***hdata3_long = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname3.Data(),datatag3.Data(),nmin3,nmax3,-1,-1,-1.0,1.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D ***hdata3_short = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname3.Data(),datatag3.Data(),nmin3,nmax3,-1,-1,2.0,4.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D* hdata3_long_zyam[ntrgptbins][nassptbins];
  TH1D* hdata3_short_zyam[ntrgptbins][nassptbins];

  TH1D ***hdata4_long = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname4.Data(),datatag4.Data(),nmin4,nmax4,-1,-1,-1.0,1.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D ***hdata4_short = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname4.Data(),datatag4.Data(),nmin4,nmax4,-1,-1,2.0,4.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D* hdata4_long_zyam[ntrgptbins][nassptbins];
  TH1D* hdata4_short_zyam[ntrgptbins][nassptbins];

  TF1* fitfunc_vn[ntrgptbins][nassptbins];
  TF1* fitfunc2_vn[ntrgptbins][nassptbins];
  TF1* fitfunc_vn_zyam[ntrgptbins][nassptbins];
  TF1* fitfunc2_vn_zyam[ntrgptbins][nassptbins];
  TF1* zyamfunc_short[ntrgptbins][nassptbins];
  TF1* zyamfunc2_short[ntrgptbins][nassptbins];

  TCanvas* cc = new TCanvas("cc","cc",250*ntrgptbins,280*2);
  makeMultiPanelCanvas(cc,ntrgptbins,2,0.01,0.02,0.25,0.18,0.14);
  TCanvas* cc1 = new TCanvas("cc1","cc1",250*ntrgptbins,280*2);
  makeMultiPanelCanvas(cc1,ntrgptbins,2,0.01,0.02,0.25,0.18,0.14);
  TCanvas* cc2 = new TCanvas("cc2","cc2",250*ntrgptbins,280*2);
  makeMultiPanelCanvas(cc2,ntrgptbins,2,0.01,0.02,0.25,0.18,0.14);
  TCanvas* cc3 = new TCanvas("cc3","cc3",250*(ntrgptbins-1),310*2);
  makeMultiPanelCanvas(cc3,ntrgptbins-1,2,0.01,0.02,0.25,0.14,0.05);

    // be careful
    int jass=1;
    hdata4_short[4][1]->SetBinContent(hdata4_short[4][1]->FindBin(3.1416),hdata4_short[4][1]->GetBinContent(hdata4_short[4][1]->FindBin(3.1416))*0.7);
    hdata4_long[4][1]->SetBinContent(hdata4_long[4][1]->FindBin(0),hdata4_long[4][1]->GetBinContent(hdata4_long[4][1]->FindBin(0))*0.85);

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      cc->cd(itrg+1);
      hdata_short[itrg][jass]->SetMarkerColor(2);     
      hdata_short[itrg][jass]->Draw("PE");
      hdata3_short[itrg][jass]->SetMarkerColor(4);
      hdata3_short[itrg][jass]->SetMarkerStyle(24);
      hdata3_short[itrg][jass]->Draw("PESAME");
      hdata2_short[itrg][jass]->SetMarkerColor(4);
      hdata2_short[itrg][jass]->SetMarkerStyle(20);
      hdata2_short[itrg][jass]->Draw("PE");
      hdata4_short[itrg][jass]->SetMarkerColor(4);
      hdata4_short[itrg][jass]->SetMarkerStyle(20);
      hdata4_short[itrg][jass]->Draw("PESAME");

      cc->cd(ntrgptbins+itrg+1);
      hdata_long[itrg][jass]->SetMarkerColor(2);
      hdata_long[itrg][jass]->Draw("PE");
      hdata3_long[itrg][jass]->SetMarkerColor(4);
      hdata3_long[itrg][jass]->SetMarkerStyle(24);
      hdata3_long[itrg][jass]->Draw("PESAME");
      hdata2_long[itrg][jass]->SetMarkerColor(4);
      hdata2_long[itrg][jass]->SetMarkerStyle(21);
      hdata2_long[itrg][jass]->Draw("PE");
      hdata4_long[itrg][jass]->SetMarkerColor(4);
      hdata4_long[itrg][jass]->SetMarkerStyle(21);
      hdata4_long[itrg][jass]->Draw("PESAME");
   
      fitfunc_vn[itrg][jass] = FitVnFunc(hdata_short[itrg][jass]);
      fitfunc2_vn[itrg][jass] = FitVnFunc(hdata2_short[itrg][jass]);

      cc1->cd(itrg+1);
      zyamfunc2_short[itrg][jass]=GetFitFunc_ZYAM_pp(hdata2_short[itrg][jass]);
      hdata2_short_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata2_short[itrg][jass]);
      hdata2_short_zyam[itrg][jass]->SetAxisRange(-0.03,0.43,"Y");
      hdata2_short_zyam[itrg][jass]->SetMarkerStyle(21);
      hdata2_short_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata2_short_zyam[itrg][jass],1.9,3.3);
      hdata2_short_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata2_short_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.1);
      hdata2_short_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata2_short_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata2_short_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata2_short_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.2);
      hdata2_short_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata2_short_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.15);
      hdata2_short_zyam[itrg][jass]->Draw("PE");
      hdata4_short_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata4_short[itrg][jass]);
      hdata4_short_zyam[itrg][jass]->SetAxisRange(-0.03,0.43,"Y");
      hdata4_short_zyam[itrg][jass]->SetMarkerStyle(25);
      hdata4_short_zyam[itrg][jass]->SetMarkerColor(4);
      hdata4_short_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata4_short_zyam[itrg][jass],1.9,3.3);
      hdata4_short_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata4_short_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.1);
      hdata4_short_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata4_short_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata4_short_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata4_short_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.2);
      hdata4_short_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata4_short_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.15);
      hdata4_short_zyam[itrg][jass]->Draw("PESAME");
      fitfunc2_vn_zyam[itrg][jass] = new TF1(Form("fitfunc_vn_zyam_%d_%d",itrg,jass),"[0]+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x)+2*[3]*TMath::Cos(3*x)",-1.2,2*PI-1.2);
      fitfunc2_vn_zyam[itrg][jass]->SetParameters(fitfunc2_vn[itrg][jass]->GetParameter(0)-zyamfunc2_short[itrg][jass]->GetParameter(0),fitfunc2_vn[itrg][jass]->GetParameter(0)*fitfunc2_vn[itrg][jass]->GetParameter(1),fitfunc2_vn[itrg][jass]->GetParameter(0)*fitfunc2_vn[itrg][jass]->GetParameter(2),fitfunc2_vn[itrg][jass]->GetParameter(0)*fitfunc2_vn[itrg][jass]->GetParameter(3));
      fitfunc2_vn_zyam[itrg][jass]->SetLineStyle(1);
      fitfunc2_vn_zyam[itrg][jass]->SetLineWidth(1);
      fitfunc2_vn_zyam[itrg][jass]->Draw("LSAME");
      cc2->cd(itrg+1);
      zyamfunc_short[itrg][jass]=GetFitFunc_ZYAM_pp(hdata_short[itrg][jass]);
      hdata_short_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata_short[itrg][jass]);
      hdata_short_zyam[itrg][jass]->SetAxisRange(-0.03,0.43,"Y");
      hdata_short_zyam[itrg][jass]->SetMarkerStyle(20);
      hdata_short_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata_short_zyam[itrg][jass],1.9,3);
      hdata_short_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata_short_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata_short_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata_short_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata_short_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata_short_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata_short_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata_short_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata_short_zyam[itrg][jass]->Draw("PE");
      hdata3_short_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata3_short[itrg][jass]);
      hdata3_short_zyam[itrg][jass]->SetMarkerStyle(24);
      hdata3_short_zyam[itrg][jass]->SetMarkerColor(2);
      hdata3_short_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata3_short_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata3_short_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata3_short_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata3_short_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata3_short_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata3_short_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata3_short_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata3_short_zyam[itrg][jass]->Draw("PESAME");
      fitfunc_vn_zyam[itrg][jass] = new TF1(Form("fitfunc_vn_zyam_%d_%d",itrg,jass),"[0]+2*[1]*TMath::Cos(x)+2*[2]*TMath::Cos(2*x)+2*[3]*TMath::Cos(3*x)",-1.2,2*PI-1.2);
      fitfunc_vn_zyam[itrg][jass]->SetParameters(fitfunc_vn[itrg][jass]->GetParameter(0)-zyamfunc_short[itrg][jass]->GetParameter(0),fitfunc_vn[itrg][jass]->GetParameter(0)*fitfunc_vn[itrg][jass]->GetParameter(1),fitfunc_vn[itrg][jass]->GetParameter(0)*fitfunc_vn[itrg][jass]->GetParameter(2),fitfunc_vn[itrg][jass]->GetParameter(0)*fitfunc_vn[itrg][jass]->GetParameter(3));
      fitfunc_vn_zyam[itrg][jass]->SetLineStyle(1);
      fitfunc_vn_zyam[itrg][jass]->SetLineWidth(1);
      fitfunc_vn_zyam[itrg][jass]->Draw("LSAME");
      cc1->cd(ntrgptbins+itrg+1);
      hdata2_long_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata2_long[itrg][jass]);
      hdata2_long_zyam[itrg][jass]->SetAxisRange(-0.07,0.85,"Y");
      hdata2_long_zyam[itrg][jass]->SetMarkerStyle(21);
      hdata2_long_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata2_long_zyam[itrg][jass],1.9,3);
      hdata2_long_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata2_long_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata2_long_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata2_long_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata2_long_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata2_long_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata2_long_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata2_long_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata2_long_zyam[itrg][jass]->Draw("PE");
      hdata4_long_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata4_long[itrg][jass]);
      hdata4_long_zyam[itrg][jass]->SetAxisRange(-0.07,0.85,"Y");
      hdata4_long_zyam[itrg][jass]->SetMarkerStyle(25);
      hdata4_long_zyam[itrg][jass]->SetMarkerColor(4);
      hdata4_long_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata4_long_zyam[itrg][jass],1.9,3);
      hdata4_long_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata4_long_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata4_long_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata4_long_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata4_long_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata4_long_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata4_long_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata4_long_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata4_long_zyam[itrg][jass]->Draw("PESAME");
      cc2->cd(ntrgptbins+itrg+1);
      hdata_long_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata_long[itrg][jass]);
      hdata_long_zyam[itrg][jass]->SetAxisRange(-0.07,0.85,"Y");
      hdata_long_zyam[itrg][jass]->SetMarkerStyle(20);
      hdata_long_zyam[itrg][jass]->SetYTitle("#frac{1}{N_{trig}} #frac{dN^{pair}}{d#Delta#phi} - C_{ZYAM}");
      fixedFontHist1D(hdata_long_zyam[itrg][jass],1.9,3);
      hdata_long_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata_long_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata_long_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata_long_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata_long_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata_long_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata_long_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata_long_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata_long_zyam[itrg][jass]->Draw("PE"); 
      hdata3_long_zyam[itrg][jass]=GetJetCorrFunc1D_ZYAM_pp(hdata3_long[itrg][jass]);
      hdata3_long_zyam[itrg][jass]->SetMarkerStyle(24);
      hdata3_long_zyam[itrg][jass]->SetMarkerColor(2);
      hdata3_long_zyam[itrg][jass]->GetXaxis()->SetTitleSize(hdata3_long_zyam[itrg][jass]->GetXaxis()->GetTitleSize()*1.2);
      hdata3_long_zyam[itrg][jass]->GetYaxis()->SetTitleSize(hdata3_long_zyam[itrg][jass]->GetYaxis()->GetTitleSize()*1.13);
      hdata3_long_zyam[itrg][jass]->GetXaxis()->SetLabelSize(hdata3_long_zyam[itrg][jass]->GetXaxis()->GetLabelSize()*1.4);
      hdata3_long_zyam[itrg][jass]->GetYaxis()->SetLabelSize(hdata3_long_zyam[itrg][jass]->GetYaxis()->GetLabelSize()*1.4);
      hdata3_long_zyam[itrg][jass]->Draw("PESAME");

      if(itrg<ntrgptbins-1)
      {
        cc3->cd(itrg+1);
        hdata4_short_zyam[itrg][jass]->SetMarkerSize(1.8);
        hdata3_short_zyam[itrg][jass]->SetMarkerSize(1.5);
        hdata4_short_zyam[itrg][jass]->SetMarkerStyle(21);
        hdata3_short_zyam[itrg][jass]->SetMarkerStyle(20);
        hdata4_short_zyam[itrg][jass]->SetAxisRange(-0.02,0.12,"Y");
        hdata4_short_zyam[itrg][jass]->Draw("PE");
        hdata3_short_zyam[itrg][jass]->Draw("PESAME");
        cc3->cd(ntrgptbins-1+itrg+1);
        hdata2_short_zyam[itrg][jass]->SetMarkerSize(1.8);
        hdata_short_zyam[itrg][jass]->SetMarkerSize(1.5);
        hdata2_short_zyam[itrg][jass]->SetAxisRange(-0.04,0.37,"Y");
        hdata2_short_zyam[itrg][jass]->Draw("PE");
        hdata_short_zyam[itrg][jass]->Draw("PESAME");
        fitfunc2_vn_zyam[itrg][jass]->SetLineStyle(7);
        fitfunc_vn_zyam[itrg][jass]->SetLineStyle(7);
        fitfunc2_vn_zyam[itrg][jass]->SetLineColor(4);
        fitfunc_vn_zyam[itrg][jass]->SetLineColor(2);
        fitfunc2_vn_zyam[itrg][jass]->SetLineWidth(2);
        fitfunc_vn_zyam[itrg][jass]->SetLineWidth(2);
        fitfunc2_vn_zyam[itrg][jass]->Draw("LSAME");
        fitfunc_vn_zyam[itrg][jass]->Draw("LSAME");
      }
    }
  
    TString pttitle[ntrgptbins] = {"0.3 < p_{T}^{trig} < 1 GeV/c","1 < p_{T}^{trig} < 2 GeV/c","2 < p_{T}^{trig} < 4 GeV/c","4 < p_{T}^{trig} < 6 GeV/c","6 < p_{T}^{trig} < 12 GeV/c"};
    ofstream fout("data/Figure3.txt");
    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout << pttitle[itrg] << ", 220 <= N_trk^offline < 260, |deta|>2, pPb" <<endl;
      for ( int i=1 ; i <=hdata_short_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata_short_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata_short_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout << fixed << setprecision(4) << hdata_short_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata_short_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata_short_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout << pttitle[itrg] << ", 0 <= N_trk^offline < 20, |deta|>2, pPb" <<endl;
      for ( int i=1 ; i <=hdata3_short_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata3_short_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata3_short_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout << fixed << setprecision(4) << hdata3_short_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata3_short_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata3_short_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout << pttitle[itrg] << ", 220 <= N_trk^offline < 260, |deta|<1, pPb" <<endl;
      for ( int i=1 ; i <=hdata_long_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata_long_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata_long_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout << fixed << setprecision(4) << hdata_long_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata_long_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata_long_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout << pttitle[itrg] << ", 0 <= N_trk^offline < 20, |deta|<1, pPb" <<endl;
      for ( int i=1 ; i <=hdata3_long_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata3_long_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata3_long_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout << fixed << setprecision(4) << hdata3_long_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata3_long_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata3_long_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout << endl<<endl<<endl<<endl;
    }
    fout.close();

    ofstream fout1("data/Figure4.txt");
    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout1 << pttitle[itrg] << ", 220 <= N_trk^offline < 260, |deta|>2, PbPb" <<endl;
      for ( int i=1 ; i <=hdata2_short_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata2_short_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata2_short_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout1 << fixed << setprecision(4) << hdata2_short_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata2_short_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata2_short_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout1 << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout1 << pttitle[itrg] << ", 0 <= N_trk^offline < 20, |deta|>2, PbPb" <<endl;
      for ( int i=1 ; i <=hdata4_short_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata4_short_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata4_short_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout1 << fixed << setprecision(4) << hdata4_short_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata4_short_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata4_short_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout1 << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout1 << pttitle[itrg] << ", 220 <= N_trk^offline < 260, |deta|<1, PbPb" <<endl;
      for ( int i=1 ; i <=hdata2_long_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata2_long_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata2_long_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout1 << fixed << setprecision(4) << hdata2_long_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata2_long_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata2_long_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout1 << endl<<endl<<endl<<endl;
    }

    for(int itrg = 0; itrg<ntrgptbins; itrg++)
    {
      fout1 << pttitle[itrg] << ", 0 <= N_trk^offline < 20, |deta|<1, PbPb" <<endl;
      for ( int i=1 ; i <=hdata4_long_zyam[itrg][jass]->GetNbinsX() ; i++) {
        if ( (hdata4_long_zyam[itrg][jass]->GetBinCenter(i) > -1.2) && ( hdata4_long_zyam[itrg][jass]->GetBinCenter(i) < 2*PI-1.2) )
          fout1 << fixed << setprecision(4) << hdata4_long_zyam[itrg][jass]->GetBinCenter(i) << " " << fixed << setprecision(4) << hdata4_long_zyam[itrg][jass]->GetBinContent(i) << " " << fixed << setprecision(4) << hdata4_long_zyam[itrg][jass]->GetBinError(i) << endl;
      }
      fout1 << endl<<endl<<endl<<endl;
    }
    fout1.close();

    cc1->cd(1);
    TLatex* cms = new TLatex();
    cms->SetNDC();
    cms->SetTextSize(1.35*cms->GetTextSize());
    cms->DrawLatex(0.29,0.75,"CMS PbPb  #sqrt{s_{NN}} = 2.76 TeV");
//    cms->DrawLatex(0.3,0.55,Form("%d #leq N_{trk}^{offline} < %d",nmin,nmax));
    cms->SetTextSize(1.1*cms->GetTextSize());
    cms->DrawLatex(0.32,0.63,"Long-range (|#Delta#eta|>2)");
    cc1->cd(6);
    cms->DrawLatex(0.32,0.75,"Short-range (|#Delta#eta|<1)");
    cc1->cd(2);
    cms->SetTextSize(cms->GetTextSize()*1.1);
    cms->DrawLatex(0.18,0.7,"1 < p_{T}^{assoc} < 2 GeV/c");

    TLatex* latex_pt = new TLatex();
    latex_pt->SetNDC();
    latex_pt->SetTextSize(1.7*latex_pt->GetTextSize());    
    cc1->cd(1);
    latex_pt->DrawLatex(0.3,0.905,"0.3 < p_{T}^{trig} < 1 GeV/c");
    cc1->cd(2);
    latex_pt->DrawLatex(0.2,0.905,"1 < p_{T}^{trig} < 2 GeV/c");
    cc1->cd(3);
    latex_pt->DrawLatex(0.2,0.905,"2 < p_{T}^{trig} < 4 GeV/c");
    cc1->cd(4);
    latex_pt->DrawLatex(0.2,0.905,"4 < p_{T}^{trig} < 6 GeV/c");
    cc1->cd(5);
    latex_pt->DrawLatex(0.18,0.905,"6 < p_{T}^{trig} < 12 GeV/c");

    cc2->cd(1);
    TLatex* cms2 = new TLatex();
    cms2->SetNDC();
    cms2->SetTextSize(1.4*cms2->GetTextSize());
    cms2->DrawLatex(0.29,0.75,"CMS pPb  #sqrt{s_{NN}} = 5.02 TeV");
    cms2->SetTextSize(1.1*cms2->GetTextSize());
    cms2->DrawLatex(0.32,0.63,"Long-range (|#Delta#eta|>2)");
    cc2->cd(6);
    cms2->DrawLatex(0.32,0.75,"Short-range (|#Delta#eta|<1)");
    cc2->cd(2);
    cms2->SetTextSize(cms2->GetTextSize()*1.1);
    cms2->DrawLatex(0.18,0.7,"1 < p_{T}^{assoc} < 2 GeV/c");

    TLatex* latex_pt2 = new TLatex();
    latex_pt2->SetNDC();
    latex_pt2->SetTextSize(1.7*latex_pt2->GetTextSize());
    cc2->cd(1);
    latex_pt2->DrawLatex(0.3,0.905,"0.3 < p_{T}^{trig} < 1 GeV/c");
    cc2->cd(2);
    latex_pt2->DrawLatex(0.2,0.905,"1 < p_{T}^{trig} < 2 GeV/c");
    cc2->cd(3);
    latex_pt2->DrawLatex(0.2,0.905,"2 < p_{T}^{trig} < 4 GeV/c");
    cc2->cd(4);
    latex_pt2->DrawLatex(0.2,0.905,"4 < p_{T}^{trig} < 6 GeV/c");
    cc2->cd(5);
    latex_pt2->DrawLatex(0.18,0.905,"6 < p_{T}^{trig} < 12 GeV/c");

    TLegend* legend = new TLegend(0.28,0.27,0.7,0.575);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.075);
    legend->AddEntry(hdata2_long_zyam[0][jass],Form("%d #leq N_{trk}^{offline} < %d",nmin,nmax),"P");
    legend->AddEntry(hdata4_long_zyam[0][jass],Form("%d #leq N_{trk}^{offline} < %d",nmin3,nmax3),"P");
    legend->AddEntry(fitfunc_vn_zyam[1][1],"Fourier fit","L");
    cc1->cd(1);
    legend->Draw("same");

    TLegend* legend = new TLegend(0.28,0.252,0.7,0.56);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.075);
    legend->AddEntry(hdata_long_zyam[0][jass],Form("%d #leq N_{trk}^{offline} < %d",nmin,nmax),"P");
    legend->AddEntry(hdata3_long_zyam[0][jass],Form("%d #leq N_{trk}^{offline} < %d",nmin3,nmax3),"P");
    legend->AddEntry(fitfunc_vn_zyam[1][1],"Fourier fit","L");
    cc2->cd(1);
    legend->Draw("same");

    TLatex* latex3_pt = new TLatex();
    latex3_pt->SetNDC();
    latex3_pt->SetTextSize(1.35*latex3_pt->GetTextSize());
    cc3->cd(1);
    latex3_pt->DrawLatex(0.35,0.86,"0.3 < p_{T}^{trig} < 1 GeV/c");
    cc3->cd(2);
    latex3_pt->SetTextSize(1.25*latex3_pt->GetTextSize());
    latex3_pt->DrawLatex(0.18,0.86,"1 < p_{T}^{trig} < 2 GeV/c");
    cc3->cd(3);
    latex3_pt->DrawLatex(0.18,0.86,"2 < p_{T}^{trig} < 4 GeV/c");
    cc3->cd(4);
    latex3_pt->DrawLatex(0.16,0.86,"4 < p_{T}^{trig} < 6 GeV/c");
    latex3_pt->SetTextSize(latex3_pt->GetTextSize()/1.25);

    cc3->cd(1);
    latex3_pt->DrawLatex(0.35,0.72,"N_{trk}^{offline} < 20");
    cc3->cd(5);
    latex3_pt->DrawLatex(0.35,0.72,"220 < N_{trk}^{offline} < 260");
    cc3->cd(2);
    latex3_pt->SetTextSize(1.25*latex3_pt->GetTextSize());
    latex3_pt->DrawLatex(0.18,0.66,"1 < p_{T}^{assoc} < 2 GeV/c");
    latex3_pt->DrawLatex(0.18,0.56,"|#Delta#eta| > 2");
    latex3_pt->SetTextSize(latex3_pt->GetTextSize()/1.25);

    TLegend* legend = new TLegend(0.28,0.31,0.7,0.575);
    legend->SetFillColor(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.065);
    legend->AddEntry(hdata3_short_zyam[1][1],"pPb #sqrt{s_{NN}} = 5.02 TeV","P");
    legend->AddEntry(hdata4_short_zyam[1][1],"PbPb #sqrt{s_{NN}} = 2.76 TeV","P");
    cc3->cd(1);
    legend->Draw("same");

    SaveCanvas(cc3,"pPb/corr","corr1D_pPb_N220260_proposal");

return;
    TCanvas* cc3 = new TCanvas("cc3","cc3",550,500);
    cc3->SetLeftMargin(0.21);
    hdata_short[1][1]->SetAxisRange(hdata_short[1][1]->GetMinimum()*0.97,hdata_short[1][1]->GetMaximum()*1.03,"Y");
    hdata_short[1][1]->Draw("PE");
    fitfunc_vn[1][1]->Draw("Lsame");
    TF1* fitfunc_v1 = new TF1("fitfunc_v1","[0]*(1+2*[1]*TMath::Cos(x))",-1.2,2*PI-1.2);
    TF1* fitfunc_v2 = new TF1("fitfunc_v2","[0]*(1+2*[1]*TMath::Cos(2*x))",-1.2,2*PI-1.2);
    TF1* fitfunc_v3 = new TF1("fitfunc_v3","[0]*(1+2*[1]*TMath::Cos(3*x))",-1.2,2*PI-1.2);
    fitfunc_v1->SetParameters(fitfunc_vn[1][1]->GetParameter(0),fitfunc_vn[1][1]->GetParameter(1));
    fitfunc_v2->SetParameters(fitfunc_vn[1][1]->GetParameter(0),fitfunc_vn[1][1]->GetParameter(2));
    fitfunc_v3->SetParameters(fitfunc_vn[1][1]->GetParameter(0),fitfunc_vn[1][1]->GetParameter(3));
    fitfunc_v1->SetLineColor(6);
    fitfunc_v2->SetLineColor(kBlue);
    fitfunc_v3->SetLineColor(kGreen+2);
    fitfunc_v1->SetLineStyle(2);
    fitfunc_v2->SetLineStyle(9);
    fitfunc_v3->SetLineStyle(5);   
    fitfunc_v1->Draw("Lsame");
    fitfunc_v2->Draw("Lsame");
    fitfunc_v3->Draw("Lsame");

    TCanvas* cc4 = new TCanvas("cc4","cc4",550,500);
    cc4->SetLeftMargin(0.21);
    hdata2_short[1][1]->SetAxisRange(hdata2_short[1][1]->GetMinimum()*0.97,hdata2_short[1][1]->GetMaximum()*1.03,"Y");
    hdata2_short[1][1]->Draw("PE");
    fitfunc2_vn[1][1]->Draw("Lsame");
    TF1* fitfunc2_v1 = new TF1("fitfunc2_v1","[0]*(1+2*[1]*TMath::Cos(x))",-1.2,2*PI-1.2);
    TF1* fitfunc2_v2 = new TF1("fitfunc2_v2","[0]*(1+2*[1]*TMath::Cos(2*x))",-1.2,2*PI-1.2);
    TF1* fitfunc2_v3 = new TF1("fitfunc2_v3","[0]*(1+2*[1]*TMath::Cos(3*x))",-1.2,2*PI-1.2);
    fitfunc2_v1->SetParameters(fitfunc2_vn[1][1]->GetParameter(0),fitfunc2_vn[1][1]->GetParameter(1));
    fitfunc2_v2->SetParameters(fitfunc2_vn[1][1]->GetParameter(0),fitfunc2_vn[1][1]->GetParameter(2));
    fitfunc2_v3->SetParameters(fitfunc2_vn[1][1]->GetParameter(0),fitfunc2_vn[1][1]->GetParameter(3));
    fitfunc2_v1->SetLineColor(6);
    fitfunc2_v2->SetLineColor(kRed);
    fitfunc2_v3->SetLineColor(kGreen+2);
    fitfunc2_v1->SetLineStyle(2);
    fitfunc2_v2->SetLineStyle(9);
    fitfunc2_v3->SetLineStyle(5);
    fitfunc2_v1->Draw("Lsame");
    fitfunc2_v2->Draw("Lsame");
    fitfunc2_v3->Draw("Lsame");

    TCanvas* cc33 = new TCanvas("cc33","cc33",550,500);
    cc33->SetLeftMargin(0.21);
    hdata_short_zyam[1][1]->GetXaxis()->SetTitleOffset(hdata_short_zyam[1][1]->GetXaxis()->GetTitleOffset()*0.6);
    hdata_short_zyam[1][1]->GetYaxis()->SetTitleOffset(hdata_short_zyam[1][1]->GetYaxis()->GetTitleOffset()*0.6);
    hdata_short_zyam[1][1]->SetAxisRange(-0.1,0.2,"Y");
    hdata_short_zyam[1][1]->Draw("PE");

    TCanvas* cc44 = new TCanvas("cc44","cc44",550,500);
    cc44->SetLeftMargin(0.21);
    hdata2_short_zyam[1][1]->GetXaxis()->SetTitleOffset(hdata2_short_zyam[1][1]->GetXaxis()->GetTitleOffset()*0.6);
    hdata2_short_zyam[1][1]->GetYaxis()->SetTitleOffset(hdata2_short_zyam[1][1]->GetYaxis()->GetTitleOffset()*0.6);
    hdata2_short_zyam[1][1]->SetAxisRange(-0.05,0.3,"Y");
    hdata2_short_zyam[1][1]->Draw("PE");

    TCanvas* ccc3 = new TCanvas("ccc3","ccc3",550,500);
    ccc3->SetLeftMargin(0.21);
//    hdata_long[1][1]->SetAxisRange(hdata_long[1][1]->GetMinimum()*0.97,hdata_long[1][1]->GetMaximum()*1.03,"Y");
    hdata_long_zyam[1][1]->GetXaxis()->SetTitleOffset(hdata_long_zyam[1][1]->GetXaxis()->GetTitleOffset()*0.6);
    hdata_long_zyam[1][1]->GetYaxis()->SetTitleOffset(hdata_long_zyam[1][1]->GetYaxis()->GetTitleOffset()*0.6);
    hdata_long_zyam[1][1]->SetAxisRange(-0.1,0.4,"Y");
    hdata_long_zyam[1][1]->Draw("PE");

    TCanvas* ccc4 = new TCanvas("ccc4","ccc4",550,500);
    ccc4->SetLeftMargin(0.21);
//    hdata2_long[1][1]->SetAxisRange(hdata2_long[1][1]->GetMinimum()*0.97,hdata2_long[1][1]->GetMaximum()*1.03,"Y");
    hdata2_long_zyam[1][1]->GetXaxis()->SetTitleOffset(hdata2_long_zyam[1][1]->GetXaxis()->GetTitleOffset()*0.6);
    hdata2_long_zyam[1][1]->GetYaxis()->SetTitleOffset(hdata2_long_zyam[1][1]->GetYaxis()->GetTitleOffset()*0.6);
    hdata2_long_zyam[1][1]->SetAxisRange(-0.05,0.3,"Y");
    hdata2_long_zyam[1][1]->Draw("PE");

    SaveCanvas(cc33,"pPb/corr","corr1D_pPb_N220260_pt11_long_zyam");
    SaveCanvas(ccc3,"pPb/corr","corr1D_pPb_N220260_pt11_short_zyam");
    SaveCanvas(cc3,"pPb/corr","corr1D_pPb_N220260_pt11_long_fourier");
    SaveCanvas(cc4,"pPb/corr","corr1D_PbPb_N220260_pt11_long_fourier");
    SaveCanvas(cc33,"pPb/corr","corr1D_pPb_N220260_pt11_long");
    SaveCanvas(cc44,"pPb/corr","corr1D_PbPb_N220260_pt11_long");
    SaveCanvas(ccc3,"pPb/corr","corr1D_pPb_N220260_pt11_short");
    SaveCanvas(ccc4,"pPb/corr","corr1D_PbPb_N220260_pt11_short");

//0.96
/*
    cc1->cd(3);
    TLegend* legend1 = new TLegend(0.0,0.55,0.9,0.8);
    legend1->SetFillColor(0);
    legend1->SetFillStyle(0);
    legend1->AddEntry(hdata_long_zyam[0][jass],"pPb  #sqrt{s_{NN}} = 5.02 TeV","P");
    legend1->AddEntry(hdata2_long_zyam[0][jass],"PbPb  #sqrt{s_{NN}} = 2.76 TeV","P");
    legend1->Draw("same");
*/
//  SaveCanvas(cc1,"pPb/corr","corr1D_PbPb_N220260_paper");
//  SaveCanvas(cc2,"pPb/corr","corr1D_pPb_N220260_paper");

    TCanvas* cc3_zyam = new TCanvas("cc3_zyam","cc3_zyam",550,500);
    cc3_zyam->SetLeftMargin(0.21);
    hdata_short_zyam[1][1]->SetAxisRange(-0.03,0.18,"Y");
    hdata_short_zyam[1][1]->GetXaxis()->SetTitleOffset(hdata_short_zyam[1][1]->GetXaxis()->GetTitleOffset()*0.6);
    hdata_short_zyam[1][1]->GetYaxis()->SetTitleOffset(hdata_short_zyam[1][1]->GetYaxis()->GetTitleOffset()*0.75);
    hdata_short_zyam[1][1]->Draw("PE");
    SaveCanvas(cc3_zyam,"pPb/corr","corr1D_pPb_N220260_pt11_long_zyam");
}


void yieldvspt()
{
  const int nptbins = 7;
/*
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag = "HM_PromptRecoNew_Reverse_INCLEFF1v6";
*/

//  TString dataname = "PythiaReco_Minbias_7TeV";
  TString dataname = "PPData_Minbias_7TeV";
//  TString datatag = "TuneZ2_MinBias_42x_INCLGENv4";
//  TString datatag = "TuneZ2_MinBias_42x_INCLGenTrack032v7";
//  TString datatag = "TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8";
//  TString datatag = "TuneZ2_MinBias_42x_INCLGenTrackP033HFMPtWeight032";
  TString datatag = "TRIGGER_All_Apr21ReReco_INCLKshortPHFMPtWeight032v1";
//  TString datatag = "MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v9";
//  TString datatag = "TRIGGER_All_Apr21ReReco_INCLKshortPHFMPtWeight032v9";
//  TString datatag = "TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032";
//  TString datatag_periph = "MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight00Infv3";
//  TString datatag_periph = "TuneZ2_MinBias_42x_INCLTrackPHFMPtWeight00InfAllPt";
  TString datatag_periph = "TuneZ2_MinBias_42x_INCLGENv4";

/*
  TString dataname = "HIData_Minbias_2760GeV";
  TString datatag = "PPRereco_INCLEFF1v6HydjetEff";
  TString datatag_periph = "PPRereco_INCLEFF1v6";
*/
/*
  TString dataname = "HIJINGGen_Minbias_5TeV";
  TString datatag = "MBpPb_INCLv61";
  TString datatag_periph = "MBpPb_INCLv61";
*/
  TH1D *hdata_short[nptbins];
  TH1D *hdata_short_zyam[nptbins];
  TH1D *hdata_long[nptbins];
  TH1D *hdata_long_zyam[nptbins];
  TF1  *fitfunc_short[nptbins];
  TF1  *fitfunc_long[nptbins];
  TF1  *fitfunc_vn[nptbins];

  TH1D *hdata_short_periph[nptbins];
  TH1D *hdata_short_zyam_periph[nptbins];
  TH1D *hdata_long_periph[nptbins];
  TH1D *hdata_long_zyam_periph[nptbins];
  TF1  *fitfunc_short_periph[nptbins];
  TF1  *fitfunc_long_periph[nptbins];
  TF1  *fitfunc_vn_periph[nptbins];

  double etatrgmin=0.0;
  double etatrgmax=2.4;
  double etaassmin=-5.0;
  double etaassmax=-3.0;

/*
  double etatrgmin=-2.4;
  double etatrgmax=0.0;
  double etaassmin=3.0;
  double etaassmax=5.0;
*/
/*
  double etatrgmin=-5.0;
  double etatrgmax=-3.0;
  double etaassmin=3.0;
  double etaassmax=5.0;
*/
  int nmultmin=120;
  int nmultmax=1000;
//  int nmultmin=10;
//  int nmultmax=35;

  int nmultmin_periph=0;
  int nmultmax_periph=20;
  int ptbin_ref=0;
  bool IsSubPeriph=0;
  double scale=1.0;

  TGraphErrors* gr_yield_short_pt = new TGraphErrors(nptbins);
  TGraphErrors* gr_yield_long_pt = new TGraphErrors(nptbins);

//  TString filename_ref = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight053",nmultmin,nmultmax,-1,-1,-5.0,-3.0,3.0,5.0);
//  TString filename_ref = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032",150,1000,-1,-1,-5.0,-3.0,3.0,5.0);
//  TString filename_ref = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TuneZ2_MinBias_42x_INCLHFPHFMPtWeight032",0,20,-1,-1,-5.0,-3.0,3.0,5.0);
//  TString filename_ref = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TRIGGER_All_Apr21ReReco_INCLTrack032v7",150,1000,-1,-1,-2.4,2.4,-2.4,2.4);
  TString filename_ref = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TRIGGER_All_Apr21ReReco_INCLKshortPHFMPtWeight032v9",120,150,-1,-1,0.0,2.4,-5.0,-3.0);

  TH1D* hdata_short_ref = (TH1D*)Get1DCFdPhiv3(filename_ref.Data(),0,0,-10,10);
  TH1D* hdata_long_ref = (TH1D*)Get1DCFdPhiv3(filename_ref.Data(),0,0,-10,10);
  TH1D* hdata_short_ref_zyam=GetJetCorrFunc1D_ZYAM_pp(hdata_short_ref);
  TH1D* hdata_long_ref_zyam=GetJetCorrFunc1D_ZYAM_pp(hdata_long_ref);
  TF1* fitfunc_short_ref=GetFitFunc_ZYAM_pp(hdata_short_ref);
  TF1* fitfunc_long_ref=GetFitFunc_ZYAM_pp(hdata_long_ref);
  TF1* fitfunc_vn_ref = FitVnFunc(hdata_long_ref);

  double yielderror_long_ref=0;
  double yield_long_ref = hdata_long_ref_zyam->IntegralAndError(hdata_long_ref_zyam->FindBin(-1.2),hdata_long_ref_zyam->FindBin(1.2),yielderror_long_ref,"width");//zyamminimumscale;
  yielderror_long_ref=sqrt(yielderror_long_ref*yielderror_long_ref+fitfunc_long_ref->GetParError(0)*fitfunc_long_ref->GetParError(0)*2.4*2.4);

  double yielderror_short_ref=0;
  double yield_short_ref = hdata_short_ref_zyam->IntegralAndError(hdata_short_ref_zyam->FindBin(-1.2),hdata_short_ref_zyam->FindBin(1.2),yielderror_short_ref,"width");//zyamminimumscale;
  yielderror_short_ref=sqrt(yielderror_short_ref*yielderror_short_ref+fitfunc_short_ref->GetParError(0)*fitfunc_short_ref->GetParError(0)*2.4*2.4);

  double V0_ref = fitfunc_vn_ref->GetParameter(0);
  double V1_ref = fitfunc_vn_ref->GetParameter(1);
  double V2_ref = fitfunc_vn_ref->GetParameter(2);
  double V3_ref = fitfunc_vn_ref->GetParameter(3);
  double V0_ref_err = fitfunc_vn_ref->GetParError(0);
  double V1_ref_err = fitfunc_vn_ref->GetParError(1);
  double V2_ref_err = fitfunc_vn_ref->GetParError(2);
  double V3_ref_err = fitfunc_vn_ref->GetParError(3);

  if(IsSubPeriph)
  {
    TH1D* hdata_short_ref_periph;
    TH1D* hdata_long_ref_periph;
    TH1D* hdata_short_zyam_ref_periph;
    TH1D* hdata_long_zyam_ref_periph;
    TF1* fitfunc_short_ref_periph;
    TF1* fitfunc_long_ref_periph;
    TF1* fitfunc_vn_ref_periph;

//    TString filename_ref_periph = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TuneZ2_MinBias_42x_INCLHFPHFMPtWeight00Inf",nmultmin_periph,nmultmax_periph,-1,-1,-5.0,-3.0,3.0,5.0);
//    TString filename_ref_periph = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"MINBIAS_All_Apr21ReReco_INCLHFPHFMv2",nmultmin_periph,nmultmax_periph,-1,-1,-5.0,-3.0,3.0,5.0);
    TString filename_ref_periph = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),"TuneZ2_MinBias_42x_INCLHFPHFMPtWeight00Inf",0,20,-1,-1,-5.0,-3.0,3.0,5.0);

    hdata_short_ref_periph = (TH1D*)Get1DCFdPhiv3(filename_ref_periph.Data(),0,0,-10,10);
    hdata_long_ref_periph = (TH1D*)Get1DCFdPhiv3(filename_ref_periph.Data(),0,0,-10,10);

    hdata_short_ref->Add(hdata_short_ref_periph,-1*scale);
    hdata_long_ref->Add(hdata_long_ref_periph,-1*scale);
    hdata_short_ref_zyam=GetJetCorrFunc1D_ZYAM_pp(hdata_short_ref);
    hdata_long_ref_zyam=GetJetCorrFunc1D_ZYAM_pp(hdata_long_ref);
    fitfunc_short_ref=GetFitFunc_ZYAM_pp(hdata_short_ref);
    fitfunc_long_ref=GetFitFunc_ZYAM_pp(hdata_long_ref);
    fitfunc_vn_ref = FitVnFunc(hdata_long_ref);

/*
    fitfunc_short_ref_periph=GetFitFunc_ZYAM_pp(hdata_short_ref_periph);
    fitfunc_long_ref_periph=GetFitFunc_ZYAM_pp(hdata_long_ref_periph);
    fitfunc_vn_ref_periph = FitVnFunc(hdata_long_ref_periph);
    hdata_short_zyam_ref_periph=GetJetCorrFunc1D_ZYAM_pp(hdata_short_ref_periph);
    hdata_long_zyam_ref_periph=GetJetCorrFunc1D_ZYAM_pp(hdata_long_ref_periph);

    double yielderror_long_ref_periph=0;
    double yield_long_ref_periph = hdata_long_zyam_ref_periph->IntegralAndError(hdata_long_zyam_ref_periph->FindBin(-1.2),hdata_long_zyam_ref_periph->FindBin(1.2),yielderror_long_ref_periph,"width");//zyamminimumscale;
    yielderror_long_ref_periph=sqrt(yielderror_long_ref_periph*yielderror_long_ref_periph+fitfunc_long_ref_periph->GetParError(0)*fitfunc_long_ref_periph->GetParError(0)*2.4*2.4);

    double yielderror_short_ref_periph=0;
    double yield_short_ref_periph = hdata_short_zyam_ref_periph->IntegralAndError(hdata_short_zyam_ref_periph->FindBin(-1.2),hdata_short_zyam_ref_periph->FindBin(1.2),yielderror_short_ref_periph,"width");//zyamminimumscale;
    yielderror_short_ref_periph=sqrt(yielderror_short_ref_periph*yielderror_short_ref_periph+fitfunc_short_ref_periph->GetParError(0)*fitfunc_short_ref_periph->GetParError(0)*2.4*2.4);

    double V2_ref_periph=fitfunc_vn_ref_periph->GetParameter(2);
    double V3_ref_periph=fitfunc_vn_ref_periph->GetParameter(3);
    double V2_ref_periph_err=fitfunc_vn_ref_periph->GetParError(2);
    double V3_ref_periph_err=fitfunc_vn_ref_periph->GetParameter(3);
    double V0_ref_periph=fitfunc_vn_ref_periph->GetParameter(0);
    double V0_ref_periph_err=fitfunc_vn_ref_periph->GetParError(0);
    double V1_ref_periph=fitfunc_vn_ref_periph->GetParameter(1);
    double V1_ref_periph_err=fitfunc_vn_ref_periph->GetParError(1);

//    double amp_factor=(yield_short_ref-yield_long_ref)/(yield_short_ref_periph-yield_long_ref_periph);
    double amp_factor=1;

    V1_ref=V1_ref-amp_factor*V1_ref_periph*V0_ref_periph/V0_ref;
    V2_ref=V2_ref-amp_factor*V2_ref_periph*V0_ref_periph/V0_ref;
    V3_ref=V3_ref-amp_factor*V3_ref_periph*V0_ref_periph/V0_ref;
    V1_ref_err=sqrt(V1_ref_err*V1_ref_err+V1_ref_periph_err*V1_ref_periph_err*V0_ref_periph*V0_ref_periph/V0_ref/V0_ref*amp_factor*amp_factor);
    V2_ref_err=sqrt(V2_ref_err*V2_ref_err+V2_ref_periph_err*V2_ref_periph_err*V0_ref_periph*V0_ref_periph/V0_ref/V0_ref*amp_factor*amp_factor);
    V3_ref_err=sqrt(V3_ref_err*V3_ref_err+V3_ref_periph_err*V3_ref_periph_err*V0_ref_periph*V0_ref_periph/V0_ref/V0_ref*amp_factor*amp_factor);
*/
    V0_ref = fitfunc_vn_ref->GetParameter(0);
    V1_ref = fitfunc_vn_ref->GetParameter(1);
    V2_ref = fitfunc_vn_ref->GetParameter(2);
    V3_ref = fitfunc_vn_ref->GetParameter(3);
    V0_ref_err = fitfunc_vn_ref->GetParError(0);
    V1_ref_err = fitfunc_vn_ref->GetParError(1);
    V2_ref_err = fitfunc_vn_ref->GetParError(2);
    V3_ref_err = fitfunc_vn_ref->GetParError(3);
  }
  double v2_ref = sqrt(fabs(V2_ref));
  double v3_ref = sqrt(fabs(V3_ref));
  double v2_ref_err = V2_ref_err/V2_ref*v2_ref;
  double v3_ref_err = V3_ref_err/V3_ref*v3_ref;

  double v2[nptbins];
  double v3[nptbins];
  double v2_err[nptbins];
  double v3_err[nptbins];
  double ptbins[nptbins];
  double ptyield[nptbins];
  double V0[nptbins];
  double V1[nptbins];
  double V2[nptbins];
  double V3[nptbins];
  double V0_err[nptbins];
  double V1_err[nptbins];
  double V2_err[nptbins];
  double V3_err[nptbins];
  double yield_short[nptbins];
  double yield_long[nptbins];
  double yielderror_short[nptbins];
  double yielderror_long[nptbins];

  for(int ibin = 0; ibin<nptbins; ibin++)
  {
    TString filename = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),datatag.Data(),nmultmin,nmultmax,-1,-1,etatrgmin,etatrgmax,etaassmin,etaassmax);
    hdata_short[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ibin,ptbin_ref,-10.0,10.0);
    hdata_long[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ibin,ptbin_ref,-10.0,10.0);
//    hdata_short[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ibin,ibin,-10.0,10.0);
//    hdata_long[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ibin,ibin,-10.0,10.0);

    fitfunc_short[ibin]=GetFitFunc_ZYAM_pp(hdata_short[ibin]);
    fitfunc_long[ibin]=GetFitFunc_ZYAM_pp(hdata_long[ibin]);
    fitfunc_vn[ibin] = FitVnFunc(hdata_long[ibin]);
    hdata_short_zyam[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_short[ibin]);
    hdata_long_zyam[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_long[ibin]);
    TH1D* hpt = (TH1D*)GetHist(filename,Form("ptcorr_signal_trg_%d",ibin));

    yielderror_long[ibin]=0;
    yield_long[ibin] = hdata_long_zyam[ibin]->IntegralAndError(hdata_long_zyam[ibin]->FindBin(-1.2),hdata_long_zyam[ibin]->FindBin(1.2),yielderror_long[ibin],"width");//zyamminimumscale;
 //   yielderror_short=yielderror_short/2;
    yielderror_long[ibin]=sqrt(yielderror_long[ibin]*yielderror_long[ibin]+fitfunc_long[ibin]->GetParError(0)*fitfunc_long[ibin]->GetParError(0)*2.4*2.4);
    gr_yield_long_pt->SetPoint(ibin,hpt->GetMean(),yield_long[ibin]);
    gr_yield_long_pt->SetPointError(ibin,0,yielderror_long[ibin]);

    yielderror_short[ibin]=0;
    yield_short[ibin] = hdata_short_zyam[ibin]->IntegralAndError(hdata_short_zyam[ibin]->FindBin(-1.2),hdata_short_zyam[ibin]->FindBin(1.2),yielderror_short[ibin],"width");//zyamminimumscale;
 //   yielderror_short=yielderror_short/2;
    yielderror_short[ibin]=sqrt(yielderror_short[ibin]*yielderror_short[ibin]+fitfunc_short[ibin]->GetParError(0)*fitfunc_short[ibin]->GetParError(0)*2.4*2.4);
    gr_yield_short_pt->SetPoint(ibin,hpt->GetMean(),(yield_short[ibin]-yield_long[ibin]));
    gr_yield_short_pt->SetPointError(ibin,0,sqrt(yielderror_short[ibin]*yielderror_short[ibin]+yielderror_long[ibin]*yielderror_long[ibin]));
//    gr_yield_short_pt->SetPoint(ibin,hpt->GetMean(),yield_short);
//    gr_yield_short_pt->SetPointError(ibin,0,yielderror_short);

    ptbins[ibin]=hpt->GetMean();
    ptyield[ibin]=hpt->Integral();

    V2[ibin]=fitfunc_vn[ibin]->GetParameter(2);
    V3[ibin]=fitfunc_vn[ibin]->GetParameter(3);
    V2_err[ibin]=fitfunc_vn[ibin]->GetParError(2);
    V3_err[ibin]=fitfunc_vn[ibin]->GetParError(3);
    V0[ibin]=fitfunc_vn[ibin]->GetParameter(0);
    V0_err[ibin]=fitfunc_vn[ibin]->GetParError(0);
    V1[ibin]=fitfunc_vn[ibin]->GetParameter(1);
    V1_err[ibin]=fitfunc_vn[ibin]->GetParError(1);

//    V0[ibin]=fitfunc_vn[ibin]->GetParameter(0)*V2[ibin];
//    V0_err[ibin]=sqrt((fitfunc_vn[ibin]->GetParError(0)/fitfunc_vn[ibin]->GetParameter(0))**2+(V2_err[ibin]/V2[ibin])**2)*V0[ibin];

    if(IsSubPeriph)
    {
      TString filename_periph = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),datatag_periph.Data(),nmultmin_periph,nmultmax_periph,-1,-1,0.0,2.4,-5.0,-3.0);
      hdata_short_periph[ibin] = (TH1D*)Get1DCFdPhiv3(filename_periph.Data(),ibin,ptbin_ref,-10,10);
      hdata_long_periph[ibin] = (TH1D*)Get1DCFdPhiv3(filename_periph.Data(),ibin,ptbin_ref,-10,10);
//      hdata_short_periph[ibin] = (TH1D*)Get1DCFdPhiv3(filename_periph.Data(),ibin,ibin,-10,10);
//      hdata_long_periph[ibin] = (TH1D*)Get1DCFdPhiv3(filename_periph.Data(),ibin,ibin,-10,10);

      hdata_short[ibin]->Add(hdata_short_periph[ibin],-1*scale);
      hdata_long[ibin]->Add(hdata_short_periph[ibin],-1*scale);
      fitfunc_short[ibin]=GetFitFunc_ZYAM_pp(hdata_short[ibin]);
      fitfunc_long[ibin]=GetFitFunc_ZYAM_pp(hdata_long[ibin]);
      fitfunc_vn[ibin] = FitVnFunc(hdata_long[ibin]);

/*
      fitfunc_short_periph[ibin]=GetFitFunc_ZYAM_pp(hdata_short_periph[ibin]);
      fitfunc_long_periph[ibin]=GetFitFunc_ZYAM_pp(hdata_long_periph[ibin]);
      fitfunc_vn_periph[ibin] = FitVnFunc(hdata_long_periph[ibin]);
      hdata_short_zyam_periph[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_short_periph[ibin]);
      hdata_long_zyam_periph[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_long_periph[ibin]);

      double yielderror_long_periph=0;
      double yield_long_periph = hdata_long_zyam_periph[ibin]->IntegralAndError(hdata_long_zyam_periph[ibin]->FindBin(-1.2),hdata_long_zyam_periph[ibin]->FindBin(1.2),yielderror_long_periph,"width");//zyamminimumscale;
      yielderror_long_periph=sqrt(yielderror_long_periph*yielderror_long_periph+fitfunc_long_periph[ibin]->GetParError(0)*fitfunc_long_periph[ibin]->GetParError(0)*2.4*2.4);

      double yielderror_short_periph=0;
      double yield_short_periph = hdata_short_zyam_periph[ibin]->IntegralAndError(hdata_short_zyam_periph[ibin]->FindBin(-1.2),hdata_short_zyam_periph[ibin]->FindBin(1.2),yielderror_short_periph,"width");//zyamminimumscale;
      yielderror_short_periph=sqrt(yielderror_short_periph*yielderror_short_periph+fitfunc_short_periph[ibin]->GetParError(0)*fitfunc_short_periph[ibin]->GetParError(0)*2.4*2.4);

      double V2_periph=fitfunc_vn_periph[ibin]->GetParameter(2);
      double V3_periph=fitfunc_vn_periph[ibin]->GetParameter(3);
      double V2_periph_err=fitfunc_vn_periph[ibin]->GetParError(2);
      double V3_periph_err=fitfunc_vn_periph[ibin]->GetParameter(3);
      double V0_periph=fitfunc_vn_periph[ibin]->GetParameter(0);
      double V0_periph_err=fitfunc_vn_periph[ibin]->GetParError(0);
      double V1_periph=fitfunc_vn_periph[ibin]->GetParameter(1);
      double V1_periph_err=fitfunc_vn_periph[ibin]->GetParError(1);

      double amp_factor=(yield_short[ibin]-yield_long[ibin])/(yield_short_periph-yield_long_periph);
//cout<<"amplification factor:"<<amp_factor<<endl;
//      double amp_factor=2;

      V1[ibin]=V1[ibin]-amp_factor*V1_periph*V0_periph/V0[ibin];
      V2[ibin]=V2[ibin]-amp_factor*V2_periph*V0_periph/V0[ibin];
      V3[ibin]=V3[ibin]-amp_factor*V3_periph*V0_periph/V0[ibin];
      V1_err[ibin]=sqrt(V1_err[ibin]*V1_err[ibin]+V1_periph_err*V1_periph_err*V0_periph*V0_periph/V0[ibin]/V0[ibin]*amp_factor*amp_factor);
      V2_err[ibin]=sqrt(V2_err[ibin]*V2_err[ibin]+V2_periph_err*V2_periph_err*V0_periph*V0_periph/V0[ibin]/V0[ibin]*amp_factor*amp_factor);
      V3_err[ibin]=sqrt(V3_err[ibin]*V3_err[ibin]+V3_periph_err*V3_periph_err*V0_periph*V0_periph/V0[ibin]/V0[ibin]*amp_factor*amp_factor);
*/
    V2[ibin]=fitfunc_vn[ibin]->GetParameter(2);
    V3[ibin]=fitfunc_vn[ibin]->GetParameter(3);
    V2_err[ibin]=fitfunc_vn[ibin]->GetParError(2);
    V3_err[ibin]=fitfunc_vn[ibin]->GetParError(3);
    V0[ibin]=fitfunc_vn[ibin]->GetParameter(0);//*V2[ibin];
    V0_err[ibin]=fitfunc_vn[ibin]->GetParError(0);//sqrt((fitfunc_vn[ibin]->GetParError(0)/fitfunc_vn[ibin]->GetParameter(0))**2+(V2_err[ibin]/V2[ibin])**2)*V0[ibin];
    V1[ibin]=fitfunc_vn[ibin]->GetParameter(1);
    V1_err[ibin]=fitfunc_vn[ibin]->GetParError(1);
    }

// Fourier decomposition
    v2[ibin] = V2[ibin]/v2_ref;
    v3[ibin] = V3[ibin]/v3_ref;
    v2_err[ibin] = sqrt((V2_err[ibin]/V2[ibin])**2+(V2_ref_err/V2_ref)**2)*fabs(v2[ibin]);
    v3_err[ibin] = sqrt((V3_err[ibin]/V3[ibin])**2+(V3_ref_err/V3_ref)**2)*fabs(v3[ibin]);
//cout<<V3_err<<" "<<V3<<" "<<V3_ref_err<<" "<<V3_ref<<endl;
  }

  TCanvas* c1D_short = new TCanvas("c1D_short","c1D_short",800,800);
  c1D_short->Divide(4,4);
  TCanvas* c1D_long = new TCanvas("c1D_long","c1D_long",800,800);
  c1D_long->Divide(4,4);

  for(int i=1;i<=nptbins;i++) 
  {
    c1D_short->cd(i);
    hdata_short[i-1]->Draw("PE");
    fitfunc_short[i-1]->SetLineColor(2);
    fitfunc_short[i-1]->Draw("Lsame");
//    hdata_short_ref_zyam->Draw("PESAME");
    c1D_long->cd(i);
    hdata_long[i-1]->Draw("PE");
    fitfunc_long[i-1]->SetLineColor(2);
    fitfunc_long[i-1]->Draw("Lsame");
    fitfunc_vn[i-1]->Draw("Lsame");
//    hdata_long_ref_zyam->Draw("PESAME");
  }

  double ptbins_err[nptbins]={0};
  TGraphErrors* gr_v0 = new TGraphErrors(nptbins,ptbins,V0,ptbins_err,V0_err);
  TGraphErrors* gr_v2 = new TGraphErrors(nptbins,ptbins,V2,ptbins_err,V2_err);
  TGraphErrors* gr_v3 = new TGraphErrors(nptbins,ptbins,V3,ptbins_err,V3_err);
  TH2D* htmpv2 = new TH2D("tmpv2",";p_{T} (GeV/c);v_{2}",1000,0,20.0,1000,-0.15,0.5);
  TH2D* htmpv3 = new TH2D("tmpv3",";p_{T} (GeV/c);v_{3}",1000,0,6.0,1000,-0.1,0.1);

  TCanvas* c = new TCanvas("c","",550,500); 
  htmpv2->Draw();
  gr_v2->Draw("PESAME");

  TCanvas* c1 = new TCanvas("c1","",550,500);
  htmpv3->Draw();
  gr_v3->Draw("PESAME");

  TFile* ffout = new TFile(Form("results/%s_%s_nmin%d_nmax%d.root",dataname.Data(),datatag.Data(),nmultmin,nmultmax),"recreate");
  gr_v0->Write("graph_v0");
  gr_v2->Write("graph_v2");
  gr_v3->Write("graph_v3");
  ffout->Close();
  ofstream fcout(Form("data/Nv2vspt_%s_N%d%d_PythiaTrackHF.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d_trigger160.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d_nvtx1.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d.txt",dataname.Data(),nmultmin,nmultmax));
  fcout << "pt    v_0*v_2    v_0*v_2 err." << endl;
  for(int i = 0; i<nptbins; i++) fcout << ptbins[i] << " "<< V0[i] << " " << V0_err[i] <<endl;

return;
  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d_TrackHF.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d_trigger160.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d_nvtx1.txt",dataname.Data(),nmultmin,nmultmax));
//  ofstream fcout(Form("data/v2v3vspt_%s_N%d%d.txt",dataname.Data(),nmultmin,nmultmax));
  fcout << "pt    v_2    v_2 err." << endl;
  for(int i = 0; i<nptbins; i++) fcout << ptbins[i] << " "<< v2[i] << " " << v2_err[i] <<endl;
  fcout<<endl<<endl;
  fcout << "pt    v_3    v_3 err." << endl;
  for(int i = 0; i<nptbins; i++) fcout << ptbins[i] << " "<< v3[i] << " " << v3_err[i] <<endl;
  fcout<<endl<<endl; 

  ofstream fcout(Form("data/NMultV2V3vspt_%s_N%d%d_TrackHF.txt",dataname.Data(),nmultmin,nmultmax));
  fcout << "pt    V_0 " << endl;
  for(int i = 0; i<nptbins; i++) fcout <<ptbins[i] << " "<< V0[i] << endl;
  fcout<<endl<<endl;
  fcout << "pt    V_1*NMult    V_1*NMult err." << endl;
  for(int i = 0; i<nptbins; i++) fcout <<ptbins[i] << " "<< V1[i]*V0[i] << " " << V1_err[i]*V0[i] <<endl;
  fcout<<endl<<endl;
  fcout << "pt    V_2*NMult    V_2*NMult err." << endl;
  for(int i = 0; i<nptbins; i++) fcout <<ptbins[i] << " "<< V2[i]*V0[i] << " " << V2_err[i]*V0[i] <<endl;
  fcout<<endl<<endl;
  fcout << "pt    V_3*NMult    V_3*NMult err." << endl;
  for(int i = 0; i<nptbins; i++) fcout <<ptbins[i] << " "<< V3[i]*V0[i] << " " << V3_err[i]*V0[i] <<endl;
  fcout<<endl<<endl;

  ofstream fcout(Form("data/JetYieldvspt_%s_N%d%d_TrackHF.txt",dataname.Data(),nmultmin,nmultmax));
  fcout << "pt    yield    yield err." << endl;
  for(int i = 0; i<nptbins; i++) fcout <<ptbins[i] << " "<< yield_short[i]-yield_long[i] << " " << sqrt(yielderror_short[i]*yielderror_short[i]+yielderror_long[i]*yielderror_long[i]) <<endl;
  fcout<<endl<<endl;

return;
  TFile* fout = new TFile(Form("results/pPb_yieldvspt_N%d%d_lambda.root",nmultmin,nmultmax),"recreate");
//  TFile* fout = new TFile(Form("results/pPb_yieldvspt_N%d%d_sub020.root",nmultmin,nmultmax),"recreate");
//  TFile* fout = new TFile(Form("results/pPb_yieldvspt_N%d%d_nvtx1.root",nmultmin,nmultmax),"recreate");
//  TFile* fout = new TFile(Form("results/pPb_yieldvspt_N%d%d.root",nmultmin,nmultmax),"recreate");
//  TFile* fout = new TFile(Form("results/PbPb_yieldvspt_N%d%d_allnew_hydjeteff.root",nmultmin,nmultmax),"recreate");
  gr_yield_short_pt->SetName("gr_yield_short_pt");
  gr_yield_long_pt->SetName("gr_yield_long_pt");
  gr_yield_short_pt->Write();
  gr_yield_long_pt->Write();
  fout->Close();
}

void yieldvsn()
{
/*
  const int nmultbin = 14;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"MBCent","MBCent","MBCent","MBCent","MBCent","MBCent","MBCent","HFHMCent","HFHMCent","HFHMCent","HFHMCent","HFHMCent","HFHMCent","HFHMCent"};
  for(int i=0;i<nmultbin;i++) 
    datatag[i].Append("_PromptRecoNew_Reverse_INCLEFF1ONEPTBINv1");
  int nmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int nmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmin[nmultbin] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};
  int centmax[nmultbin] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
*/
/*
  const int nmultbin = 10;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"MB","MB","MB","MB","HM","HM","HM","HM","HM","HM"};
  for(int i=0;i<nmultbin;i++) 
    if(i<4) datatag[i].Append("_PromptRecoNew_All_INCLEFF1ONEPTBINv1Tight");
    else datatag[i].Append("_PromptRecoNew_INCLEFF1ONEPTBINv1Tight");
  int nmin[nmultbin] = {0,20,50,80,120,150,185,220,260,300};
  int nmax[nmultbin] = {20,50,80,120,150,185,220,260,300,350};
*/
/*
  const int nmultbin = 1;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"HFHM"};
  for(int i=0;i<nmultbin;i++) 
    datatag[i].Append("_PromptRecoNew_Reverse_INCLEFF1ONEPTBINv1");       
  int nmin[nmultbin] = {-1};
  int nmax[nmultbin] = {-1};
*/
/*
  const int nmultbin = 6;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"HM","HM","HM","HM","HM","HM220"};
  for(int i=0;i<nmultbin;i++) 
    if(i==(nmultbin-1)) datatag[i].Append("_PromptRecoNew_All_INCLEFF1ONEPTBINv1VZ03");       
    else datatag[i].Append("_PromptRecoNew_INCLEFF1ONEPTBINv1VZ03");
  int nmin[nmultbin] = {120,150,185,220,260,300};
  int nmax[nmultbin] = {150,185,220,260,300,350};
*/

/*
  const int nmultbin = 16;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"MB","MB","MB","MB","MB","MB","MB","MB","MB","MB","HM","HM","HM","HM","HM","HM220"};
  for(int i=0;i<nmultbin;i++) 
  {
      datatag[i].Append("_PromptRecoNew_Reverse_INCLEFF1ONEPTBINv1");
  }      
  int nmin[nmultbin] = {0,5,10,20,30,40,50,60,80,100,120,150,185,220,260,300};
  int nmax[nmultbin] = {5,10,20,30,40,50,60,80,100,120,150,185,220,260,300,350};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
*/
/*
  const int nmultbin = 20;
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag[nmultbin] = {"MB","MB","MB","MB","MB","MB","MB","MB","MB","MB","HM","HM","HM","HM","HM","HM","HM","HM220","HM220","MB"};
  for(int i=0;i<nmultbin;i++) 
  {
//    if(i<nmultbin-9) datatag[i].Append("_PromptRecoNew_INCLEFF1ONEPTBINv1");                
//      datatag[i].Append("_PromptRecoNew_Reverse_INCLEFF1ONEPTBINv1");
    datatag[i].Append("_PromptRecoNew_All_INCLEFF1ONEPTBINv11");
//    datatag[i].Append("_PromptRecoNew_INCLEFF1ONEPTBINv12");
  }      
  int nmin[nmultbin] = {0,5,10,20,30,40,50,60,80,100,120,150,185,220,240,260,280,300,320,0};
  int nmax[nmultbin] = {5,10,20,30,40,50,60,80,100,120,150,185,220,240,260,280,300,320,350,35};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
*/
/*
  const int nmultbin = 20;
  TString dataname = "HIData_Minbias_2760GeV";
  TString datatag[nmultbin];
  for(int i=0;i<nmultbin;i++) datatag[i] = "PPRereco_INCLEFF1ONEPTBINv1HydjetEff";
  int nmin[nmultbin] = {0,5,10,20,30,40,50,60,80,100,120,150,185,220,240,260,280,300,320,0};
  int nmax[nmultbin] = {5,10,20,30,40,50,60,80,100,120,150,185,220,240,260,280,300,320,350,20};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
*/
/*
  const int nmultbin = 9;
  TString dataname = "PythiaReco_Minbias_7TeV";
  TString datatag[nmultbin];
//  for(int i=0;i<nmultbin;i++) datatag[i] = "TuneZ2_MinBias_42x_INCLGenTrackP033HFMPtWeight032";
  for(int i=0;i<nmultbin;i++) datatag[i] = "TuneZ2_MinBias_42x_INCLHFPHFMPtWeight032";

  int nmin[nmultbin] = {0,10,20,30,40,50,60,70,80};
  int nmax[nmultbin] = {10,20,30,40,50,60,70,80,120};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1};
*/

  const int nmultbin = 8;
  TString dataname = "PPData_Minbias_7TeV";
  TString datatag[nmultbin];
  for(int i=0;i<nmultbin;i++) 
  {  
//    datatag[i] = "MINBIAS_All_Apr21ReReco_INCLTrackP033HFMPtWeight032";
//    if(i>=4) datatag[i] = "TRIGGER_All_Apr21ReReco_INCLTrackP033HFMPtWeight032";
//    datatag[i] = "MINBIAS_All_Apr21ReReco_INCLTrackM13HFPPtWeight032";
//    if(i>=13) datatag[i] = "TRIGGER_All_Apr21ReReco_INCLTrackM13HFPPtWeight032";
    datatag[i] = "MINBIAS_All_Apr21ReReco_INCLHFPHFMPtWeight032";
    if(i>=4) datatag[i] = "TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032";
  }
  int nmin[nmultbin] = {0,10,35,60,85,100,120,140};
  int nmax[nmultbin] = {10,35,60,85,100,120,140,1000};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1};
/*
  double etatrgmin=0.0;
  double etatrgmax=2.4;
  double etaassmin=-5.0;
  double etaassmax=-3.0;
*/
/*
  double etatrgmin=-2.4;
  double etatrgmax=0.0;
  double etaassmin=3.0;
  double etaassmax=5.0;
*/
/*
  double etatrgmin=-5.2;
  double etatrgmax=-3.0;
  double etaassmin=3.0;
  double etaassmax=5.2;
*/
/*
  double etatrgmin=3.0;
  double etatrgmax=5.0;
  double etaassmin=3.0;
  double etaassmax=5.0;
*/
/*
  double etatrgmin=-5.0;
  double etatrgmax=-3.0;
  double etaassmin=-5.0;
  double etaassmax=-3.0;
*/

  double etatrgmin=-5.0;
  double etatrgmax=-3.0;
  double etaassmin=3.0;
  double etaassmax=5.0;

/*
  const int nmultbin = 12;
  TString dataname = "HIJINGGen_Minbias_5TeV";
  TString datatag[nmultbin];
  for(int i=0;i<nmultbin;i++) datatag[i] = "MBpPb_INCLONEPTBINv11ETA45";
//  int nmin[nmultbin] = {0,5,10,20,30,40,50,60,80,100,120,150,0};
//  int nmax[nmultbin] = {5,10,20,30,40,50,60,80,100,120,150,200,20};
  int nmin[nmultbin] = {0,2,5,7,10,15,20,25,30,35,40,0};
  int nmax[nmultbin] = {2,5,7,10,15,20,25,30,35,40,50,5};
  int centmin[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
  int centmax[nmultbin] = {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1};
*/
  TH1D *hdata_short[nmultbin];
  TH1D *hdata_short_zyam[nmultbin];
  TH1D *hdata_long[nmultbin];
  TH1D *hdata_long_zyam[nmultbin];
  TF1  *fitfunc_short[nmultbin];
  TF1  *fitfunc_long[nmultbin];
  TF1  *fitfunc_vn[nmultbin]; 

  double V2[nmultbin];
  double V3[nmultbin];
  double V2_err[nmultbin];
  double V3_err[nmultbin];
  double V0[nmultbin];
  double V0_err[nmultbin];
  double v2[nmultbin];
  double v3[nmultbin];
  double v2_err[nmultbin];
  double v3_err[nmultbin];
  double nmult[nmultbin];

  double yield_short[nmultbin];
  double yielderror_short[nmultbin];
  double yield_long[nmultbin];
  double yielderror_long[nmultbin];

  bool IsSubPeriph=0;
  int ptbin=0;
  int bin_periph=nmultbin-1;
//  int bin_periph=4;
  TGraphErrors* gr_yield_short_mult = new TGraphErrors(nmultbin-1);
  TGraphErrors* gr_yield_long_mult = new TGraphErrors(nmultbin-1);

  for(int ibin = nmultbin-1; ibin>=0; ibin--)
  { 
    TString filename = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),datatag[ibin].Data(),nmin[ibin],nmax[ibin],centmin[ibin],centmax[ibin],etatrgmin,etatrgmax,etaassmin,etaassmax);
    hdata_short[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ptbin,ptbin,-1,1);
    hdata_long[ibin] = (TH1D*)Get1DCFdPhiv3(filename.Data(),ptbin,ptbin,-12,12);
//    hdata_short[ibin]->Add(hdata_long[ibin],-1);
    hdata_short_zyam[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_short[ibin]);
    hdata_long_zyam[ibin]=GetJetCorrFunc1D_ZYAM_pp(hdata_long[ibin]);
    fitfunc_short[ibin]=GetFitFunc_ZYAM_pp(hdata_short[ibin]);
    fitfunc_long[ibin]=GetFitFunc_ZYAM_pp(hdata_long[ibin]);
// for pPb low stat.
//    if(ibin>=(nmultbin-2)) hdata_long_zyam[ibin]=SimpleZYAM1D(hdata_long[ibin],1.4);
//    if(ibin>=(nmultbin-1)) hdata_long_zyam[ibin]=SimpleZYAM1D(hdata_long[ibin],1.2);
//    if(ibin>=(nmultbin-1)) hdata_short_zyam[ibin]=SimpleZYAM1D(hdata_short[ibin],1.2);

    fitfunc_vn[ibin] = FitVnFunc(hdata_long[ibin]);
    TH1D* hmult = (TH1D*)GetHist(filename,Form("multrawall"));
    
    yield_long[ibin] = hdata_long_zyam[ibin]->IntegralAndError(hdata_long_zyam[ibin]->FindBin(-1.2),hdata_long_zyam[ibin]->FindBin(1.2),yielderror_long[ibin],"width");//zyamminimumscale;
    yielderror_long[ibin]=sqrt(yielderror_long[ibin]*yielderror_long[ibin]+fitfunc_long[ibin]->GetParError(0)*fitfunc_long[ibin]->GetParError(0)*2.4*2.4);
    gr_yield_long_mult->SetPoint(ibin,hmult->GetMean(),yield_long[ibin]);
    gr_yield_long_mult->SetPointError(ibin,0,yielderror_long[ibin]);

    yield_short[ibin] = hdata_short_zyam[ibin]->IntegralAndError(hdata_short_zyam[ibin]->FindBin(-1.2),hdata_short_zyam[ibin]->FindBin(1.2),yielderror_short[ibin],"width");//zyamminimumscale;
    yielderror_short[ibin]=sqrt(yielderror_short[ibin]*yielderror_short[ibin]+fitfunc_short[ibin]->GetParError(0)*fitfunc_short[ibin]->GetParError(0)*2.4*2.4);
    gr_yield_short_mult->SetPoint(ibin,hmult->GetMean(),yield_short[ibin]-yield_long[ibin]);
    gr_yield_short_mult->SetPointError(ibin,0,sqrt(yielderror_short[ibin]*yielderror_short[ibin]+yielderror_long[ibin]*yielderror_long[ibin]));


//    gr_yield_short_mult->SetPoint(ibin,hmult->GetMean(),yield_short);
//    gr_yield_short_mult->SetPointError(ibin,0,yielderror_short);

    nmult[ibin] = hmult->GetMean();

// Fourier decomposition
    V2[ibin]=fitfunc_vn[ibin]->GetParameter(2);//*fitfunc_vn[ibin]->GetParameter(0);
    V3[ibin]=fitfunc_vn[ibin]->GetParameter(3);//*fitfunc_vn[ibin]->GetParameter(0);
    V2_err[ibin]=fitfunc_vn[ibin]->GetParError(2);//*fitfunc_vn[ibin]->GetParameter(0);
    V3_err[ibin]=fitfunc_vn[ibin]->GetParError(3);//*fitfunc_vn[ibin]->GetParameter(0);
    V0[ibin]=fitfunc_vn[ibin]->GetParameter(0);
    V0_err[ibin]=fitfunc_vn[ibin]->GetParError(0);

    if(IsSubPeriph && ibin!=bin_periph) 
    {
      double amp_factor = (yield_short[ibin]-yield_long[ibin])/(yield_short[bin_periph]-yield_long[bin_periph]);
 //     double amp_factor = 1.0;

      V2[ibin]=V2[ibin]-amp_factor*V2[bin_periph]*V0[bin_periph]/V0[ibin];
      V3[ibin]=V3[ibin]-amp_factor*V3[bin_periph]*V0[bin_periph]/V0[ibin];
      V2_err[ibin]=sqrt(V2_err[ibin]*V2_err[ibin]+amp_factor*amp_factor*V2_err[bin_periph]*V2_err[bin_periph]*V0[bin_periph]*V0[bin_periph]/V0[ibin]/V0[ibin]);
      V3_err[ibin]=sqrt(V3_err[ibin]*V3_err[ibin]+amp_factor*amp_factor*V3_err[bin_periph]*V3_err[bin_periph]*V0[bin_periph]*V0[bin_periph]/V0[ibin]/V0[ibin]);
    }

    v2[ibin] = V2[ibin]/sqrt(fabs(V2[ibin]));
    v3[ibin] = V3[ibin]/sqrt(fabs(V3[ibin]));
    v2_err[ibin] = V2_err[ibin]/V2[ibin]*v2[ibin];
    v3_err[ibin] = V3_err[ibin]/V3[ibin]*v3[ibin];
  }

  TCanvas* c1D_short = new TCanvas("c1D_short","c1D_short",1000,800);
  c1D_short->Divide(5,4);
  TCanvas* c1D_long = new TCanvas("c1D_long","c1D_long",1000,800);
  c1D_long->Divide(4,3);

  for(int i=1;i<=nmultbin;i++) 
  {
    c1D_short->cd(i);
    hdata_short_zyam[i-1]->Add(hdata_long_zyam[i-1],-1);
    hdata_short_zyam[i-1]->SetAxisRange(-0.05,0.3,"Y");
    hdata_short_zyam[i-1]->Draw("PE");
    fitfunc_short[i-1]->SetLineColor(2);
//    fitfunc_short[i-1]->Draw("Lsame");
    c1D_long->cd(i);
//    hdata_long[i-1]->SetAxisRange(-0.02,0.8,"Y");
    hdata_long[i-1]->Draw("PE");
    fitfunc_long[i-1]->SetLineColor(2);
//    fitfunc_long[i-1]->Draw("Lsame");
    fitfunc_vn[i-1]->SetLineColor(2);
    fitfunc_vn[i-1]->Draw("Lsame");
  }
  
  TCanvas* c11D = new TCanvas("c11D","c11D",600,500);
  hdata_short_zyam[bin_periph]->Draw("PE");  
  hdata_short_zyam[bin_periph-1]->Draw("PESAME");

  double nmult_err[nmultbin]={0};
  TGraphErrors* gr_v0 = new TGraphErrors(nmultbin,nmult,V0,nmult_err,V0_err);
  TGraphErrors* gr_v2 = new TGraphErrors(nmultbin,nmult,V2,nmult_err,V2_err);
  TGraphErrors* gr_v3 = new TGraphErrors(nmultbin,nmult,V3,nmult_err,V3_err);
  TH2D* htmpv2 = new TH2D("tmpv2",";N;v_{2}",200,0,200,1000,-0.005,0.06*10);
  TH2D* htmpv3 = new TH2D("tmpv3",";N;v_{3}",200,0,200,1000,-0.01,0.05*2);
  
  TCanvas* c = new TCanvas("c","",550,500); 
  htmpv2->Draw();
  gr_v2->Draw("PESAME");

  TCanvas* c1 = new TCanvas("c1","",550,500);
  htmpv3->Draw();
  gr_v3->Draw("PESAME");

  TFile* ffout = new TFile(Form("results/%s_%sN120140_v2vsnmult.root",dataname.Data(),datatag[nmultbin-1].Data()),"recreate");
  gr_v0->Write("graph_v0");
  gr_v2->Write("graph_v2");
  gr_v3->Write("graph_v3");
  ffout->Close();

//  SaveCanvas(c1D_long,"pPb/corr","corr1D_pPb_allHFbins_longrange");
//  SaveCanvas(c1D_short,"pPb/corr","corr1D_pPb_allHFbins_shortrange");
return;
  ofstream fcout(Form("data/v2v3vsN_%s_pt1-2.txt",dataname.Data()));
  fcout << "N    v_2    v_2 err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< v2[i] << " " << v2_err[i] <<endl;
  fcout<<endl<<endl;
  fcout << "N    v_3    v_3 err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< v3[i] << " " << v3_err[i] <<endl;
  fcout<<endl<<endl;    

  ofstream fcout(Form("data/V2V3vsN_%s_pt1-2.txt",dataname.Data()));
  fcout << "N    V_2    V_2 err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< V2[i] << " " << V2_err[i] <<endl;
  fcout<<endl<<endl;
  fcout << "N    V_3    V_3 err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< V3[i] << " " << V3_err[i] <<endl;
  fcout<<endl<<endl;

  ofstream fcout(Form("data/NMultV2V3vsN_%s_pt1-2.txt",dataname.Data()));
  fcout << "N    V_2*NMult    V_2*NMult err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< V2[i]*V0[i] << " " << V2_err[i]*V0[i] <<endl;
  fcout<<endl<<endl;
  fcout << "N    V_3*NMult    V_3*NMult err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< V3[i]*V0[i] << " " << V3_err[i]*V0[i] <<endl;
  fcout<<endl<<endl;

  ofstream fcout(Form("data/JetYieldvsN_%s_pt1-2.txt",dataname.Data()));
  fcout << "N    yield    yield err." << endl;
  for(int i = 0; i<nmultbin-1; i++) fcout << nmult[i] << " "<< yield_short[i]-yield_long[i] << " " << sqrt(yielderror_short[i]*yielderror_short[i]+yielderror_long[i]*yielderror_long[i]) <<endl;
  fcout<<endl<<endl;

//  TFile* fout = new TFile("results/pPb_yieldvsn_pt1-2_hfall.root","recreate");
  TFile* fout = new TFile("results/PbPb_yieldvsn_pt1-2_all_hydjeteff.root","recreate");
//  TFile* fout = new TFile("results/hijing_yieldvsn_pt03-3.root","recreate");
  gr_yield_short_mult->SetName("gr_yield_short_mult");
  gr_yield_long_mult->SetName("gr_yield_long_mult");
  gr_yield_short_mult->Write();
  gr_yield_long_mult->Write();
  fout->Close();
}

void compare_v2vspt()
{

  TFile* f[10];
  f[0] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin0_nmax20.root");
  f[1] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin20_nmax40.root");
  f[2] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin40_nmax60.root");
  f[3] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin60_nmax85.root");
  f[4] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin85_nmax100.root");
  f[5] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin100_nmax120.root");
  f[6] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin120_nmax150.root");
  f[7] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin150_nmax1000.root");
  
  TFile* ff[10];
  ff[0] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin0_nmax20.root");
  ff[1] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin20_nmax40.root");
  ff[2] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin40_nmax60.root");
  ff[3] =  new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin60_nmax85.root");
  ff[4] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin85_nmax100.root");
  ff[5] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin100_nmax120.root");
  ff[6] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin120_nmax150.root");
  ff[7] =  new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin150_nmax1000.root");

/*
  TFile* f[10];
  f[0] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8_nmin0_nmax20.root");
  f[1] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8_nmin20_nmax40.root");
  f[2] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8_nmin40_nmax60.root");
  f[3] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8_nmin60_nmax80.root");
  f[4] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackPHFMPtWeight032v8_nmin80_nmax120.root");

  TFile* ff[10];
  ff[0] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v8_nmin0_nmax20.root");
  ff[1] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v8_nmin20_nmax40.root");
  ff[2] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v8_nmin40_nmax60.root");
  ff[3] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v8_nmin60_nmax80.root");
  ff[4] =  new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v8_nmin80_nmax120.root");
*/
  double nmult[]={8.3,27.4,46.7,67.0,91.1,107.5,128.4,156.7}; // pp data
//  double nmult[]={6.9,27.6,46.8,65.9,86.1}; // pythia Gen

  TGraphErrors* grp_v2[10];
  TGraphErrors* grm_v2[10];
  TGraphErrors* grp_v0[10];
  TGraphErrors* grm_v0[10];
  TGraphErrors* grp_v2v0[10];
  TGraphErrors* grm_v2v0[10];
  TGraphErrors* grpmult_v2v0[11];
  TGraphErrors* grmmult_v2v0[11];  

  for(int i=0;i<11;i++)
  {
    grpmult_v2v0[i] = new TGraphErrors(7);
    grmmult_v2v0[i] = new TGraphErrors(7);
  }

  for(int i=0;i<8;i++)
  {
    grp_v2[i] = (TGraphErrors*)f[i]->Get("graph_v2");
    grm_v2[i] = (TGraphErrors*)ff[i]->Get("graph_v2");
    grp_v0[i] = (TGraphErrors*)f[i]->Get("graph_v0");
    grm_v0[i] = (TGraphErrors*)ff[i]->Get("graph_v0");

    grp_v2v0[i] = GraphMultiplyv2(grp_v2[i],grp_v0[i]);
    grm_v2v0[i] = GraphMultiplyv2(grm_v2[i],grm_v0[i]);
    
    for(int j=0;j<grp_v2[i]->GetN();j++)
    {
      double x1,y1,x2,y2,x10,y10,x20,y20,y1err,y2err,y10err,y20err;
      grp_v2[i]->GetPoint(j,x1,y1);
      grm_v2[i]->GetPoint(j,x2,y2);
      y1err = grp_v2[i]->GetErrorY(j);
      y2err = grm_v2[i]->GetErrorY(j);

      grpmult_v2v0[j]->SetPoint(i,nmult[i],y1);
      grmmult_v2v0[j]->SetPoint(i,nmult[i],y2);
      grpmult_v2v0[j]->SetPointError(i,0,y1err);
      grmmult_v2v0[j]->SetPointError(i,0,y2err);
    }
  }
  TH2D* htmpv2 = new TH2D("tmpv2",";p_{T} (GeV/c);v_{2}",1000,0,10.0,1000,-0.01,0.15);
  TH2D* htmpv2mult = new TH2D("tmpv2mult",";N;v_{2}",1000,0,250.0,1000,-0.01,0.15);

  grp_v2v0[0]->SetMarkerColor(1);
  grp_v2v0[1]->SetMarkerColor(2);
  grp_v2v0[2]->SetMarkerColor(3);
  grp_v2v0[3]->SetMarkerColor(4);
  grp_v2v0[4]->SetMarkerColor(5);
//  grp_v2v0[5]->SetMarkerColor(6);
//  grp_v2v0[6]->SetMarkerColor(7);
  grm_v2v0[0]->SetMarkerColor(1);
  grm_v2v0[1]->SetMarkerColor(2);
  grm_v2v0[2]->SetMarkerColor(3);
  grm_v2v0[3]->SetMarkerColor(4);
  grm_v2v0[4]->SetMarkerColor(5);
//  grm_v2v0[5]->SetMarkerColor(6);
//  grm_v2v0[6]->SetMarkerColor(7);

  TCanvas* c = new TCanvas("c","c",500,500);
  htmpv2->Draw();
  grp_v2v0[0]->Draw("PSAME");
  grp_v2v0[1]->Draw("PSAME");
  grp_v2v0[2]->Draw("PSAME");
  grp_v2v0[3]->Draw("PSAME");
  grp_v2v0[4]->Draw("PSAME");
//  grp_v2v0[5]->Draw("PSAME");
//  grp_v2v0[6]->Draw("PSAME");

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  htmpv2->Draw();
  grm_v2v0[0]->Draw("PSAME");
  grm_v2v0[1]->Draw("PSAME");
  grm_v2v0[2]->Draw("PSAME");
  grm_v2v0[3]->Draw("PSAME");
  grm_v2v0[4]->Draw("PSAME");
//  grm_v2v0[5]->Draw("PSAME");
//  grm_v2v0[6]->Draw("PSAME");

  grpmult_v2v0[0]->SetMarkerColor(1);
  grpmult_v2v0[1]->SetMarkerColor(2);
  grpmult_v2v0[2]->SetMarkerColor(3);
  grpmult_v2v0[3]->SetMarkerColor(4);
  grpmult_v2v0[4]->SetMarkerColor(5);
  grpmult_v2v0[5]->SetMarkerColor(6);
  grpmult_v2v0[6]->SetMarkerColor(7);
  grmmult_v2v0[0]->SetMarkerColor(1);
  grmmult_v2v0[1]->SetMarkerColor(2);
  grmmult_v2v0[2]->SetMarkerColor(3);
  grmmult_v2v0[3]->SetMarkerColor(4);
  grmmult_v2v0[4]->SetMarkerColor(5);
  grmmult_v2v0[5]->SetMarkerColor(6);
  grmmult_v2v0[6]->SetMarkerColor(7);

  TCanvas* cc = new TCanvas("cc","cc",500,500);
  htmpv2mult->Draw();
  grpmult_v2v0[0]->Draw("PSAME");
  grpmult_v2v0[1]->Draw("PSAME");
  grpmult_v2v0[2]->Draw("PSAME");
  grpmult_v2v0[3]->Draw("PSAME");
  grpmult_v2v0[4]->Draw("PSAME");
  grpmult_v2v0[5]->Draw("PSAME");

  TCanvas* cc1 = new TCanvas("cc1","cc1",500,500);
  htmpv2mult->Draw();
  grmmult_v2v0[0]->Draw("PSAME");
  grmmult_v2v0[1]->Draw("PSAME");
  grmmult_v2v0[2]->Draw("PSAME");
  grmmult_v2v0[3]->Draw("PSAME");
  grmmult_v2v0[4]->Draw("PSAME");
  grmmult_v2v0[5]->Draw("PSAME");
}

void MCv2vspt()
{
//  TFile* f1 = new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v7_nmin20_nmax40.root");
//  TFile* f2 = new TFile("results/PythiaReco_Minbias_7TeV_TuneZ2_MinBias_42x_INCLGenTrackMHFPPtWeight032v7_nmin40_nmax60.root");
  
  TFile* f1 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v7Eta45_nmin150_nmax1000.root");
  TFile* f2 = new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v7Eta45_nmin20_nmax50.root");   

  TGraphErrors* gr1 = f1->Get("graph_v2");
  TGraphErrors* gr2 = f2->Get("graph_v2");
  TGraphErrors* gr10 = f1->Get("graph_v0");
  TGraphErrors* gr20 = f2->Get("graph_v0");

  TGraphErrors* grsub = new TGraphErrors(gr1->GetN());

  for(int i=0;i<gr1->GetN();i++)
  {
    double x1,y1,x2,y2,x10,y10,x20,y20,y1err,y2err,y10err,y20err;
    gr2->GetPoint(i,x2,y2);
    y2err = gr2->GetErrorY(i);
    gr10->GetPoint(i,x10,y10);
    gr20->GetPoint(i,x20,y20);
//    y10err = gr10->GetErrorY(i);
//    y20err = gr20->GetErrorY(i);

    double v2 = y2*y20/y10;
    grsub->SetPoint(i,x2,v2);
    grsub->SetPointError(i,0,y2err/y2*v2);
cout<<v2<<endl;
  }

  TH2D* htmpv2 = new TH2D("tmpv2",";p_{T} (GeV/c);v_{2}",1000,0,10.0,1000,-0.0,0.15);

  TCanvas* c = new TCanvas("c","c",500,500);
  htmpv2->Draw();
  gr1->Draw("PSAME");
  grsub->SetMarkerColor(2);
  grsub->Draw("PSAME");
}

void v2vsn()
{
/*
  TFile* f1 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackP033HFMPtWeight032Eta47_v2vsnmult.root");
  TFile* f2 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackM033HFPPtWeight032Eta47_v2vsnmult.root");
  TFile* f3 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032Eta47_v2vsnmult.root");
*/

  TFile* f1 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackP033HFMPtWeight032_v2vsnmult.root");
  TFile* f2 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackM033HFPPtWeight032_v2vsnmult.root");
  TFile* f3 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032_v2vsnmult.root");
 
  TGraphErrors* gr1 = f1->Get("graph_v2");
  TGraphErrors* gr2 = f2->Get("graph_v2");
  TGraphErrors* gr3 = f3->Get("graph_v2");
  TGraphErrors* gr10 = f1->Get("graph_v0");
  TGraphErrors* gr20 = f2->Get("graph_v0");
  TGraphErrors* gr30 = f3->Get("graph_v0");

  int binmb = 1;
  double scale = 1.5;
  double x1mb,y1mb,x2mb,y2mb,x10mb,y10mb,x20mb,y20mb,y1errmb,y2errmb,y10errmb,y20errmb,x0mb,y0mb,x00mb,y00mb,y0errmb,y00errmb;
  gr1->GetPoint(binmb,x1mb,y1mb);
  gr2->GetPoint(binmb,x2mb,y2mb);
  gr3->GetPoint(binmb,x0mb,y0mb);
  gr10->GetPoint(binmb,x10mb,y10mb);
  gr20->GetPoint(binmb,x20mb,y20mb);
  gr30->GetPoint(binmb,x00mb,y00mb);
  y1errmb = gr1->GetErrorY(binmb);
  y2errmb = gr2->GetErrorY(binmb);
  y0errmb = gr3->GetErrorY(binmb);
  y10errmb = gr10->GetErrorY(binmb);
  y20errmb = gr20->GetErrorY(binmb);
  y00errmb = gr30->GetErrorY(binmb);

  TGraphErrors* gr = new TGraphErrors(gr1->GetN());
  TGraphErrors* grsub = new TGraphErrors(gr1->GetN());

  for(int i=0;i<gr->GetN();i++)
  {
    double x1,y1,x2,y2,x10,y10,x20,y20,y1err,y2err,y10err,y20err,x0,y0,x00,y00,y0err,y00err;
    gr1->GetPoint(i,x1,y1);
    gr2->GetPoint(i,x2,y2);
    gr3->GetPoint(i,x0,y0);
    gr10->GetPoint(i,x10,y10);
    gr20->GetPoint(i,x20,y20);
    gr30->GetPoint(i,x00,y00);
    y1err = gr1->GetErrorY(i);
    y2err = gr2->GetErrorY(i);
    y0err = gr3->GetErrorY(i);
    y10err = gr10->GetErrorY(i);
    y20err = gr20->GetErrorY(i);
    y00err = gr30->GetErrorY(i);

//    double v2 = y2/sqrt(y0);
    double v2 = sqrt(y1*y2/y0);
    gr->SetPoint(i,x1,v2);
    gr->SetPointError(i,0,v2*0.5*sqrt((y0err/y0)*(y0err/y0)+(y1err/y1)*(y1err/y1)+(y2err/y2)*(y2err/y2)));
//    gr->SetPointError(i,0,v2*0.5*sqrt((y1err/y1)*(y1err/y1)+(y2err/y2)*(y2err/y2)));
    
    double v2sub = sqrt((y1-scale*y1mb*y10mb/y10)*(y2-scale*y2mb*y20mb/y20)/(y0-scale*y0mb*y00mb/y00));
    grsub->SetPoint(i,x1,v2sub);
    grsub->SetPointError(i,0,v2sub*0.5*sqrt((y0err/y0)*(y0err/y0)+(y1err/y1)*(y1err/y1)+(y2err/y2)*(y2err/y2)));
  }

  TH2D* htmpv2 = new TH2D("tmpv2",";N;v_{2}",200,0,200,1000,-0.0,0.25);

  TCanvas* c = new TCanvas("c","c",500,500);
  htmpv2->Draw();
  gr->Draw("PSAME");
  grsub->SetMarkerColor(2);
  grsub->Draw("PSAME");
}
void v2vspt()
{
  TFile* f1 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin100_nmax120.root");
  TFile* f2 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin100_nmax120.root");
  TFile* f3 = new TFile("results/PPData_Minbias_7TeV_TRIGGER_All_Apr21ReReco_INCLHFPHFMPtWeight032_nmin100_nmax120.root");

  TGraphErrors* gr1 = f1->Get("graph_v2");
  TGraphErrors* gr2 = f2->Get("graph_v2");
  TGraphErrors* gr3 = f3->Get("graph_v2");
  TGraphErrors* gr10 = f1->Get("graph_v0");
  TGraphErrors* gr20 = f2->Get("graph_v0");
  TGraphErrors* gr30 = f3->Get("graph_v0");

  TFile* f1mb = new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackMHFPPtWeight032v8_nmin10_nmax35.root");
  TFile* f2mb = new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLTrackPHFMPtWeight032v8_nmin10_nmax35.root");
  TFile* f3mb = new TFile("results/PPData_Minbias_7TeV_MINBIAS_All_Apr21ReReco_INCLHFPHFMPtWeight032_nmin10_nmax35.root");

  TGraphErrors* gr1mb = f1mb->Get("graph_v2");
  TGraphErrors* gr2mb = f2mb->Get("graph_v2");
  TGraphErrors* gr3mb = f3mb->Get("graph_v2");
  TGraphErrors* gr10mb = f1mb->Get("graph_v0");
  TGraphErrors* gr20mb = f2mb->Get("graph_v0");
  TGraphErrors* gr30mb = f3mb->Get("graph_v0");

  double x0,y0,x00,y00,y0err,y00err;
  gr3->GetPoint(0,x0,y0);
  gr30->GetPoint(0,x00,y00);
  y0err = gr3->GetErrorY(0);
  y00err = gr30->GetErrorY(0);

  double x0mb,y0mb,x00mb,y00mb,y0mberr,y00mberr;
  gr3mb->GetPoint(0,x0mb,y0mb);
  gr30mb->GetPoint(0,x00mb,y00mb);
  y0mberr = gr3mb->GetErrorY(0);
  y00mberr = gr30mb->GetErrorY(0);

  TGraphErrors* gr = new TGraphErrors(gr1->GetN());
  TGraphErrors* grsub = new TGraphErrors(gr1->GetN());
  
  for(int i=0;i<gr->GetN();i++)
  {
    double x1,y1,x2,y2,x10,y10,x20,y20,y1err,y2err,y10err,y20err;
    gr1->GetPoint(i,x1,y1);
    gr2->GetPoint(i,x2,y2);
    gr10->GetPoint(i,x10,y10);
    gr20->GetPoint(i,x20,y20);
    y1err = gr1->GetErrorY(i);
    y2err = gr2->GetErrorY(i);
    y10err = gr10->GetErrorY(i);
    y20err = gr20->GetErrorY(i);

    double v2 = sqrt(y1*y2/y0);
    gr->SetPoint(i,x1,v2);
    gr->SetPointError(i,0,v2*0.5*sqrt((y0err/y0)*(y0err/y0)+(y1err/y1)*(y1err/y1)+(y2err/y2)*(y2err/y2)));

//cout<<y1err/y1<<" "<<y2err/y2<<" "<<y10err/y10<<" "<<y20err/y20<<endl;
    double x1mb,y1mb,x2mb,y2mb,x10mb,y10mb,x20mb,y20mb,y1mberr,y2mberr,y10mberr,y20mberr;
    gr1mb->GetPoint(i,x1mb,y1mb);
    gr2mb->GetPoint(i,x2mb,y2mb);
    gr10mb->GetPoint(i,x10mb,y10mb);
    gr20mb->GetPoint(i,x20mb,y20mb);
    y1mberr = gr1mb->GetErrorY(i);
    y2mberr = gr2mb->GetErrorY(i);
    y10mberr = gr10mb->GetErrorY(i);
    y20mberr = gr20mb->GetErrorY(i);

    double v2sub = sqrt((y1-y10mb*y1mb/y10)*(y2-y20mb*y2mb/y20)/(y0-y0mb*y00mb/y00));
    grsub->SetPoint(i,x1,v2sub);
    grsub->SetPointError(i,0,v2sub*0.5*sqrt((y0err/y0)*(y0err/y0)+(y00err/y00)*(y00err/y00)+(y1err/y1)*(y1err/y1)+(y10err/y10)*(y10err/y10)+(y2err/y2)*(y2err/y2)+(y20err/y20)*(y20err/y20)));
  }
 
  TH2D* htmpv2 = new TH2D("tmpv2",";p_{T} (GeV/c);v_{2}",1000,0,10.0,1000,-0.0,0.15);

  TCanvas* c = new TCanvas("c","c",500,500);
  htmpv2->Draw();
  gr->Draw("PSAME");
  grsub->SetMarkerColor(2);
  grsub->Draw("PSAME");
//  grmb->Draw("PSAME");
}
