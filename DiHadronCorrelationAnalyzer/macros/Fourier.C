#include <fstream>

void Fourier()
{
/*
  TString dataname = "HIData_Minbias_2760GeV";
  TString datatag = "PPRereco_INCLEFF1v4";
*/
  TString dataname = "PAData_Minbias_5TeV";
  TString datatag = "HM_Express_INCLEFF1v5";
//  TString datatag = "MB_INCLEFF1v4";

  int nmin=220;
  int nmax=260;
  const int ntrgptbins=11;
  const int nassptbins=11;

  TH1D ***hdata = compare_flow_pp_phi_multipad_pt<TH1D**>(dataname.Data(),datatag.Data(),nmin,nmax,-1,-1,2.0,4.0,ntrgptbins,nassptbins,-2.4,2.4,-2.4,2.4);
  TH1D* hdata_zyam[ntrgptbins][nassptbins];
  TF1* fitfunc_vn[ntrgptbins][nassptbins];

  for(int itrg = 0; itrg<ntrgptbins; itrg++)
  {
    for(int jass = 0; jass<nassptbins; jass++)
    {
      fitfunc_vn[itrg][jass] = FitVnFunc(hdata[itrg][jass]);
    }
  }

  int bin_ref=1;
  double V2_ref = fitfunc_vn[bin_ref][bin_ref]->GetParameter(2);
  double V3_ref = fitfunc_vn[bin_ref][bin_ref]->GetParameter(3);
  double V2_ref_err = fitfunc_vn[bin_ref][bin_ref]->GetParError(2);
  double V3_ref_err = fitfunc_vn[bin_ref][bin_ref]->GetParError(3);

  double v2[ntrgptbins];
  double v3[ntrgptbins];
  double v2_err[ntrgptbins];
  double v3_err[ntrgptbins];
  double pttrg[ntrgptbins] = {0.0};
  double pttrg_err[ntrgptbins] = {0.0};

  TH1D* hpt[ntrgptbins];
  for(int itrg = 0; itrg<ntrgptbins; itrg++)
  {
    TString filename = Form("/net/hisrv0001/home/davidlw/scratch1/DiHadronCorrelations/outputs_312/%s/merged/%s_%s_nmin%d_nmax%d_etatrg%.1f-%.1f_etaass%.1f-%.1f_centmin%d_centmax%d.root",dataname.Data(),dataname.Data(),datatag.Data(),nmin,nmax,-1,-1,-2.4,2.4,-2.4,2.4);    
    hpt[itrg] = (TH1D*)GetHist(filename,Form("ptcorr_signal_trg_%d",itrg));
    pttrg[itrg] = hpt[itrg]->GetMean();
  }

  for(int itrg = 0; itrg<ntrgptbins; itrg++)
  {
    v2[itrg] = fitfunc_vn[itrg][bin_ref]->GetParameter(2)/sqrt(fabs(V2_ref));
    v3[itrg] = fitfunc_vn[itrg][bin_ref]->GetParameter(3)/sqrt(fabs(V3_ref));
    v2_err[itrg] = fitfunc_vn[itrg][bin_ref]->GetParError(2)/fitfunc_vn[itrg][bin_ref]->GetParameter(2)*v2[itrg];
    v3_err[itrg] = fitfunc_vn[itrg][bin_ref]->GetParError(3)/fitfunc_vn[itrg][bin_ref]->GetParameter(3)*v3[itrg];
  }

  ofstream fcout(Form("data/%s_%s_nmin%d_nmax%d_centmin%d_centmax%d_ptass%d.txt",dataname.Data(),datatag.Data(),nmin,nmax,-1,-1,bin_ref));
  fcout << "pT    v_2    v_2 err." << endl;
  for(int itrg = 0; itrg<ntrgptbins; itrg++) fcout << pttrg[itrg] << " "<< v2[itrg] << " " << v2_err[itrg] <<endl;
  fcout<<endl<<endl;
  fcout << "pT    v_3    v_3 err." << endl;
  for(int itrg = 0; itrg<ntrgptbins; itrg++) fcout << pttrg[itrg] << " "<< v3[itrg] << " " << v3_err[itrg] <<endl;
  fcout<<endl<<endl;

  TGraphErrors* g_v2 = new TGraphErrors(ntrgptbins,pttrg,v2,pttrg_err,v2_err);
  TGraphErrors* g_v3 = new TGraphErrors(ntrgptbins,pttrg,v3,pttrg_err,v3_err);
  g_v3->SetMarkerStyle(24);

  TCanvas* cv = new TCanvas("cv","cv",600,550);
  TH2D* hdummy = new TH2D("dummy",";p_{T} (GeV/c);v_{2}",100,0,12.0,100,0,0.25);
  hdummy->Draw();
  g_v2->Draw("PESAME");
  g_v3->Draw("PESAME");
}
