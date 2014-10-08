#ifndef setupDiHadronCorrelation_
#define setupDiHadronCorrelation_

#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TList.h>
#include <TString.h>
#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FlowCorrAna/DiHadronCorrelationAnalyzer/interface/CutParameters.h"
#include "RunAna.C"
#endif

int gCheckSign;
int gFfrom;
int gFto;
int gNCentBins;
int gCentRunnum;
int gNEtaBins;
int gNPhiBins;
int gBkgFactor;
int gEPIndex;
CutParameters gCut;
TString gFilesetname;
TString gOutputDir;
TString gCentfilename;
TString gCenttablename;
TString gEffhistname;
TString gTrghistname;
TString gPileupdistfunchistname;
TString gEventtype;
TString gTreename;
TString gTag;
/*
void setup(int ffrom=0, int fto=1,
      	   double zvtxmin=-1.0, double zvtxmax=1.0,
           int nmin=-1, int nmax=-1,
           int centmin=-1, int centmax=-1);
*/
void setup(int ffrom, int fto,
           double zvtxmin, double zvtxmax,
           int nmin, int nmax,
           int centmin, int centmax)
{ 
  if(gSystem->Getenv("FILESETNAME")) 
    gFilesetname = gSystem->Getenv("FILESETNAME");

  if(gSystem->Getenv("OUTPUTDIR"))
    gOutputDir = gSystem->Getenv("OUTPUTDIR");

  if(gSystem->Getenv("CENTFILENAME"))
    gCentfilename = gSystem->Getenv("CENTFILENAME");
  
  if(gSystem->Getenv("CENTTABLENAME"))
    gCenttablename = gSystem->Getenv("CENTTABLENAME");

  if(gSystem->Getenv("EFFHISTNAME"))
    gEffhistname = gSystem->Getenv("EFFHISTNAME");  
  
  if(gSystem->Getenv("TRGHISTNAME"))
    gTrghistname = gSystem->Getenv("TRGHISTNAME");

  if(gSystem->Getenv("PILEUPDISTFUNCHISTNAME"))
    gPileupdistfunchistname = gSystem->Getenv("PILEUPDISTFUNCHISTNAME");

  if(gSystem->Getenv("EVENTTYPE"))
    gEventtype = gSystem->Getenv("EVENTTYPE");
  
  if(gSystem->Getenv("TAG")) 
    gTag = gSystem->Getenv("TAG");

  TString stretatrgmin = "-2.0";
  if(gSystem->Getenv("ETATRGMIN")) stretatrgmin = gSystem->Getenv("ETATRGMIN");
  gCut.etatrgmin = atof(stretatrgmin);
  
  TString stretatrgmax = "2.0";
  if(gSystem->Getenv("ETATRGMAX")) stretatrgmax = gSystem->Getenv("ETATRGMAX");
  gCut.etatrgmax = atof(stretatrgmax);

  TString stretaassmin = "-2.0";
  if(gSystem->Getenv("ETAASSMIN")) stretaassmin = gSystem->Getenv("ETAASSMIN");
  gCut.etaassmin = atof(stretaassmin);

  TString stretaassmax = "2.0";
  if(gSystem->Getenv("ETAASSMAX")) stretaassmax = gSystem->Getenv("ETAASSMAX");
  gCut.etaassmax = atof(stretaassmax);

  TString stretamultmin = "-2.0";
  if(gSystem->Getenv("ETAMULTMIN")) stretamultmin = gSystem->Getenv("ETAMULTMIN");
  gCut.etamultmin = atof(stretamultmin);

  TString stretamultmax = "2.0";
  if(gSystem->Getenv("ETAMULTMAX")) stretamultmax = gSystem->Getenv("ETAMULTMAX");
  gCut.etamultmax = atof(stretamultmax);

  TString strptmultmin = "0.0";
  if(gSystem->Getenv("PTMULTMIN")) strptmultmin = gSystem->Getenv("PTMULTMIN");
  gCut.ptmultmin = atof(strptmultmin);

  TString strptmultmax = "10000.0";
  if(gSystem->Getenv("PTMULTMAX")) strptmultmax = gSystem->Getenv("PTMULTMAX");
  gCut.ptmultmax = atof(strptmultmax);
  
  TString strncentbins = "40";
  if(gSystem->Getenv("NCENTBINS")) strncentbins = gSystem->Getenv("NCENTBINS");
  gNCentBins = atoi(strncentbins);

  TString strcentrunnum = "181502";
  if(gSystem->Getenv("CENTRUNNUM")) strcentrunnum = gSystem->Getenv("CENTRUNNUM");
  gCentRunnum = atoi(strcentrunnum);

  TString strnetabins = "40";
  if(gSystem->Getenv("NETABINS")) strnetabins = gSystem->Getenv("NETABINS");
  gNEtaBins = atoi(strnetabins);
 
  TString strnphibins = "32";
  if(gSystem->Getenv("NPHIBINS")) strnphibins = gSystem->Getenv("NPHIBINS");
  gNPhiBins = atoi(strnphibins);

  TString strbkgfactor = "10";
  if(gSystem->Getenv("BKGFACTOR")) strbkgfactor = gSystem->Getenv("BKGFACTOR");
  gBkgFactor = atoi(strbkgfactor);

  TString strmasstrg = "0.140";
  if(gSystem->Getenv("MASSTRG")) strmasstrg = gSystem->Getenv("MASSTRG");
  gCut.mass_trg = atof(strmasstrg);

  TString strmassass = "0.140";
  if(gSystem->Getenv("MASSASS")) strmassass = gSystem->Getenv("MASSASS");
  gCut.mass_ass = atof(strmassass);

  TString strpdgidtrg = "-999999";
  if(gSystem->Getenv("PDGIDTRG")) strpdgidtrg = gSystem->Getenv("PDGIDTRG");
  gCut.genpdgId_trg = atoi(strpdgidtrg);

  TString strpdgidass = "-999999";
  if(gSystem->Getenv("PDGIDASS")) strpdgidass = gSystem->Getenv("PDGIDASS");
  gCut.genpdgId_ass = atoi(strpdgidass);

  TString strisstabletrg = "1";
  if(gSystem->Getenv("ISSTABLETRG")) strisstabletrg = gSystem->Getenv("ISSTABLETRG");
  gCut.isstable_trg = atoi(strisstabletrg);

  TString strisstableass = "1";
  if(gSystem->Getenv("ISSTABLEASS")) strisstableass = gSystem->Getenv("ISSTABLEASS");
  gCut.isstable_ass = atoi(strisstableass);

  TString strischargetrg = "1";
  if(gSystem->Getenv("ISCHARGETRG")) strischargetrg = gSystem->Getenv("ISCHARGETRG");
  gCut.ischarge_trg = atoi(strischargetrg);

  TString strischargeass = "1";
  if(gSystem->Getenv("ISCHARGEASS")) strischargeass = gSystem->Getenv("ISCHARGEASS");
  gCut.ischarge_ass = atoi(strischargeass);

  TString strisgenmult = "0";
  if(gSystem->Getenv("ISGENMULT")) strisgenmult = gSystem->Getenv("ISGENMULT");
  gCut.IsGenMult = atoi(strisgenmult);

  TString strisgenvtx = "0";
  if(gSystem->Getenv("ISGENVTX")) strisgenvtx = gSystem->Getenv("ISGENVTX");
  gCut.IsGenVtx = atoi(strisgenvtx);

  TString strisvtxsel = "0";
  if(gSystem->Getenv("ISVTXSEL")) strisvtxsel = gSystem->Getenv("ISVTXSEL");
  gCut.IsVtxSel = atoi(strisvtxsel);

  TString striseventengineer = "0";
  if(gSystem->Getenv("ISEVENTENGINEER")) striseventengineer = gSystem->Getenv("ISEVENTENGINEER");
  gCut.IsEventEngineer = atoi(striseventengineer);

  TString striscorr = "1";
  if(gSystem->Getenv("ISCORR")) striscorr = gSystem->Getenv("ISCORR");
  gCut.IsCorr = atoi(striscorr);

  TString strisinvmass = "0";
  if(gSystem->Getenv("ISINVMASS")) strisinvmass = gSystem->Getenv("ISINVMASS");
  gCut.IsInvMass = atoi(strisinvmass);

  TString strishi = "0";
  if(gSystem->Getenv("ISHI")) strishi = gSystem->Getenv("ISHI");
  gCut.IsHI = atoi(strishi);

  TString strisdoeventshape = "0";
  if(gSystem->Getenv("ISDOEVENTSHAPE")) strisdoeventshape = gSystem->Getenv("ISDOEVENTSHAPE");
  gCut.IsDoEventShape = atoi(strisdoeventshape);

  TString stristrackntuple = "0";
  if(gSystem->Getenv("ISTRACKNTUPLE")) stristrackntuple = gSystem->Getenv("ISTRACKNTUPLE");
  gCut.IsTrackNtuple = atoi(stristrackntuple);

  TString strishigeninfontuple = "0";
  if(gSystem->Getenv("ISHIGENINFONTUPLE")) strishigeninfontuple = gSystem->Getenv("ISHIGENINFONTUPLE");
  gCut.IsHIGenInfoNtuple = atoi(strishigeninfontuple);

  TString strisleadtrack = "0";
  if(gSystem->Getenv("ISLEADTRACK")) strisleadtrack = gSystem->Getenv("ISLEADTRACK");
  gCut.IsLeadTrack = atoi(strisleadtrack);

  TString strissubjet = "0";
  if(gSystem->Getenv("ISSUBJET")) strissubjet = gSystem->Getenv("ISSUBJET");
  gCut.IsSubJet = atoi(strissubjet);

  TString strismonojet = "0";
  if(gSystem->Getenv("ISMONOJET")) strismonojet = gSystem->Getenv("ISMONOJET");
  gCut.IsMonoJet = atoi(strismonojet);

  TString strisfullmatrix = "0";
  if(gSystem->Getenv("ISFULLMATRIX")) strisfullmatrix = gSystem->Getenv("ISFULLMATRIX");
  gCut.IsFullMatrix = atoi(strisfullmatrix);

  TString strchecksign = "-1";
  if(gSystem->Getenv("CHECKSIGN")) strchecksign = gSystem->Getenv("CHECKSIGN");
  gCheckSign = atoi(strchecksign);

  TString strissymmetrize = "0";
  if(gSystem->Getenv("ISSYMMETRIZE")) strissymmetrize = gSystem->Getenv("ISSYMMETRIZE");
  gCut.IsSymmetrize = atoi(strissymmetrize);

  TString strisptweighttrg = "0";
  if(gSystem->Getenv("ISPTWEIGHTTRG")) strisptweighttrg = gSystem->Getenv("ISPTWEIGHTTRG");
  gCut.IsPtWeightTrg = atoi(strisptweighttrg);

  TString strisptweightass = "0";
  if(gSystem->Getenv("ISPTWEIGHTASS")) strisptweightass = gSystem->Getenv("ISPTWEIGHTASS");
  gCut.IsPtWeightAss = atoi(strisptweightass);

  TString strisharmonics = "0";
  if(gSystem->Getenv("ISHARMONICS")) strisharmonics = gSystem->Getenv("ISHARMONICS");
  gCut.IsHarmonics = atoi(strisharmonics);

  TString stristrkquality = "0";
  if(gSystem->Getenv("ISTRKQUALITY")) stristrkquality = gSystem->Getenv("ISTRKQUALITY");
  gCut.IsTrkQuality = atoi(stristrkquality);

  TString strepindex = "-1";
  if(gSystem->Getenv("EPINDEX")) strepindex = gSystem->Getenv("EPINDEX");
  gEPIndex = atof(strepindex);

  string strtriggernamesarray = "";
  if(gSystem->Getenv("TRIGGERNAMES")) strtriggernamesarray = gSystem->Getenv("TRIGGERNAMES");
  stringstream ss_triggernames(strtriggernamesarray,ios_base::out | ios_base::in);
  char buffer_triggername[1024];
  while(ss_triggernames>>buffer_triggername)
  {
    TString string_triggername(buffer_triggername); 
    gCut.triggernames.push_back(string_triggername);
  }

  string strvetonamesarray = "";
  if(gSystem->Getenv("VETONAMES")) strvetonamesarray = gSystem->Getenv("VETONAMES");
  stringstream ss_vetonames(strvetonamesarray,ios_base::out | ios_base::in);
  char buffer_vetoname[1024];
  while(ss_vetonames>>buffer_vetoname)
  {
    TString string_vetoname(buffer_vetoname);
    gCut.vetonames.push_back(string_vetoname);
  }

  string strtrgtrackcollection = "generalTracks";
  if(gSystem->Getenv("TRGTRACKCOLLECTION")) strtrgtrackcollection = gSystem->Getenv("TRGTRACKCOLLECTION");
  gCut.trgtrackCollection = strtrgtrackcollection;

  string strasstrackcollection = "generalTracks";
  if(gSystem->Getenv("ASSTRACKCOLLECTION")) strasstrackcollection = gSystem->Getenv("ASSTRACKCOLLECTION");
  gCut.asstrackCollection = strasstrackcollection;

  string strvertexcollection = "offlinePrimaryVertices";
  if(gSystem->Getenv("VERTEXCOLLECTION")) strvertexcollection = gSystem->Getenv("VERTEXCOLLECTION");
  gCut.vertexCollection = strvertexcollection;

  string strjetcollection = "selectedPatJets";
  if(gSystem->Getenv("JETCOLLECTION")) strjetcollection = gSystem->Getenv("JETCOLLECTION");
  gCut.jetCollection = strjetcollection;

  string strgenparticlecollection = "genParticles";
  if(gSystem->Getenv("GENPARTICLECOLLECTION")) strgenparticlecollection = gSystem->Getenv("GENPARTICLECOLLECTION");
  gCut.genParticleCollection = strgenparticlecollection;

  string strgenjetcollection = "ak5GenJets";
  if(gSystem->Getenv("GENJETCOLLECTION")) strgenjetcollection = gSystem->Getenv("GENJETCOLLECTION");
  gCut.genJetCollection = strgenjetcollection;

  string strv0candidatecollection = "generalV0Candidates";
  if(gSystem->Getenv("V0CANDIDATECOLLECTION")) strv0candidatecollection = gSystem->Getenv("V0CANDIDATECOLLECTION");
  gCut.v0CandidateCollection = strv0candidatecollection;

  string strsuperclustercollection = "correctedHybridSuperClusters";
  if(gSystem->Getenv("SUPERCLUSTERCOLLECTION")) strsuperclustercollection = gSystem->Getenv("SUPERCLUSTERCOLLECTION");
  gCut.superClusterCollection = strsuperclustercollection;

  string strsuperclusterfwdcollection = "correctedMulti5x5SuperClusterWithPreshower";
  if(gSystem->Getenv("SUPERCLUSTERFWDCOLLECTION")) strsuperclusterfwdcollection = gSystem->Getenv("SUPERCLUSTERFWDCOLLECTION");
  gCut.superClusterFwdCollection = strsuperclusterfwdcollection;

  string strcentralitycollection = "hiCentrality";
  if(gSystem->Getenv("CENTRALITYCOLLECTION")) strcentralitycollection = gSystem->Getenv("CENTRALITYCOLLECTION");
  gCut.centralityCollection = strcentralitycollection;

  TString strpthard1 = "1.5";
  if(gSystem->Getenv("PTHARD1")) strpthard1 = gSystem->Getenv("PTHARD1");
  gCut.pthard1 = atof(strpthard1);

  TString strpthard2 = "0.7";
  if(gSystem->Getenv("PTHARD2")) strpthard2 = gSystem->Getenv("PTHARD2");
  gCut.pthard2 = atof(strpthard2);

  TString strdelta = "0.5";
  if(gSystem->Getenv("DELTA")) strdelta = gSystem->Getenv("DELTA");
  gCut.delta = atof(strdelta);

  TString strrunmin = "1";
  if(gSystem->Getenv("RUNMIN")) strrunmin = gSystem->Getenv("RUNMIN");
  gCut.runmin = atoi(strrunmin);
  
  TString strrunmax = "1000000";
  if(gSystem->Getenv("RUNMAX")) strrunmax = gSystem->Getenv("RUNMAX");
  gCut.runmax = atoi(strrunmax);

  TString strxvtxcenter = "0.0";
  if(gSystem->Getenv("XVTXCENTER")) strxvtxcenter = gSystem->Getenv("XVTXCENTER");
  gCut.xvtxcenter = atof(strxvtxcenter);

  TString stryvtxcenter = "0.0";
  if(gSystem->Getenv("YVTXCENTER")) stryvtxcenter = gSystem->Getenv("YVTXCENTER");
  gCut.yvtxcenter = atof(stryvtxcenter);

  TString strzvtxcenter = "0.0";
  if(gSystem->Getenv("ZVTXCENTER")) strzvtxcenter = gSystem->Getenv("ZVTXCENTER");
  gCut.zvtxcenter = atof(strzvtxcenter);

  TString strzvtxbin = "0.0";
  if(gSystem->Getenv("ZVTXBIN")) strzvtxbin = gSystem->Getenv("ZVTXBIN");
  gCut.zvtxbin = atof(strzvtxbin)/10;

  TString strvtxsepmin = "0.0";
  if(gSystem->Getenv("VTXSEPMIN")) strvtxsepmin = gSystem->Getenv("VTXSEPMIN");
  gCut.vtxsepmin = atof(strvtxsepmin);

  TString strnvtxmax = "100";
  if(gSystem->Getenv("NVTXMAX")) strnvtxmax = gSystem->Getenv("NVTXMAX");
  gCut.nvtxmax = atoi(strnvtxmax);

  TString strrhomin = "-10000.0";
  if(gSystem->Getenv("RHOMIN")) strrhomin = gSystem->Getenv("RHOMIN");
  gCut.rhomin = atof(strrhomin);
  
  TString strrhomax = "10000.0";
  if(gSystem->Getenv("RHOMAX")) strrhomax = gSystem->Getenv("RHOMAX");
  gCut.rhomax = atof(strrhomax);
  
  TString strnjetsmin = "-1";
  if(gSystem->Getenv("NJETSMIN")) strnjetsmin = gSystem->Getenv("NJETSMIN");
  gCut.njetsmin = atof(strnjetsmin);

  TString strnjetsmax = "-1";
  if(gSystem->Getenv("NJETSMAX")) strnjetsmax = gSystem->Getenv("NJETSMAX");
  gCut.njetsmax = atof(strnjetsmax);

  TString strjetetmin = "0.0";
  if(gSystem->Getenv("JETETMIN")) strjetetmin = gSystem->Getenv("JETETMIN");
  gCut.jetetmin = atof(strjetetmin);

  TString strjetetmax = "10000.0";
  if(gSystem->Getenv("JETETMAX")) strjetetmax = gSystem->Getenv("JETETMAX");
  gCut.jetetmax = atof(strjetetmax);

  TString strjetetfraction = "0.05";
  if(gSystem->Getenv("JETETFRACTION")) strjetetfraction = gSystem->Getenv("JETETFRACTION");
  gCut.jetetfraction = atof(strjetetfraction);

  TString strleadjetetmin = "120";
  if(gSystem->Getenv("LEADJETETMIN")) strleadjetetmin = gSystem->Getenv("LEADJETETMIN");
  gCut.leadjetetmin = atof(strleadjetetmin);

  TString strsubjetetmin = "50";
  if(gSystem->Getenv("SUBJETETMIN")) strsubjetetmin = gSystem->Getenv("SUBJETETMIN");
  gCut.subjetetmin = atof(strsubjetetmin);
 
  string strpttrgminarray = "20";
  if(gSystem->Getenv("PTTRGMINARRAY")) strpttrgminarray = gSystem->Getenv("PTTRGMINARRAY");
  stringstream ss_pttrgmin(strpttrgminarray,ios_base::out | ios_base::in);
  double buffer_pttrgmin=0;
  while(ss_pttrgmin>>buffer_pttrgmin) gCut.pttrgmin.push_back(buffer_pttrgmin);

  string strpttrgmaxarray = "40";
  if(gSystem->Getenv("PTTRGMAXARRAY")) strpttrgmaxarray = gSystem->Getenv("PTTRGMAXARRAY");
  stringstream ss_pttrgmax(strpttrgmaxarray,ios_base::out | ios_base::in);
  double buffer_pttrgmax=0;
  while(ss_pttrgmax>>buffer_pttrgmax) gCut.pttrgmax.push_back(buffer_pttrgmax);

  string strptassminarray = "6 8 10 15";
  if(gSystem->Getenv("PTASSMINARRAY")) strptassminarray = gSystem->Getenv("PTASSMINARRAY");
  stringstream ss_ptassmin(strptassminarray,ios_base::out | ios_base::in);
  double buffer_ptassmin=0;
  while(ss_ptassmin>>buffer_ptassmin) gCut.ptassmin.push_back(buffer_ptassmin);

  string strptassmaxarray = "8 10 15 20";
  if(gSystem->Getenv("PTASSMAXARRAY")) strptassmaxarray = gSystem->Getenv("PTASSMAXARRAY");
  stringstream ss_ptassmax(strptassmaxarray,ios_base::out | ios_base::in);;
  double buffer_ptassmax=0;
  while(ss_ptassmax >> buffer_ptassmax) gCut.ptassmax.push_back(buffer_ptassmax);

  gFfrom=ffrom;
  gFto=fto; 
  gCut.nmin=nmin;
  gCut.nmax=nmax;
  gCut.centmin=centmin;
  gCut.centmax=centmax;
  gCut.zvtxmin=zvtxmin/10.;
  gCut.zvtxmax=zvtxmax/10.;
}

#endif //setupDiHadronCorrelation_
