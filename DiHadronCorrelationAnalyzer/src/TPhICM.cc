#include "../interface/TPhICM.h"
#include <iostream>

#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TF1.h"
#include "TF2.h"
#include "TList.h"
#include "TFile.h"
#include "TTree.h"
#include "TDatime.h"
#include "TString.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TClonesArray.h"

using namespace std;

TPhICM::TPhICM(Int_t nevents, Int_t seed, Int_t multiplicity, Int_t modelnum):
   fNevents(nevents),
   fMultiplicity(multiplicity),
   fInitMultiplicity(multiplicity),
   fNuclR(8.0),
   fB(0.0), 
   fEccA(0),
   fEccB(0),
   fAbsPar(0),
   fBetaRadial(0.0),
   fErrorflag(0),
   fConserveMom(0),
   fOutTree(0),
   fOutTreef(0),
   EventKind(1),
   MonteCarlo(0),
   fEtaShift(0),
   fEtaBroaden(0)
{
   gRandom=new TRandom3(seed);
   FixModelParameters(modelnum);
}

void TPhICM::FixModelParameters(Int_t modelnum) {

   for(Int_t model=0;model<6;model++) {
      fClusProb[model] = 0.0;
      fClusNchild[model] = 0;
      fClusM[model] = 0.0;
      fClusPzDist[model] = 0;
//      fClusPtDist[model] = new TF1("ptdist","x*exp(-x/0.45)",0,10);
      fClusPtDist[model] = 0; 
      fClusMf[model]=0;
      fClusEtaBroad[model]=0;
   }
   if(modelnum==0) {
      fConserveMom=0;      
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","8.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
   }
   else if(modelnum==1) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","8.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
   }
   else if(modelnum==2) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200); 

      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
   }
   else if(modelnum==3) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.0/cosh(x/1.3)+3.0/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200); 
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
   }
   else if(modelnum==4) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200); 
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
   }
   else if(modelnum==5) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","4.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
      fEtaBroaden = 0.7;
   }
   else if(modelnum==6) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.8/cosh(x/1.3)+2.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 0.7;
   }
   else if(modelnum==7) { // CuCu 20%-30%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.5/cosh(x/1.3)+1.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = 0.7;
   }
   else if(modelnum==8) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.7/cosh(x/1.3)+1.8/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = 0.7;
   }
   else if(modelnum==9) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10); 
      fClusPzDist[0] = new TF1("pzdist","4.0/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
      fEtaBroaden = 0.9;
   }
   else if(modelnum==10) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.2/cosh(x/1.3)+2.2/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 0.9;
   }  
   else if(modelnum==11) {   
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.1/cosh(x/1.3)+1.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = 0.9;
   }  
   else if(modelnum==12) {      
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.4/cosh(x/1.3)+1.8/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = 0.9;
   }  
   else if(modelnum==13) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.0/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
      fEtaBroaden = 0.5;
   }
   else if(modelnum==14) { // CuCu 0%-10%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 0.5;
   }  
   else if(modelnum==15) {   
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.0/cosh(x/1.3)+1.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = 0.5;
   }  
   else if(modelnum==16) {      
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.0/cosh(x/1.3)+2.0/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = 0.5;
   } 
   else if(modelnum==17) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","6.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
      fEtaBroaden = 0.3;
   }
   else if(modelnum==18) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","4.0/cosh(x/1.3)+3.0/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 0.3;
   }
   else if(modelnum==19) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = 0.3;
   }
   else if(modelnum==20) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.0/cosh(x/1.3)+2.0/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = 0.3;
   }
   else if(modelnum==21) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.0/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;   
      fClusNchild[0] = 2;
      fClusM[0] = 0.7; 
      fEtaBroaden = 1.1;
   }  
   else if(modelnum==22) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","1.4/cosh(x/1.3)+2.2/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;      
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 1.1;
   }  
   else if(modelnum==23) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.6/cosh(x/1.3)+1.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = 1.1;
   }
   else if(modelnum==24) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.4/cosh(x/1.3)+1.3/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = 1.1;
   }
   else if(modelnum==25) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.370)",0,10);
      fClusPzDist[0] = new TF1("pzdist","10.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
      fEtaBroaden = -0.4;
   }
   else if(modelnum==26) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.550)",0,10);
      fClusPzDist[0] = new TF1("pzdist","6.5/cosh(x/1.3)+3.5/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = -0.4;
   }
   else if(modelnum==27) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.750)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.5/cosh(x/1.3)+3.0/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
      fEtaBroaden = -0.4;
   }
   else if(modelnum==28) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.920)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.5/cosh(x/1.3)+2.7/cosh(x/5.5)+0.4/cosh(x/40.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      fEtaBroaden = -0.4;
   }
   else if(modelnum==33) { // pp 200GeV
      fConserveMom=1;

      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/1.80)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.8/cosh(x/15.0)",-200,200);
      fClusProb[0] = 0.04;
      fClusNchild[0] = 2;
      fClusM[0] = 0.5;

      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.600)",0,10);     
      fClusPzDist[1] = new TF1("pzdist","2/cosh(x/0.5)+1.8/cosh(x/2.5)"
                               "+0.2/cosh(x/12)+0.3/cosh(x/18)",-200,200);
      fClusProb[1] = 1.0;
      fClusNchild[1] = 3;
      fClusM[1] = 1.05;
      fEtaBroaden = 0.27;
   }
   else if(modelnum==34) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.450)",0,10); 
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.450)",0,10); 
      fClusPzDist[0] = new TF1("pzdist","3/cosh(x/0.5)+1.8/cosh(x/2.5)"
                               "+0.2/cosh(x/12)+0.3/cosh(x/18)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","3/cosh(x/0.5)+1.8/cosh(x/2.5)"
                               "+0.2/cosh(x/12)+0.3/cosh(x/18)",-200,200);
      fClusProb[0] = 1.4;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
      
      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 0.75;
   }
   else if(modelnum==35) { // CuCu 40%-50%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 2.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.32;
   }   
   else if(modelnum==36) { // CuCu 0%-10%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.520)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.0/cosh(x/1.3)+2.5/cosh(x/5.0)+0.8/cosh(x/20.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
      fEtaBroaden = 0.65;
   }  
   else if(modelnum==37) { // CuCu 40%-50% k+5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 6.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.32;
   }
   else if(modelnum==38) { // CuCu 40%-50% k-5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.1;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.32;
   }
   else if(modelnum==39) { // CuCu 40%-50% d+5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 2.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.44;
   }
   else if(modelnum==40) { // CuCu 40%-50% d-5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 2.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.20;
   }
   else if(modelnum==41) { // CuCu 40%-50% k+5%, d+5% 
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 6.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.44;
   }
   else if(modelnum==42) { // CuCu 40%-50% k-5%, d+5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.1;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.44;
   }
   else if(modelnum==43) { // CuCu 40%-50% k+5%, d-5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 6.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.0;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.20;
   }
   else if(modelnum==44) { // CuCu 40%-50% k-5%, d-5%
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","0.15/cosh(x/1.3)+2.4/cosh(x/4.8)+0.8/cosh(x/26.0)",-200,200);
      fClusPzDist[1] = new TF1("pzdist","0.35/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/24.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;

      fClusProb[1] = 1.1;
      fClusNchild[1] = 4;
      fClusM[1] = 1.4;
      fEtaBroaden = 1.20;
   }   
   else if(modelnum==45) { // K=1, momentum conserved
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.450)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.0/cosh(x/1.3)+2.5/cosh(x/5.0)+0.8/cosh(x/20.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 1;
      fClusM[0] = 0.35;
   }
   else if(modelnum==46) { //K=2, no momentum conserved
      fConserveMom=0;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.350)",0,10);
      fClusPzDist[0] = new TF1("pzdist","14.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.8/cosh(x/15.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 2;
      fClusM[0] = 0.7;
   }
   else if(modelnum==47) { //K=3, no momentum conserved
      fConserveMom=0;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.520)",0,10);
      fClusPzDist[0] = new TF1("pzdist","7.5/cosh(x/1.3)+2.5/cosh(x/5.0)+0.8/cosh(x/18.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 3;
      fClusM[0] = 1.05;
   }
   else if(modelnum==48) { //K=4, no momentum conserved
      fConserveMom=0;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.690)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.8/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/20.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
   }
   else if(modelnum==49) { //K=5, no momentum conserved
      fConserveMom=0;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.860)",0,10);
      fClusPzDist[0] = new TF1("pzdist","2.8/cosh(x/1.3)+2.4/cosh(x/4.7)+0.8/cosh(x/22.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 5;
      fClusM[0] = 1.75;
   }
   else if(modelnum==50) { // pp 200GeV, no momentum conserved
      fConserveMom=0;

      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/1.80)",0,10);
      fClusPzDist[0] = new TF1("pzdist","5.5/cosh(x/1.3)+2.5/cosh(x/5.5)+0.8/cosh(x/15.0)",-200,200);
      fClusProb[0] = 0.015;
      fClusNchild[0] = 2;
      fClusM[0] = 0.5;

      fClusPtDist[1] = new TF1("ptdist","x*exp(-x/0.430)",0,10);
      fClusPzDist[1] = new TF1("pzdist","2/cosh(x/0.5)+1.8/cosh(x/2.5)"
                               "+0.2/cosh(x/12)+0.3/cosh(x/18)",-200,200);
      fClusProb[1] = 1.0;
      fClusNchild[1] = 3;
      fClusM[1] = 1.05;
      fEtaBroaden = 0.27;
   }
   else if(modelnum==60) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.450)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.8/cosh(x/1.5)+2.5/cosh(x/4.3)+0.8/cosh(x/20.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
   }
   else if(modelnum==70) {
      fConserveMom=1;
      fClusPtDist[0] = new TF1("ptdist","x*exp(-x/0.450)",0,10);
      fClusPzDist[0] = new TF1("pzdist","3.8/cosh(x/(1.5+0.8-0.8))+2.5/cosh(x/(4.3+1.5-1.5))+0.8/cosh(x/20.0)",-200,200);
//      fClusPzDist[0] = new TF1("pzdist","3.8/cosh(x/[0])+2.5/cosh(x/[1])+0.8/cosh(x/20.0)",-200,200);
      fClusProb[0] = 1.0;
      fClusNchild[0] = 4;
      fClusM[0] = 1.4;
   }
   else{
      Info("FixModelParameters","undefined model num %d",modelnum);
      fErrorflag=1;
   }
}

Int_t TPhICM::GenClusType() { 
   Double_t random=gRandom->Rndm();
   for(Int_t i=0;i<6;i++) {
      if(random<=fClusProb[i]) return i;
   }
   return 0;
}

Int_t TPhICM::GenClusNChild(Int_t clustype) {
   return fClusNchild[clustype];
}

Double_t TPhICM::GenClusMass(Int_t clustype) {
//   if(fClusMf) return fClusMf[clustype]->GetRandom();
   return fClusM[clustype];
}

//Wei's New functions ------------------------------------------------
Double_t TPhICM::ProbPathLength(Double_t PathLength)
{
  return TMath::Exp(-PathLength/fAbsPar);
}

void TPhICM::GenInitEcc()
{
   fEccB=TMath::Sqrt(fNuclR*fNuclR-(fB/2.0)*(fB/2.0));
   fEccA=fNuclR-(fB/2.0);
}

void TPhICM::GenClusPos(Double_t& clusX, Double_t& clusY) {
  static TF1* func = new TF1("func","pol0(0)",-1,1);
  func->SetParameter(0,1);
  while(1)
  {
    clusX = func->GetRandom()*fEccA; 
    clusY = func->GetRandom()*fEccB;
    if(((clusX/fEccA)*(clusX/fEccA)+(clusY/fEccB)*(clusY/fEccB))<=1) break;
  }
}

//--------------------------------------------------------------------
void TPhICM::Run() {
   if(fErrorflag) return;

   for(Int_t ievent=1; ievent<=fNevents; ievent++) {
      if(!(ievent%100)) cout<<ievent<<" of "<<fNevents<<endl;
      RunEvent(ievent);
   }

      fOutTreef->cd("UEAnalysisRootple");
      fOutTree->Write();
      fOutTree->ls();
      fOutTreef->ls();
      delete fOutTreef;
}

void TPhICM::RunEvent(Int_t evtno) {
   if(fNuclR!=0 && fB<=2*fNuclR) GenInitEcc(); 

   vector<Cluster_t> clusterList;

   EventKind=1;
   Int_t nparticles=0;
   Int_t iMonteCarlo=0;
   fMultiplicity=fInitMultiplicity;

// only for dN/deta fluctuation study --------------------
//     fClusPzDist[0]->SetParameter(0,gRandom->Gaus(1.5,0.8));
//     fClusPzDist[0]->SetParameter(1,gRandom->Gaus(4.3,1.5));
// --------------------------------------------------------

   while(nparticles<fMultiplicity) {
      Int_t clustype=GenClusType();
      Int_t nchild=GenClusNChild(clustype);
      nparticles+=nchild;
      Double_t clusmass=GenClusMass(clustype);  
      Double_t cluspt=fClusPtDist[clustype]->GetRandom();
      Double_t cluspz=fClusPzDist[clustype]->GetRandom();
//Wei's
      Double_t clusX=0;
      Double_t clusY=0;
      Double_t clusZ=0;
      if(fNuclR!=0 && fB<=2*fNuclR) GenClusPos(clusX,clusY); 
      TVector3 dummyV(cluspt,0,cluspz);
      Double_t cluspx=cluspt*clusX/TMath::Sqrt(clusX*clusX+clusY*clusY);
      Double_t cluspy=cluspt*clusY/TMath::Sqrt(clusX*clusX+clusY*clusY);

      Cluster_t cluster={clusX,clusY,clusZ,cluspx,cluspy,cluspz,clusmass,nchild};
      clusterList.push_back(cluster);
   } 

   if(fConserveMom) {
      ConserveMomentumInEtaRange(clusterList,-8,8);
   }
   Int_t nclusters=clusterList.size();
   for(Int_t icluster=0;icluster<nclusters;icluster++) {
      Cluster_t cluster=clusterList[icluster];
      TLorentzVector *V=new TLorentzVector();
      V->SetXYZM(cluster.px,
         cluster.py,
         cluster.pz,
         cluster.m);

/*  
      if(fBetaRadial!=0) 
      {
         V->Boost(fBetaRadial*cluster.x/TMath::Sqrt(cluster.x*cluster.x+cluster.y*cluster.y),fBetaRadial*cluster.y/TMath::Sqrt(cluster.x*cluster.x+cluster.y*cluster.y),0);
      }
*/
      Int_t nchild=cluster.n;
      if(nchild==1) {
         new((*MonteCarlo)[iMonteCarlo]) TLorentzVector(*V);
         iMonteCarlo++;
     continue;
      }
      TGenPhaseSpace *clusterDecayer = new TGenPhaseSpace();
      
      Double_t masses[5]={0.14,0.14,0.14,0.14,0.14};
      if(nchild>5) {
     Info("TPhICM::ProduceCluster","nchild>5 not allowed");
     return;
      }
      clusterDecayer->SetDecay(*V,nchild,masses);
      Double_t wtmax=clusterDecayer->GetWtMax();
      while(1) {
     Double_t weight=clusterDecayer->Generate();
     Double_t wtlim=gRandom->Rndm()*wtmax;
     if(weight>wtlim) break;
      }
     
      Double_t etasum=0;
      for(Int_t ipart=0;ipart<nchild;ipart++) {
     TLorentzVector *particle = (TLorentzVector*)clusterDecayer->GetDecay(ipart)->Clone();
     if(fEtaBroaden) {
            Double_t etamod = (particle->Eta()-V->PseudoRapidity())*fEtaBroaden+particle->Eta(); 
        particle->SetPtEtaPhiE(particle->Pt(),etamod,
                   particle->Phi(),particle->E());
     }
         etasum = etasum + particle->Eta()/nchild;
//         htmp->Fill(particle->Eta()-V->PseudoRapidity());
         if(fAbsPar!=0)
         {
           Double_t k=TMath::Tan(particle->Phi());
           Double_t A = (1.0/fEccA)*(1.0/fEccA)+(k/fEccB)*(k/fEccB);
           Double_t B = (2*cluster.y*k-2*k*k*cluster.x)/(fEccB*fEccB);
           Double_t C = (cluster.y-k*cluster.x)*(cluster.y-k*cluster.x)/(fEccB*fEccB)-1;
           Double_t x1 = (-B+TMath::Sqrt(B*B-4*A*C))/(2*A);
           Double_t x2 = (-B-TMath::Sqrt(B*B-4*A*C))/(2*A);
           Double_t y1 = k*(x1-cluster.x)+cluster.y;
           Double_t y2 = k*(x2-cluster.x)+cluster.y;
  
           Double_t PathLength=0;
           Double_t xtrue=0;
           Double_t ytrue=0;
           if((y1-cluster.y)*particle->Py()>=0 && (x1-cluster.x)*particle->Px()>=0)
           {
              xtrue=x1;
              ytrue=y1;
           }
           else if((y2-cluster.y)*particle->Py()>=0 && (x2-cluster.x)*particle->Px()>=0)
           {
              xtrue=x2;
              ytrue=y2;
           }
           else cout<<"PathLength calculation error!"<<endl;
           PathLength=TMath::Sqrt((xtrue-cluster.x)*(xtrue-cluster.x)+(ytrue-cluster.y)*(ytrue-cluster.y));
           if(gRandom->Rndm()>ProbPathLength(PathLength)) {fMultiplicity--; continue;}   
         }
         new((*MonteCarlo)[iMonteCarlo]) TLorentzVector(*particle);
         iMonteCarlo++;
      }
      delete clusterDecayer;
      delete V;
   }

   WriteTreeEvent();
}

void TPhICM::ConserveMomentumInEtaRange(vector<Cluster_t> &clusterList,
                    Double_t etamin,Double_t etamax) {
   Int_t nclusters=clusterList.size();
   Double_t pxsum=0,pysum=0;
   Int_t n=0;
   for(Int_t icluster=0;icluster<nclusters;icluster++) {
      Cluster_t cluster=clusterList[icluster];
      TVector3 vector(cluster.px,cluster.py,cluster.pz);
      Double_t eta=vector.Eta();
      if(eta<etamin || eta>etamax) continue;
      pxsum+=cluster.px;
      pysum+=cluster.py;
      n++;
   }
   Double_t pxaver=pxsum/n;
   Double_t pyaver=pysum/n;
   
   for(Int_t icluster=0;icluster<nclusters;icluster++) {
      Cluster_t cluster=clusterList[icluster];
      TVector3 vect(cluster.px,cluster.py,cluster.pz);
      Double_t eta=vect.Eta();
      if(eta<etamin || eta>etamax) continue;
      vect.SetXYZ(cluster.px-pxaver,cluster.py-pyaver,cluster.pz);
      /*
      Double_t phi=vect.Phi();
      Double_t pt=vect.Pt();
      vect.SetPtEtaPhi(pt,phi,eta);
      */
      Cluster_t newcluster={cluster.x,
                            cluster.y,
                            cluster.z,
                            vect.x(),
                vect.y(),
                vect.z(),
                cluster.m,
                cluster.n};
      clusterList[icluster]=newcluster;
   }
}

void TPhICM::SetWriteTree(TString outfname) {

   if(!outfname.CompareTo("")) return;
   fOutTreef = new TFile(outfname.Data(),"recreate");
   if(!fOutTreef) {
      Info("SetWriteTree","could not create file");
      return;
   }

   fOutTreef->mkdir("UEAnalysisRootple");

   fOutTree=new TTree("AnalysisTree","MBUE Analysis Tree");
   fOutTree->Branch("EventKind",&EventKind,"EventKind/I");

   MonteCarlo = new TClonesArray("TLorentzVector", 10000);
   fOutTree->Branch("MonteCarlo", "TClonesArray", &MonteCarlo, 128000, 0);
}

void TPhICM::WriteTreeEvent() {
   fOutTree->Fill();
}
