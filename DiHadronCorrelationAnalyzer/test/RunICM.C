#include "../interface/TPhICM.h"
#include "TString.h"

void RunICM(int nevents, int seed, int multiplicity, int modnum){
TPhICM a(nevents,seed,multiplicity,modnum);
TString outputname = Form("/net/hisrv0001/home/davidlw/scratch/DiHadronCorrelations/icm/mod%d/icm_run%d.root",modnum,seed);
a.SetWriteTree(outputname);
a.Run();
}
