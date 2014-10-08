#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include <TEnv.h>
#include <TChain.h>
#include <TList.h>
#include <TString.h>
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>
#endif


#ifndef _RUNANA_
#define _RUNANA_
vector<string> GetDataVector(TList *files, Int_t ffrom=0, Int_t fto=-1);
TChain* GetDataChain(TList *files, Int_t ffrom=0,
	      Int_t fto=-1, const Char_t *treetype="AnT");
TList* GetFileList(const Char_t* fname);
TList* GetFileListFromFile(const Char_t* filename);
TList* GetFileListFromDir(const Char_t* dirname);
#endif

TChain* GetDataChain(TList *files, Int_t ffrom, 
	      Int_t fto, const Char_t *treetype)
{
   TChain *dataChain=new TChain(treetype); 
   if (dataChain==0) {
      cerr << "Could not create dataChain." << endl;
      return 0;
   }

   if (files==0) {
      cerr << "File list is zero." << endl;
      return 0;
   }

   Int_t filefrom=ffrom;
   Int_t fileto=files->GetEntries();
   if ((fto>0) && (fto<fileto)) fileto=fto;
   // parse file list
   for(Int_t i=filefrom;i<fileto;i++) {
      if(!(files->At(i))){
	 cerr << "Problem at object number: " << i << endl;
	 break;
      }
      TObjString* objstr = (TObjString*) files->At(i); 
      if (objstr!=0) {
	 TString filename = objstr->GetString();
         cout<<filename<<endl;
	 dataChain->Add(filename);
	 if (gDebug) cout << "Added: " << filename.Data() << endl;
      }
   }
   return dataChain;   
}

vector<string> GetDataVector(TList *files, Int_t ffrom, Int_t fto)
{
   vector<string> dataVector;

   Int_t filefrom=ffrom;
   Int_t fileto=files->GetEntries();
   if ((fto>0) && (fto<fileto)) fileto=fto;
   // parse file list
   for(Int_t i=filefrom;i<fileto;i++) {
      if(!(files->At(i))){
         cerr << "Problem at object number: " << i << endl;
         break;
      }
      TObjString* objstr = (TObjString*) files->At(i);
      if (objstr!=0) {
         TString filename = objstr->GetString();
         cout<<filename<<endl;
         string a(filename.Data());
         dataVector.push_back(a);
         if (gDebug) cout << "Added: " << filename.Data() << endl;
      }
   }
   return dataVector;
}

TList* GetFileList(const Char_t* fname)
{
   TString filesetname(fname);
   if (filesetname.Contains(":")){
      TObjArray* obj=filesetname.Tokenize(":");
      TString a=((TObjString*)obj->At(0))->GetString();
      TString b=((TObjString*)obj->At(1))->GetString();
      delete obj;
      if(a.CompareTo("file")==0)  
	 return GetFileListFromFile(b.Data());
      else if(a.CompareTo("dir")==0)  
	 return GetFileListFromDir(b.Data());
      else return 0; 
   } else {
      return 0;
   }
}

TList* GetFileListFromFile(const Char_t* filename)
{
   TString name = gSystem->ExpandPathName(filename);
   ifstream *in=new ifstream(name);

   if(!in) {
      cerr<< "File " << name.Data() << " could not be openned!" << endl;
      return 0;
   }

   if (!in->is_open()){
      cerr<< "File " << name.Data() << " could not be openned!" << endl;
      delete in;
      return 0;
   }

   TList *list = new TList;
   Char_t buf[4096];
   *in>>buf;
   while(!in->eof()){ 
      TString ifn = buf; 
      if (ifn.Contains(".root")) 
      {
	 TObjString *str= new TObjString(ifn);
	 list->Add(str);
	 if(gDebug) cout << "Added file " << buf << endl;
      }  
      *in>>buf;
   }

   delete in;
   return list;
}

TList* GetFileListFromDir(const Char_t* dirname)
{
   TString indir = gSystem->ExpandPathName(dirname);
   void* dir = gSystem->OpenDirectory(indir);
   if (dir==0) return 0;

   TList *list = new TList;
   Char_t* fn;
   while ( (fn=(Char_t*)gSystem->GetDirEntry(dir)) ) {
      TString ifn = indir; 
      ifn += fn;
      if (ifn.Contains(".root")) 
      {
	 TObjString *str= new TObjString(ifn);
	 list->Add(str);
      }  
   }
   return list;
}
