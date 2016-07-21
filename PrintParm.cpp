#include<iostream>
#include<fstream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <TClonesArray.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
using namespace std;

class ExRootTreeReader;
class ExRootResult;


struct Data {
  unsigned int nLep;
  float lep1_pt, lep1_iso;
  float lep2_pt, lep2_iso;
  float lep3_pt, lep3_iso;
  float lep4_pt, lep4_iso;

  float dilep1_pt, dilep1_mass;
  float dilep2_pt, dilep2_mass;
  float dilep3_pt, dilep3_mass;
  float dilep4_pt, dilep4_mass;

  float met, met_phi;

  unsigned int nJet;
  float jet1_pt, jet1_eta, jet1_phi;
  float jet2_pt, jet2_eta, jet2_phi;
  float jet3_pt, jet3_eta, jet3_phi;
  float jet4_pt, jet4_eta, jet4_phi;
  unsigned int nbJet;
  float bjet1_pt, bjet1_eta, bjet1_phi, bjet1_mva;
  float bjet2_pt, bjet2_eta, bjet2_phi, bjet2_mva;

  int type;

  Data(){
    reset();
  }
  void reset() {
    nLep =0;
    lep1_pt=0.0; lep1_iso=0.0; 
    lep2_pt=0.0; lep2_iso=0.0;
    lep3_pt=0.0; lep3_iso=0.0;
    lep4_pt=0.0; lep4_iso=0.0;

    dilep1_pt =0.0; dilep1_mass =0.0;
    dilep2_pt =0.0; dilep2_mass =0.0;
    dilep3_pt =0.0; dilep3_mass =0.0;
    dilep4_pt =0.0; dilep4_mass =0.0;

    met=0.0; met_phi=0.0;

    nJet=0;
    jet1_pt=0.0, jet1_eta=-999., jet1_phi=0.0;
    jet2_pt=0.0, jet2_eta=-999., jet2_phi=0.0;
    jet3_pt=0.0, jet3_eta=-999., jet3_phi=0.0;
    jet4_pt=0.0, jet4_eta=-999., jet4_phi=0.0;

    nbJet=0;
    bjet1_pt=0.0, bjet1_eta=-999, bjet1_phi=0.0, bjet1_mva=0.0;
    bjet2_pt=0.0, bjet2_eta=-999, bjet2_phi=0.0, bjet2_mva=0.0;
    type = 0;
  }
}; 
  


void BookingTree(TTree* tree, Data& data) {
  tree->Branch("data",&data,"nLep/i:lep1_pt/F:lep1_iso:lep2_pt:lep2_iso:lep3_pt:lep3_iso:lep4_pt:lep4_iso:dilep1_pt:dilep1_mass:dilep2_pt:dilep2_mass:dilep3_pt:dilep3_mass:dilep4_pt:dilep4_mass");
}


void testTree(TTree* tree){
  Data data;
  BookingTree(tree, data);
  for(int i=0; i<10 ; i++) {
    data.nLep=i;
    data.lep1_pt= (float)i;
    data.lep1_iso= (float)i;
    data.lep2_pt= (float)i;
    data.lep2_iso= (float)i;
    data.lep3_pt= (float)i;
    data.lep3_iso= (float)i;
    data.lep4_pt= (float)i;
    data.lep4_iso= (float)i;
    tree->Fill();
  }
}


void PrintParameters(ofstream& outfile, TTree* tree)
{



}

int FindWboson(TClonesArray* genParticles, int baseIdx )
{
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 24 ) return baseIdx;

  // First, D1
  if ( base->D1 !=-1 ) {
    int nextResult = FindWboson( genParticles, base->D1);
    if ( nextResult !=-1 ) return nextResult;
  }
  // Next, D2
  if ( base->D2 !=-1 ) {
    int nextResult = FindWboson( genParticles, base->D2);
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}
int FindLepton(TClonesArray* genParticles, int baseIdx )
{
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 11 || absPID == 13 || absPID == 15 ) return baseIdx;

  // First, D1
  if ( base->D1 !=-1 ) {
    int nextResult = FindLepton( genParticles, base->D1);
    if ( nextResult !=-1 ) return nextResult;
  }
  // Next, D2
  if ( base->D2 !=-1 ) {
    int nextResult = FindLepton( genParticles, base->D2);
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}

int channelSelection( TClonesArray* genParticles) {
  // searcing Top or anti top
  int nGenParticle = genParticles->GetEntriesFast();

  int top_idx =0, antitop_idx=0;
  
  for (int  i =0 ; i < genParticles->GetEntriesFast() ; ++i) {
    GenParticle* genParticle = (GenParticle*) genParticles->At(i);
    if ( top_idx==0 && genParticle->PID == 6) top_idx= i;
    if ( antitop_idx==0 && genParticle->PID == -6 ) antitop_idx= i;
    if ( top_idx !=0 && antitop_idx != 0 ) break;
  }

  int WbosonIdx1 = FindWboson(genParticles,     top_idx);
  int WbosonIdx2 = FindWboson(genParticles, antitop_idx);

  bool     topToLepton = false;
  bool antitopToLepton = false;

  int lep1Idx = FindLepton( genParticles, WbosonIdx1);
  int lep2Idx = FindLepton( genParticles, WbosonIdx2);
  int lep1=-1, lep2=-1;
  if ( lep1Idx != -1) lep1 = ((GenParticle*) genParticles->At( lep1Idx) )->PID;
  if ( lep2Idx != -1) lep2 = ((GenParticle*) genParticles->At( lep2Idx) )->PID;

  std::cout<<std::endl;
  int trackingIdx;
  if ( lep1Idx != -1 ) {
    trackingIdx = lep1Idx;
    int upper_count = 0;
    while( 1 ) {
      GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
      std::cout<< base->PID;
      if ( abs(base->PID) == 24 ) break;
      std::cout<<">>";
      trackingIdx = base->M1;
      upper_count++;
    }
    if ( upper_count >1 ) { std::cout<<"oops"<<std::endl; lep1 = -1; }
  } 
  std::cout<<std::endl;
  if ( lep2Idx != -1 ) {
    trackingIdx = lep2Idx;
    int upper_count = 0;
    while( 1 ) {
      GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
      std::cout<< base->PID;
      if ( abs(base->PID) == 24 ) break;
      std::cout<<">>";
      trackingIdx = base->M1;
      upper_count++;
    }
    if ( upper_count >1 ) { std::cout<<"oops"<<std::endl; lep2 = -1; }
  } 
  std::cout<<std::endl;

  int mulValue = lep1*lep2;

  if ( abs(mulValue) > 100 ) { 
    if ( abs(mulValue) %15 ==0 ) {
      std::cout<<"Dilepton tau"<<std::endl;
      return 4;
    }
    std::cout<<"Dilepton"<<std::endl;
    std::cout<<"lep1 :"<<lep1<<"  lep2 : "<<lep2<<std::endl;
    return 1;
  }
  else if ( abs(mulValue) > 10 ) {
    if ( abs(mulValue) %15 ==0 ) {
      std::cout<<"Semilepton tau"<<std::endl;
      return 5;
    }
    std::cout<<"Semilepton"<<std::endl;
    return 2;
  }
  else {
    std::cout<<"Hardronic"<<std::endl;
    return 3; 
  }

}


void AnalyseEvents(ExRootTreeReader *treeReader, bool isTT = false )
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;


  int channel[6]={0,0,0,0,0,0};
  //allEntries = 1;
  //Int_t i, j, entry;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    channel[0]++;
    channel[ channelSelection( branchParticle )]++;
  }
  for( int i= 0 ; i < 6 ; i++) {
    std::cout<<"Channel : "<<i<<" "<<channel[i]<<"\t"<< (float)channel[i] / (float)channel[0]*100<<"%"<<std::endl;
  }
}


int main(int argc, char* argv[])
{
  bool isTT = false;

  if ( argc != 3 && argc !=4  ) {
    std::cout<<"Wrong argument!"<<std::endl;
    exit(-1);
  }
  
  std::cout<<argv[0]<<std::endl;
  if ( argc ==4 ) std::cout<<argv[3]<<std::endl;
  std::string inputFile(argv[1]);
  std::string outFile(argv[2]);
  
  if ( argc == 4 && std::atoi(argv[3]) == 1 ) isTT = true;

  std::cout<<"Input : "<<inputFile << "\tOut : "<<outFile<<std::endl;
  
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile.c_str());

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TFile* file = new TFile(outFile.c_str(), "RECREATE");
  TTree* tree = new TTree("delphes","delphes");
  AnalyseEvents(treeReader, isTT );
  testTree(tree);

  file->Write();
  file->Close();
  

  cout << "** Exiting..." << endl;

}

