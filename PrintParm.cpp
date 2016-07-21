#ifndef NDEBUG
#define D(x) x 
#else
#define D(x) 
#endif

#include<iostream>
#include<fstream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <TClonesArray.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include "DecayChannel.h"
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
  tree->Branch("data",&data,"nLep/i:lep1_pt/F:lep1_iso:lep2_pt:lep2_iso:lep3_pt:lep3_iso:lep4_pt:lep4_iso:dilep1_pt:dilep1_mass:dilep2_pt:dilep2_mass:dilep3_pt:dilep3_mass:dilep4_pt:dilep4_mass:met:met_phi:nJet/i:jet1_pt/F:jet1_eta:jet1_phi:jet2_pt:jet2_eta:jet2_phi:jet3_pt:jet3_eta:jet3_phi:jet4_pt:jet4_eta:jet4_phi:nbJet/i:bjet1_pt/F:bjet1_eta:bjet1_phi:bjet1_mva:bjet2_pt/F:bjet2_eta:bjet2_phi:bjet2_mva:type/i");
}


void putData(Data& data, int i){
  data.nLep=i;
  data.lep1_pt= (float)i;
  data.lep1_iso= (float)i;
  data.lep2_pt= (float)i;
  data.lep2_iso= (float)i;
  data.lep3_pt= (float)i;
  data.lep3_iso= (float)i;
  data.lep4_pt= (float)i;
  data.lep4_iso= (float)i;
  data.bjet2_eta = (float)i;
  data.bjet1_mva = (float)i;
}


void PrintParameters(ofstream& outfile, TTree* tree)
{



}

void AnalyseEvents(ExRootTreeReader *treeReader, TTree* tree, int dataType )
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Data data;
  BookingTree(tree, data);

  int channel[6]={0,0,0,0,0,0};
  //allEntries = 1;
  //Int_t i, j, entry;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    data.reset();
    treeReader->ReadEntry(entry);
    channel[0]++;
    if ( dataType < 4 ) { 
      DecayChannel dc(branchParticle);
      data.type = dc.channel();
      channel[ dc.channel()]++;
    }

    putData(data, entry);
    tree->Fill();
  }

  if ( dataType < 4 ) {
    for( int i= 0 ; i < 6 ; i++) {
      std::cout<<"Channel : "<<i<<" "<<channel[i]<<"\t"<< (float)channel[i] / (float)channel[0]*100<<"%"<<std::endl;
    }
  }
}


int main(int argc, char* argv[])
{

  int dataType=0;
  if ( argc != 3 && argc !=4  ) {
    std::cerr<<"Wrong argument!"<<std::endl;
    exit(-1);
  }

  D(
  std::cout<<argv[0]<<std::endl;
  if ( argc ==4 ) std::cout<<argv[3]<<std::endl;
  )
  std::string inputFile(argv[1]);
  std::string outFile(argv[2]);


  if ( argc == 4 ) dataType = std::atoi(argv[3]) ;

  D(std::cout<<"Input : "<<inputFile << "\tOut : "<<outFile<<std::endl;)

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile.c_str());

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TFile* file = new TFile(outFile.c_str(), "RECREATE");
  TTree* tree = new TTree("delphes","delphes");
  AnalyseEvents(treeReader, tree, dataType );

  file->Write();
  file->Close();


  cout << "** Exiting..." << endl;

}

