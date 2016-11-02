#ifndef NDEBUG
#define D(x) x 
#else
#define D(x)  
#endif

#include<iostream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootTreeWriter.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <TClonesArray.h>
#include <TH1.h>
#include <TH2.h>
#include <DecayChannel.h>
#include <KinematicSolvers.h>
#include <analysisUtils.h>
#include <TRandom3.h>
#include <TCanvas.h>
using namespace std;

class ExRootTreeReader;


//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");

  Long64_t allEntries = treeReader->GetEntries();

  D(cout << "** Chain contains " << allEntries << " events" << endl;)


  TH1F* h1 = new TH1F("top_mass","top_mass",200,160,180);
  //Int_t i, j, entry;
  int percent=0;
  for(int entry = 0; entry < allEntries; ++entry) {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if (entry % (allEntries/100)==0 ) std::cout<<"Event : "<<entry<<"("<<percent++<<"%)"<<std::endl;

    // Calculate Weight!
    for( unsigned int i = 0 ; i < branchParticle->GetEntries() ; ++i) {
        GenParticle* particle = (GenParticle*)branchParticle->At(i);
        if ( abs(particle->PID) ==6) std::cout<<particle->Mass<<std::endl;
        h1->Fill(particle->Mass);

    }
  }
  TCanvas* c1 = new TCanvas("c1","c1",600,600);
  h1->Draw();
  c1->SaveAs("partonTop_mass.png");
}




  

int main(int argc, char* argv[])
{
  if ( argc !=2 ) {
    std::cout<<"Argument is wrong. Please, run \"JetCharge Delphes_input.root\""<<std::endl;
    return -1;
  }
  const char* inputFile = argv[1];

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  AnalyseEvents(treeReader);
  cout << "** Exiting..." << endl;

  //delete result;
  delete treeReader;
  delete chain;

}
