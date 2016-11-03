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


  //Int_t i, j, entry;
  int percent=0;
  for(int entry = 0; entry < allEntries; ++entry) {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if (entry % (allEntries/100)==0 ) std::cout<<"Event : "<<entry<<"("<<percent++<<"%)"<<std::endl;

    // Calculate Weight!
    for( unsigned int i = 0 ; i < branchParticle->GetEntries() ; ++i) {
        GenParticle* particle = (GenParticle*)branchParticle->At(i);
        if ( abs(particle->PID) ==421) {
          GenParticle* mother = (GenParticle*)branchParticle->At(particle->M1); 
          if ( abs(mother->PID) ==413 ) {
            int nDau = particle->D2-particle->D1+1;
            std::cout<<"nDau : "<<nDau<<std::endl;
            std::cout<<mother->PID<<" to ("<<particle->PID<<") to ";
            for(int i=0 ; i< nDau ; i++) {
              if ( particle->D1+i != -1) {
                GenParticle* dau = (GenParticle*)branchParticle->At(particle->D1+i); 
                std::cout<<dau->PID<<" + ";
              }
            }
            std::cout<<std::endl;
          }
          /*
          if ( abs(mother->PID)>500 && abs(mother->PID)<600 ) {
            if ( mother->D2 != -1 && mother->D1 != -1) {
              std::cout<<"nDau : "<<abs(mother->D2-mother->D1)+1<<"  from "<<mother->PID<<std::endl;
              std::cout<<mother->PID<<" to ("<<particle->PID<<")";
              for(int i=0 ; i< mother->D2-mother->D1+1 ; i++) {
                GenParticle* dau = (GenParticle*)branchParticle->At(mother->D1+i); 
                std::cout<<dau->PID<<" + ";
              }
              std::cout<<std::endl;
            }
          }
          */
        }

    }
  }
}




  

int main(int argc, char* argv[])
{
  if ( argc !=2 ) {
    std::cout<<"Argument is wrong. Please, run \"checkD0 Delphes_input.root\""<<std::endl;
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
