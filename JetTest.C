/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example3.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

//------------------------------------------------------------------------------

struct TestPlots
{

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");

  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;


  //Int_t i, j, entry;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    std::cout<<"Event : "<<entry<<std::endl;
    for(int i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      Jet* jet = (Jet*) branchJet->At(i);
      bool notice = false;
      int genPID =0;
      for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
      {
        TObject* object = jet->Constituents.At(j);
        // Check if the constituent is accessible
        if(object == 0) continue;

        if(object->IsA() == GenParticle::Class())
        {
          particle = (GenParticle*) object;
          if( abs(particle->PID) == 11 || abs(particle->PID==13) ) { 
            if ( !notice ) { std::cout<<"Jet Index : "<<i<<" Jet pT:"<<jet->PT<<"  Jet Eta "<<jet->Eta<<"  Jet Phi "<<jet->Phi<<std::endl; notice = true; }
            cout << "    GenPart pt: " << particle->PT << ", eta: " << particle->Eta << ", phi: " << particle->Phi << ", PID:"<<particle->PID<<endl;
          }
        }
        else if(object->IsA() == Track::Class())
        {
          Track* track = (Track*) object;
          GenParticle* trackGen = (GenParticle*)track->Particle.GetObject();
          genPID = trackGen->PID;
          if ( abs(genPID) != 11 && abs( genPID) != 13) continue;
          if ( track->PT < 1 ) continue;
          if ( !notice ) { std::cout<<"Jet Index : "<<i<<" Jet pT:"<<jet->PT<<"  Jet Eta "<<jet->Eta<<"  Jet Phi "<<jet->Phi<<std::endl; notice = true; }
          cout << "    Track pt: " << track->PT << ", eta: " << track->Eta << ", phi: " << track->Phi << " , PID: "<<track->PID<<" ,  GenPID: "<<genPID<<endl;
        }
        else if(object->IsA() == Tower::Class())
        {
          Tower* tower = (Tower*) object;
          //cout << "    Tower pt: " << tower->ET << ", eta: " << tower->Eta << ", phi: " << tower->Phi << endl;
          for( int k = 0 ; k< tower->Particles.GetEntriesFast() ; ++k) {
            GenParticle* towerGen = (GenParticle*)tower->Particles->At(k);
            genPID = towerGen->PID;
            if ( abs(genPID) != 11 && abs( genPID) != 13) continue;
            //if ( !notice ) { std::cout<<"Jet Index : "<<i<<" Jet pT:"<<jet->PT<<"  Jet Eta "<<jet->Eta<<"  Jet Phi "<<jet->Phi<<std::endl; notice = true; }
            //cout << "       Tower GENPart : " << towerGen->PT << ", eta: " << towerGen->Eta << ", phi: " << towerGen->Phi << ", GenPID: "<<towerGen->PID<<endl;
          }
        }
        else cout<<"Other"<<std::endl;
      }
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void JetTest(const char *inputFile="Delphes.root")
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
