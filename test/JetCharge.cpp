#ifndef NDEBUG
#define D(x) x 
#else
#define D(x)  
#endif

#include<iostream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <TClonesArray.h>
#include <TH1.h>
#include <DecayChannel.h>
using namespace std;

class ExRootTreeReader;
class ExRootResult;

struct TestPlots
{
  TH1 *fgenJpsiPT;
  TH1 *fgenJpsiEta;
  TH1 *fgenJpsiPhi;
  TH1 *fgenJpsiM;

  TH1 *fJpsiPT;
  TH1 *fJpsiEta;
  TH1 *fJpsiPhi;

};

//------------------------------------------------------------------------------

int getJetCharge(Jet* jet) {
  double charge =0.;
  for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
  {
    TObject* object = jet->Constituents.At(j);
    // Check if the constituent is accessible
    if(object == 0) continue;

    if (object->IsA() == Track::Class())
    {
      Track* track = (Track*) object;
      charge += track->Charge;
      if ( track->PT < 1 ) continue;
    }
  }
  return charge;
}
void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  plots->fgenJpsiPT = result->AddHist1D(
    "genjpsi_pt", "p_{T}",
    "p_{T}", "Entries",
    100, 0., 200. );
  plots->fgenJpsiEta = result->AddHist1D(
    "genjpsi_eta", "eta",
    "eta", "number of jpsis",
    100, -5., 5. );
  plots->fgenJpsiPhi = result->AddHist1D(
    "genjpsi_phi", "phi",
    "phi", "number of jpsis",
    100, -3.14, +3.14 );
  plots->fgenJpsiM = result->AddHist1D(
    "genjpsi_mass", "mass",
    "mass", "number of jpsis",
    100, 0., 5. );

  plots->fJpsiPT = result->AddHist1D(
    "jpsi_pt", "p_{T}",
    "p_{T}", "number of Jpsi",
    100, 0, 200 );
  plots->fJpsiEta = result->AddHist1D(
    "jpsi_eta", "eta",
    "eta", "number of jpsis",
    100, -5, 5 );
  plots->fJpsiPhi = result->AddHist1D(
    "jpsi_phi", "phi",
    "phi", "number of jpsis",
    100, -TMath::Pi(), TMath::Pi() );

}

//------------------------------------------------------------------------------
int FindMother(TClonesArray* genParticles, GenParticle* base, int motherPID )
{
  int value = 0;
  while(base->M1 != -1 ) {
    GenParticle* mother = (GenParticle*) genParticles->At(base->M1);
    if ( abs(base->PID) == motherPID ) value = base->PID;
    base = mother;
  }
  return value;
}
void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  //TClonesArray *branchTrack = treeReader->UseBranch("JetSubTracks");

  TClonesArray *branchJetSubTrack = treeReader->UseBranch("JetSubTracks");
  TClonesArray *branchJetSubPhoton = treeReader->UseBranch("JetSubPhotons");
  TClonesArray *branchJetSubNeutralHadron = treeReader->UseBranch("JetSubNeutralHadrons");
 
  /* 
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  */
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;


  //Int_t i, j, entry;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if ( allEntries < 1000 || entry%1000) std::cout<<"Event : "<<entry<<std::endl;
    D(cout<<"Process for GenParticle"<<endl;)
    DecayChannel dc(branchParticle);
    if ( dc.channel() !=0 ) continue;
    D( std::cout<<"This is dilepton events!"<<std::endl; )

    for(int i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      Jet* jet = (Jet*) branchJet->At(i);
      double jetCharge = getJetCharge( jet) ;
      std::cout<<"Jet Charge : "<<jetCharge<<std::endl;
      
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------


int main()
{
  const char* inputFile = "Delphes.root";

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
