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

/*
double getWeightfromCMSKIN(MissingET met,Lepton l1, Lepton l2, Jet j1, Jet j2)
{
  double quality=0.0;
  

  return quality;
} 
*/

int printJetConstituents(Jet* jet) {
  int genParticleIdx, trackIdx, towerIdx, candidateIdx, unknownIdx=0;
  std::cout<<"Jet nConstituents : "<< jet->Constituents.GetEntriesFast()<<std::endl;
  std::cout<<"track PID : ";
  for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
  {
    TObject* object = jet->Constituents.At(j);
    // Check if the constituent is accessible
    if(object == 0) continue;
    if ( object->IsA() == GenParticle::Class()) {
      GenParticle* particle = (GenParticle*) object;
      if ( particle->PT < 1 ) continue;
      D(std::cout<<"particle charge : "<<particle->Charge<<std::endl;)
      genParticleIdx++;
    }
    else if (object->IsA() == Track::Class())
    {
      Track* track = (Track*) object;
      if ( track->PT < 1 ) continue;
      trackIdx++;
      std::cout<< track->PID<<", ";
    }
    else if(object->IsA() == Tower::Class()) { 
      D(std::cout<<"Tower"<<std::endl;) 
      towerIdx++;
     }
    else if(object->IsA() == Candidate::Class()) { 
      D(std::cout<<"Candidate"<<std::endl; )
      candidateIdx++;
    }
    else {
      D(std::cout<<"Something diff."<<std::endl;)
      unknownIdx++;
    }
  }
  std::cout<<std::endl;
  return trackIdx;
}

void printJetParticles(Jet* jet) {
  std::cout<<"Jet nParticles : "<<jet->Particles.GetEntriesFast()<<std::endl;
  std::cout<<"Particle PID : "; 
  for(int i=0 ; i< jet->Particles.GetEntriesFast(); ++i) {
    GenParticle* cand = dynamic_cast<GenParticle*>(jet->Particles.At(i));
    if ( cand ==nullptr)  continue;
    if ( cand->PT < 1 ) continue;
    std::cout<< cand->PID<<"("<<cand->Status<<") , ";
  }
  std::cout<<std::endl;


  return ;
}
/*
double getJetCharge(Jet* jet,bool weight) {
  double charge =0.;
  double weightSum =0.0;
  int trackIdx=0;
  for(int j = 0; j < jet->Constituents.GetEntriesFast(); ++j)
  {
    TObject* object = jet->Constituents.At(j);
    // Check if the constituent is accessible
    if(object == 0) continue;


    if ( object->IsA() == GenParticle::Class()) {
      GenParticle* particle = (GenParticle*) object;
      if ( particle->PT < 1 ) continue;
      D(std::cout<<"particle charge : "<<particle->Charge<<std::endl;)
      if(weight) {
        weightSum += particle->PT;
        charge += particle->PT*particle->Charge;
      }
      else {
        weightSum +=1.0;
        charge += particle->Charge;
      }
    }
    else if (object->IsA() == Track::Class())
    {
      Track* track = (Track*) object;
      if ( track->PT < 1 ) continue;
      trackIdx++;
      if ( trackIdx==1) {
        std::cout<<std::endl;
        std::cout<<"Jet PT : "<<jet->PT<<" Jet eta: "<<jet->Eta<<" Jet phi : "<<jet->Phi<<std::endl;
      }
  
      std::cout<<trackIdx<<"th's track charge : "<<track->Charge<<std::endl;
      std::cout<<"track PT : "<<track->PT<<" track Eta : "<<track->Eta<<"  track phi : "<<jet->Phi<<std::endl;
      if( weight) {
        weightSum += track->PT;
        charge += track->PT*track->Charge;
      }
      else {
        weightSum += 1.0;
        charge += 1.0*track->Charge;
      }
    }
  }
  return charge/weightSum;
}
*/
struct TestPlots
{
 
  /* 
  TH1 *fgenJpsiPT;
  TH1 *fgenJpsiEta;
  TH1 *fgenJpsiPhi;
  TH1 *fgenJpsiM;

  TH1 *fJpsiPT;
  TH1 *fJpsiEta;
  TH1 *fJpsiPhi;
  */
  TH1 *fJetCharge;
  TH1 *fnJetTracks;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  /*
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
  */
  plots->fJetCharge = result->AddHist1D(
    "jetCharge", "charge",
    "Jet Charge", "Jet Charge",
    10, -5.0, 5.0 );
  plots->fnJetTracks = result->AddHist1D(
    "nJetTracks", "Track distribution",
    "Number of Track", "Entries",
    10, 0.0, 10.0 );
  /*
  plots->fJetChargeSub = result->AddHist1D(
    "jetChargeSub", "Jet Charge using subTrack",
    "Jet Charge", "Entries",
    10, -5.0, 5.0 );
  plots->fJetChargeWSub = result->AddHist1D(
    "jetChargeWeightedSub", "Jet weighted Charge using subtrack",
    "Jet Charge", "Entries",
    10, -5.0, 5.0 );
  plots->fnJetTracksSub = result->AddHist1D(
    "nJetTracksSub", "Track distribution using subTrack",
    "Number of Track", "Entries",
    10, 0.0, 10.0 );
  */
}

//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  TClonesArray *branchJetSubTrack = treeReader->UseBranch("JetSubTracks");
  TClonesArray *branchJetSubPhoton = treeReader->UseBranch("JetSubPhotons");
  TClonesArray *branchJetSubNeutralHadron = treeReader->UseBranch("JetSubNeutralHadrons");

  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  D(cout << "** Chain contains " << allEntries << " events" << endl;)

  GenParticle *particle;


  //Int_t i, j, entry;
  int percent=0;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if (entry % (allEntries/100)==0 ) std::cout<<"Event : "<<entry<<"("<<percent++<<"%)"<<std::endl;
    DecayChannel dc(branchParticle);
    if ( dc.channel() !=0 ) continue;
    D( std::cout<<"This is dilepton events!"<<std::endl; )

    for(int i = 0; i < branchJet->GetEntriesFast(); ++i)
    {
      Jet* jet = (Jet*) branchJet->At(i);
      if ( jet->PT < 30 || abs(jet->Eta)>2.4 ) continue;
      double jetCharge = jet->Charge;
      int nJetTracks = jet->NCharged;
      D(std::cout<<"Jet Charge : "<<jetCharge<<std::endl;)
      plots->fJetCharge->Fill( jetCharge ) ;
      D(std::cout<<"njetTrack : "<<nJetTracks<<std::endl; )
      plots->fnJetTracks->Fill( nJetTracks ) ;
      //if ( jet->Flavor !=5 ) continue;
      //std::cout<<"This is bjet!"<<std::end;
      printJetConstituents(jet); 
      printJetParticles(jet);    
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------


int main(int argc, char* argv[])
{
  if ( argc !=3 ) {
    std::cout<<"Argument is wrong. Please, run \"JetCharge Delphes_input.root result_output.root\""<<std::endl;
    return -1;
  }
  const char* inputFile = argv[1];
  const char* outputFile = argv[2];

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write(outputFile);

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

}
