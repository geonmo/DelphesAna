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
  TH1 *fgenJpsiPT;
  TH1 *fgenJpsiEta;
  TH1 *fgenJpsiPhi;
  TH1 *fgenJpsiM;

  TH1 *fJpsiPT;
  TH1 *fJpsiETA;
  TH1 *fJpsiPHI;
  TH1 *fJpsiM;

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  //TLegend *legend;
  //TPaveText *comment;

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


void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;

  Jet *jet;
  TObject *object;

  //TLorentzVector momentum;

  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    bool isJpsiExisted = false;
    for ( i =0 ; i < branchParticle->GetEntriesFast() ; ++i) {
      particle = (GenParticle*) branchParticle->At(i);
      if ( abs( particle->PID ) == 443 ) {
        GenParticle* base  = particle; 
        while(base->M1 != -1 ) {
          GenParticle* mother = (GenParticle*) branchParticle->At(base->M1);
          if ( abs(base->PID) == 6 ) { 
            //std::cout<<"Find Top from J/psi's mother"<<stD::endl; 
            isJpsiExisted = true;
            std::cout<<particle->PT<<" "<<particle->Eta<<" "<<particle->Phi<<" "<<particle->Mass<<std::endl;
            plots->fgenJpsiPT->Fill(particle->PT);
            plots->fgenJpsiEta->Fill(particle->Eta);
            plots->fgenJpsiPhi->Fill(particle->Phi);
            plots->fgenJpsiM->Fill(particle->Mass);
            break;   
          }
          base = mother ; 
        }
      }
    }
    if ( !isJpsiExisted )  continue; 
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void AnaJpsi(const char *inputFile="Delphes.root")
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
