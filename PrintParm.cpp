#ifndef NDEBUG
#define D(x) x 
#else
#define D(x)  
#endif

#include<iostream>
#include<fstream>
#include<sstream>
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
  int lep_q[4];
  float lep_pt[4];
  float lep_iso[4];

  float dilep_pt[4], dilep_mass[4];

  float met, met_eta, met_phi;

  unsigned int nJet;
  float jet_pt[4], jet_eta[4], jet_phi[4];
  int jet_btag[4];
  unsigned int nbJet;

  int type;

  Data(){
    reset();
  }
  void reset() {
    nLep =0;
    for(int i=0 ; i< 4 ; ++i) {
      lep_q[i]=0;
      lep_pt[i]=0.0f;
      lep_iso[i]=0.0f;
  
      dilep_pt[i] = 0.0f;
      dilep_mass[i] = 0.0f;

      jet_pt[i] = 0.0f;
      jet_eta[i] = 0.0f;
      jet_phi[i] = 0.0f;
      jet_btag[i] = 0;
    }

    met=0.0f; met_eta=0.0f; met_phi=0.0f;

    nJet=0;

    nbJet=0;
    type = 0;
  }

  std::string print() {
    ostringstream ost;
    ost<<nLep<<",";
    for ( int i=0 ; i< 4 ; ++i) {
      ost<<lep_pt[i]<<","<<lep_iso[i]<<","<<lep_q[i]<<",";
    }
    for ( int i=0 ; i< 4 ; ++i) {
      ost<<dilep_pt[i]<<","<<dilep_mass[i]<<",";
    }

    ost<<met<<","<<met_eta<<","<<met_phi<<",";

    ost<<nJet<<",";
    for ( int i=0 ; i< 4 ; ++i) {
      ost<<jet_pt[i]<<","<<jet_eta[i]<<","<<jet_phi[i]<<","<<jet_btag[i]<<",";
    }

    ost<<nbJet<<",";

    D( std::cout<<"type : "<<type<<std::endl; )
    for( int i=0 ; i< type ; ++i) ost<<"0,";
    ost<<"1,";
    for( int i=0 ; i< 10-type ; ++i) ost<<"0,";
    ost<<"0";
    D( std::cout<<"4"<<std::endl; )
    return ost.str();
  }
}; 



void BookingTree(TTree* tree, Data& data) {
  tree->Branch("data",&data,"nLep/i:lep_q[4]/I:lep_pt[4]/F:lep_iso[4]/F:dilep_pt[4]/F:dilep_mass[4]/F:met/F:met_eta:met_phi:nJet/i:jet_pt[4]/F:jet_eta[4]/F:jet_phi[4]/F:jet_btag[4]/I:nbJet/i:type/I");
}


void putData(Data& data, int i){
  data.nLep=i;
  data.lep_q[0] = i;
  data.lep_pt[0]  = (float)i;
  data.lep_iso[0] = (float)i;
  data.lep_pt[1]  = (float)i;
  data.lep_iso[1] = (float)i;
  data.lep_pt[2] = (float)i;
  data.lep_iso[2]= (float)i;
  data.lep_pt[3]  = (float)i;
  data.lep_iso[3] = (float)i;
  data.jet_btag[1]= i%2;
  data.jet_btag[2] = i%2;
}


void PrintParameters(ofstream& outfile, TTree* tree)
{



}


void AnalyseEvents(ExRootTreeReader *treeReader, TTree* tree, int dataType, std::string csvOutFile )
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");

  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Data data;
  BookingTree(tree, data);

  int channel[9]={0,0,0,0,0,0,0,0,0};
  //allEntries = 1;
  //Int_t i, j, entry;
  ofstream os( csvOutFile, ios::out);
  for(int entry = 0; entry < allEntries; ++entry)
  {
    D( std::cout<<"Event : "<<entry<<std::endl; )
    data.reset();
    treeReader->ReadEntry(entry);
    channel[0]++;
    if ( dataType < 4 ) { 
      DecayChannel dc(branchParticle);
      data.type = dc.channel();
      D(std::cout<<"channel : "<<data.type<<std::endl;)
      channel[ dc.channel()+1]++;
    }


    // Lepton Part
    typedef std::pair<int, float> QI;
    typedef std::pair<TLorentzVector, QI > PQI;
    vector<PQI> leptons;
    for( int i=0 ; i< branchMuon->GetEntriesFast() ; ++i) {
      Muon* mu = (Muon*)( branchMuon->At(i));
      leptons.push_back( make_pair(mu->P4(), make_pair(mu->Charge, mu->IsolationVar)));
      D( std::cout<<"Put muon"<<std::endl; 
         std::cout<<"Muon iso: "<<mu->IsolationVar<<std::endl;
      )
    }
    for( int i=0 ; i< branchElectron->GetEntriesFast() ; ++i) {
      if ( i>=2 ) break;
      Electron* el = (Electron*)( branchElectron->At(i));
      leptons.push_back( make_pair(el->P4(), make_pair(el->Charge, el->IsolationVar)));
      D( std::cout<<"Put Electron"<<std::endl; 
         std::cout<<"Electron iso: "<<el->IsolationVar<<std::endl;
      ) 
    }
    std::sort(leptons.begin(), leptons.end(), [](PQI &left, PQI &right) {
      return left.first.Pt() < right.first.Pt();  }); 

    D( std::cout<<"Lepton size : "<< leptons.size()<<std::endl; )
    data.nLep = leptons.size(); 
    if ( leptons.size() >=2 ) channel[6]++;
    else if ( leptons.size() ==1 ) channel[7]++;
    else if ( leptons.size() ==0 ) channel[8]++;

    D( std::cout<<"LL"<<std::endl; )
    // If sorting is required, you need to FIX IT.
    //std::sort(leptons.begin(), leptons.end(), );
    for(int idx=0 ; idx< leptons.size() ; ++idx ) {
      if ( idx >= 4) break;
      data.lep_q[idx]   = leptons[idx].second.first;
      data.lep_iso[idx] = leptons[idx].second.second;
      data.lep_pt[idx]  = leptons[idx].first.Pt();
    }
    int dilepton_idx=0;
    if ( leptons.size() >2) {
      for(int i = 0 ; i < leptons.size()-1 ; ++i){
        for(int j = i+1; j< leptons.size() ; ++j ) {
          if( leptons[i].second.first * leptons[j].second.first != -1 ) continue;
          data.dilep_pt[dilepton_idx]   = (leptons[i].first+leptons[j].first).Pt();
          data.dilep_mass[dilepton_idx] = (leptons[i].first+leptons[j].first).M();
          dilepton_idx++;
        }
      }
    }
    for( int i=0 ; i< branchMissingET->GetEntriesFast() ; ++i) {
      D( if ( i>2 ) std::cout<<"Did you have 2nd MET? : size -> "<<branchMissingET->GetEntriesFast()<<std::endl;)
      MissingET* met = (MissingET*)( branchMissingET->At(i));
      data.met=met->MET;
      data.met_eta = met->Eta;
      data.met_phi = met->Phi;
    }
    data.nJet = branchJet->GetEntriesFast();
    for( int i=0 ; i< branchJet->GetEntriesFast() ; ++i) {
      if ( i>=4 ) break;
      Jet* jet = (Jet*)( branchJet->At(i));
      data.jet_pt[i] = jet->PT;
      data.jet_eta[i] = jet->Eta;
      data.jet_phi[i] = jet->Phi;
      data.jet_btag[i] = jet->BTag;
      if( jet->BTag ) data.nbJet++;
    }
     
    D(std::cout<<data.print()<<std::endl;)
    if ( !csvOutFile.empty() ) {
      D(std::cout<<"csv file will be written."<<std::endl;)
      os<<data.print()<<"\n";
    }
    tree->Fill();
  }
  os.close(); 

  if ( dataType < 4 ) {
    for( int i= 0 ; i < 9 ; i++) {
      std::cout<<"Channel : "<<i<<" "<<channel[i]<<"\t"<< (float)channel[i] / (float)channel[0]*100<<"%"<<std::endl;
    }
  }
}


int main(int argc, char* argv[])
{

  int dataType=0;
  std::string csvOutFile="";
  if ( argc != 3 && argc !=4 && argc !=5 ) {
    std::cerr<<"Wrong argument!"<<std::endl;
    exit(-1);
  }

  D(
  std::cout<<argv[0]<<std::endl;
  if ( argc ==4 ) std::cout<<argv[3]<<std::endl;
  if ( argc ==5 ) std::cout<<argv[4]<<std::endl;
  )
  std::string inputFile(argv[1]);
  std::string outFile(argv[2]);


  if ( argc == 4 ) dataType = std::atoi(argv[3]) ;
  if ( argc == 5 ) csvOutFile = std::string(argv[4]);

  D(std::cout<<"Input : "<<inputFile << "\tOut : "<<outFile<<std::endl;)

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile.c_str());

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  TFile* file = new TFile(outFile.c_str(), "RECREATE");
  TTree* tree = new TTree("delphes","delphes");
  AnalyseEvents(treeReader, tree, dataType, csvOutFile );

  file->Write();
  file->Close();


  cout << "** Exiting..." << endl;

}

