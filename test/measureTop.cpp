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
using namespace std;

class ExRootTreeReader;

class JetTree {
public :
  float weight;
  float jet_pt[2]; 
  float jet_eta[2]; 
  float jet_phi[2];
  float top_mass[2]; 
  float ttbar_mass;
  int bjet_charge[2];
  int bjet_partonPdgId[2]; 
  int lep_charge[2];
  int bjet_btag[2];
  int bjet_nCharged[2];
  JetTree() {
    init();
  }
  void reset() { init(); }
  void init() {
    weight = -1e5;
    for (unsigned int i=0 ; i< 2 ; ++i) {
      jet_pt[i]=0.f;
      jet_eta[i]=-9.f;
      jet_phi[i]=-9.f;
      top_mass[i]=0.f;
      bjet_charge[i]=-9;
      bjet_partonPdgId[i] = 0 ;
      bjet_btag[i]= 0;
      bjet_charge[i]=-9;
      bjet_nCharged[i]=0;
    }
    ttbar_mass=0.f;
  }
};

class OutFileClass {
private :
  TFile* file_;
  std::map<std::string, TTree*> tree_;
  std::map<std::string, TH1*> th1_;
  std::map<std::string, TH2*> th2_;
public :
  TRandom3* random;
  OutFileClass(TFile* file) {
    SetFile(file);
    random = new TRandom3();
  }
  JetTree data, data1;
  void SetFile(TFile* file) { file_ = file; }
  TFile* GetFile() { return file_;}
  void AddTH1(TH1* h1) { th1_.insert( pair<std::string, TH1*>(std::string(h1->GetName()), h1)) ; }
  void AddTH1(std::string str, TH1* h1) { th1_.insert( pair<std::string, TH1*>(str, h1)) ; }
  void AddTree(TTree* tree) { AddTree(tree->GetName(), tree); }
  void AddTree(std::string str, TTree* tree) { tree_.insert( pair<std::string, TTree*>(str, tree)) ; }
  TH1* GetTH1(std::string str) {
    auto value = th1_.find(str);
    if ( value != th1_.end() ) return value->second;
    else return nullptr;
  }
  TTree* GetTree(std::string str) {
    auto value = tree_.find(str);
    if ( value != tree_.end() ) return value->second;
    else return nullptr;
  }


  void AllHistWrite() { 
    for( auto i = th1_.begin() ; i != th1_.end(); ++i) {
      i->second->Write();
    } 
    for( auto i = th2_.begin() ; i != th2_.end(); ++i) {
      i->second->Write();
    }
  } 
  void AllTreeWrite() {
    for( auto i = tree_.begin() ; i != tree_.end(); ++i) {
      i->second->Write();
    } 
  }

};

typedef std::pair<int, LV> QLV;
typedef std::pair<int, Jet*> QJET;

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
        D(std::cout<<std::endl;
        std::cout<<"Jet PT : "<<jet->PT<<" Jet eta: "<<jet->Eta<<" Jet phi : "<<jet->Phi<<std::endl;)
      }
      D(  
      std::cout<<trackIdx<<"th's track charge : "<<track->Charge<<std::endl;
      std::cout<<"track PT : "<<track->PT<<" track Eta : "<<track->Eta<<"  track phi : "<<jet->Phi<<std::endl;
      )
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
  return charge/weightSum*100;
}


//------------------------------------------------------------------------------


std::vector<QLV> selectedLepton( TClonesArray* muons, TClonesArray* electrons) {
  std::vector<QLV> selLeptons;
  for( int i= 0 ; i< muons->GetEntriesFast() ; ++i) {
    Muon* muon = (Muon*) muons->At(i);
    if ( muon->PT<20 || abs(muon->Eta)>2.5) continue;
    selLeptons.push_back( make_pair(muon->Charge, common::TLVtoLV(muon->P4()) ));
  }
  for( int i= 0 ; i< electrons->GetEntriesFast() ; ++i) {
    Electron* electron = (Electron*) electrons->At(i);
    if ( electron->PT<20 || abs(electron->Eta)>2.5) continue;
    selLeptons.push_back( make_pair(electron->Charge, common::TLVtoLV(electron->P4()) ));
  }
  return selLeptons; 
}
std::vector<QJET> selectedJet( TClonesArray* jets ) {
  std::vector<QJET> selJets;
  for(int i = 0; i < jets->GetEntriesFast(); ++i)
  {
    Jet* jet = (Jet*) jets->At(i);
    if ( jet->PT < 30 || abs(jet->Eta)>2.5 ) continue;
    selJets.push_back( make_pair( getJetCharge(jet,true), jet) );
  }
  return selJets; 
}

int nBJet( TClonesArray* branchParticle, std::vector<QJET> jets) {
  int nbjet=0;
  for( unsigned int i=0 ; i< jets.size() ; ++i) {
    Jet* jet = jets[i].second;
    if ( jet->Constituents.At(0)->IsA() == GenParticle::Class()) { 
      int parton_idx = DecayChannel::FindJetParton( branchParticle, jet);
      int parton =0 ; 
      if ( parton_idx != -1 ) parton = ((GenParticle*)branchParticle->At(parton_idx))->PID;
      if ( abs(parton)==5) {
        nbjet++;
      }
    } 
    else if ( jet->BTag >0 ) nbjet++;
  }
  return nbjet;
}


//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader, OutFileClass& ofc, const char* jetType)
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

  TClonesArray *branchMET = treeReader->UseBranch("MissingET");
  
  
  TClonesArray *branchJet = treeReader->UseBranch(jetType);

  Long64_t allEntries = treeReader->GetEntries();

  D(cout << "** Chain contains " << allEntries << " events" << endl;)

  GenParticle *particle;

  //Int_t i, j, entry;
  int percent=0;
  TH1* nevt = ofc.GetTH1((TString("nEvent_")+jetType).Data());
  for(int entry = 0; entry < allEntries; ++entry) {
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if (entry % (allEntries/100)==0 ) std::cout<<"Event : "<<entry<<"("<<percent++<<"%)"<<std::endl;
    nevt->Fill(1);

    /*
    // Event Selection using GenParticles.
    DecayChannel dc(branchParticle);
    if ( dc.channel() !=0 ) continue; // Only choose dilepton channel except tau decay.
    */
    //std::cout<<"Pass GEN EventSelector cut"<<std::endl;
    nevt->Fill(2);

    // Least POG Numbers, // Simple event Selection.
    //D( std::cout<<"This is dilepton events!"<<std::endl; )
    std::vector<QLV> selLepton = selectedLepton( branchMuon, branchElectron);
    if ( selLepton.size()<2  || selLepton[0].first*selLepton[1].first != -1 ) continue;
    //std::cout<<"Pass Lepton cut"<<std::endl;
    nevt->Fill(3);
    if ( abs((selLepton[0].second+selLepton[1].second).mass()-91.)<15) continue;
    nevt->Fill(4);
    std::vector<QJET> selJet = selectedJet(branchJet);
    if ( selJet.size()<2 ) continue;
    //std::cout<<"Pass Jet cut"<<std::endl;
    nevt->Fill(5);
    MissingET* met = (MissingET*) branchMET->At(0);
    if ( met->MET<40 ) continue;
    //std::cout<<"Pass MET cut"<<std::endl;
    nevt->Fill(6);
    ofc.GetTH1((TString("nJet_")+jetType).Data())->Fill( selJet.size() );

    if ( nBJet(branchParticle, selJet) <1) continue;
    nevt->Fill(7);

    ofc.data.reset();
    ofc.data1.reset();
    // Calculate Weight!
    float before_weight = -1.0;
    float before_weight_charged = -1.0;
    cat::DESYSmearedSolver* cmskin_solver = new cat::DESYSmearedSolver();
    //cat::CMSKinSolver* cmskin_solver = new cat::CMSKinSolver();
    for( unsigned int i = 0 ; i < selJet.size() ; ++i) {
      for( unsigned int j = 0 ; j < selJet.size() ; ++j) { 
        if (i==j) continue;
 
        LV jet1 = common::TLVtoLV(selJet[i].second->P4());
        LV jet2 = common::TLVtoLV(selJet[j].second->P4()); 
        const LV POG[] = { common::TLVtoLV(met->P4()), selLepton[0].second, selLepton[1].second, jet1, jet2}; 
        cmskin_solver->solve(POG);
        double quality = cmskin_solver->quality(); 
        if (quality<0) continue;
        int jet1_parton_idx = DecayChannel::FindJetParton( branchParticle, selJet[i].second);
        int jet2_parton_idx = DecayChannel::FindJetParton( branchParticle, selJet[j].second);

        if ( before_weight < quality ) { 
          before_weight = quality;

          ofc.data.weight = (float)quality;

          ofc.data.jet_pt[0] = jet1.pt();
          ofc.data.jet_eta[0] = jet1.eta();
          ofc.data.jet_phi[0] = jet1.phi();
          ofc.data.bjet_charge[0] = selJet[i].first;
          ofc.data.bjet_nCharged[0] = selJet[i].second->NCharged;
          if ( jet1_parton_idx != -1 ) ofc.data.bjet_partonPdgId[0] = ((GenParticle*)branchParticle->At(jet1_parton_idx))->PID; 
          else ofc.data.bjet_partonPdgId[0] = 0;
          ofc.data.lep_charge[0] = selLepton[0].first;
          ofc.data.bjet_btag[0] = selJet[i].second->BTag;
          
          ofc.data.jet_pt[1] = jet2.pt();
          ofc.data.jet_eta[1] = jet2.eta();
          ofc.data.jet_phi[1] = jet2.phi();
          ofc.data.bjet_charge[1] = selJet[j].first;
          ofc.data.bjet_nCharged[1] = selJet[j].second->NCharged;
          if ( jet2_parton_idx != -1 ) ofc.data.bjet_partonPdgId[1] = ((GenParticle*)branchParticle->At(jet2_parton_idx))->PID; 
          else ofc.data.bjet_partonPdgId[1] = 0;
          ofc.data.lep_charge[1] = selLepton[1].first;
          ofc.data.bjet_btag[1] = selJet[i].second->BTag;

          if ( selJet[i].second->Constituents.At(0)->IsA() == GenParticle::Class() && abs(ofc.data.bjet_partonPdgId[0])==5 ) { 
            ofc.data.bjet_btag[0] = 1;  
          }
          if ( selJet[j].second->Constituents.At(0)->IsA() == GenParticle::Class() && abs(ofc.data.bjet_partonPdgId[1])==5 ) { 
            ofc.data.bjet_btag[1] = 1; 
          } 
          ofc.data.top_mass[0]   = (float)cmskin_solver->t1().mass();
          ofc.data.top_mass[1]   = (float)cmskin_solver->t2().mass();
          ofc.data.ttbar_mass  = (float)(cmskin_solver->t1()+cmskin_solver->t2()).mass();
        }

        if ( (selLepton[0].first * selJet[i].first+selLepton[1].first* selJet[j].first)  < -50 ) {
          if ( before_weight_charged > quality ) continue; 
          before_weight_charged = quality;

          ofc.data1.weight = (float)quality;

          ofc.data1.jet_pt[0] = jet1.pt();
          ofc.data1.jet_eta[0] = jet1.eta();
          ofc.data1.jet_phi[0] = jet1.phi();
          ofc.data1.bjet_charge[0] = selJet[i].first;
          ofc.data1.bjet_nCharged[0] = selJet[i].second->NCharged;

       
          if ( jet1_parton_idx != -1 ) ofc.data1.bjet_partonPdgId[0] = ((GenParticle*)branchParticle->At(jet1_parton_idx))->PID; 
          else ofc.data1.bjet_partonPdgId[0] = 0;
          ofc.data1.lep_charge[0] = selLepton[0].first;
          ofc.data1.bjet_btag[0] = selJet[i].second->BTag;
          
          ofc.data1.jet_pt[1] = jet2.pt();
          ofc.data1.jet_eta[1] = jet2.eta();
          ofc.data1.jet_phi[1] = jet2.phi();
          ofc.data1.bjet_charge[1] = selJet[j].first;
          ofc.data1.bjet_nCharged[1] = selJet[j].second->NCharged;
          if ( jet2_parton_idx != -1 ) ofc.data1.bjet_partonPdgId[1] = ((GenParticle*)branchParticle->At(jet2_parton_idx))->PID; 
          else ofc.data1.bjet_partonPdgId[1] = 0;
          ofc.data1.lep_charge[1] = selLepton[1].first;
          ofc.data1.bjet_btag[1] = selJet[i].second->BTag;
          if ( selJet[j].second->Constituents.At(0)->IsA() == GenParticle::Class() && abs(ofc.data1.bjet_partonPdgId[1])==5 ) ofc.data1.bjet_btag[1] = 1;  
          ofc.data1.top_mass[0]   = (float)cmskin_solver->t1().mass();
          ofc.data1.top_mass[1]   = (float)cmskin_solver->t2().mass();
          ofc.data1.ttbar_mass = (float)(cmskin_solver->t1()+cmskin_solver->t2()).mass();
        } 
      }
    }
    if ( ofc.data.weight > 0 ) { 
      std::cout<<"Weight : "<<ofc.data.weight<<std::endl;
      std::cout<<"top1 mass : "<<ofc.data.top_mass[0]<<"  top2 mass : "<<ofc.data.top_mass[1]<<std::endl;
      ofc.GetTree( (TString(jetType)+"Tree").Data() )->Fill();
    }
    if ( ofc.data1.weight >0 ) { 
      std::cout<<"Weight with Charged info : "<<ofc.data1.weight<<std::endl;
      std::cout<<"top1 mass : "<<ofc.data1.top_mass[0]<<"  top2 mass : "<<ofc.data1.top_mass[1]<<std::endl;
      ofc.GetTree( (TString(jetType)+"TreeCharge").Data() )->Fill();
    }
  }
}

//------------------------------------------------------------------------------

void Write(OutFileClass& ofc)
{
  ofc.GetFile()->cd();
  //result->Print("png");
  ofc.AllHistWrite();
  ofc.AllTreeWrite();
}

//------------------------------------------------------------------------------


void BookingTree(OutFileClass& ofc, std::string treeName, const char* treeTypes) {
  ofc.GetFile()->cd();
  TTree* tree = new TTree(treeName.c_str(), treeName.c_str());
  tree->Branch(treeName.c_str(), &ofc.data, treeTypes );
  ofc.AddTree(tree);
}

void BookingHist(OutFileClass& ofc)
{
  ofc.GetFile()->cd();
  TH1F* h1 = new TH1F("nEvent_Jet","Number of Events",7,1,8);
  auto xaxis = h1->GetXaxis();
  xaxis->SetBinLabel(1,"Total");
  xaxis->SetBinLabel(2,"GenEvt");
  xaxis->SetBinLabel(3,"nlep2");
  xaxis->SetBinLabel(4,"ZVeto");
  xaxis->SetBinLabel(5,"nJet2");
  xaxis->SetBinLabel(6,"MET40");
  xaxis->SetBinLabel(7,"nbJet1");
  ofc.AddTH1(h1);

  TH1F* h2 = new TH1F("nJet_Jet","Number of Jets",10,0,10);
  ofc.AddTH1(h2);

  TH1F* h3 = new TH1F("nEvent_GenJet","Number of Events",7,1,8);
  xaxis = h3->GetXaxis();
  xaxis->SetBinLabel(1,"Total");
  xaxis->SetBinLabel(2,"GenEvt");
  xaxis->SetBinLabel(3,"nlep2");
  xaxis->SetBinLabel(4,"ZVeto");
  xaxis->SetBinLabel(5,"nJet2");
  xaxis->SetBinLabel(6,"MET40");
  xaxis->SetBinLabel(7,"nbJet1");
  ofc.AddTH1(h3);

  TH1F* h4 = new TH1F("nJet_GenJet","Number of Jets",10,0,10);
  ofc.AddTH1(h4);
}
  

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
  TFile* file = TFile::Open(outputFile,"RECREATE"); 
  auto ofc = OutFileClass(file);
  BookingHist(ofc);
  TString treeBranch("quality/F:jet_pt[2]/F:jet_eta[2]/F:jet_phi[2]/F:top_mass[2]/F:ttbar_mass/F:bjet_charge[2]/I:bjet_partonPdgId[2]/I:lep_charge[2]/I:bjet_btag[2]/I:bjet_nCharged[2]/I");
  BookingTree(ofc,"JetTree",       treeBranch.Data());
  BookingTree(ofc,"JetTreeCharge", treeBranch.Data());
  BookingTree(ofc,"GenJetTree",       treeBranch.Data());
  BookingTree(ofc,"GenJetTreeCharge", treeBranch.Data());
  AnalyseEvents(treeReader,ofc, "GenJet");
  AnalyseEvents(treeReader,ofc, "Jet");

  Write(ofc);
  cout << "** Exiting..." << endl;

  file->Close();
  //delete result;
  delete treeReader;
  delete chain;

}
