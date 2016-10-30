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
#include <TLorentzVector.h>
using namespace std;

class ExRootTreeReader;

typedef std::pair<int, LV> QLV;
typedef std::pair<int, Jet*> QJET;
typedef TLorentzVector TLV;



float d0_massWindow = 0.05;
float dstar_massWindow = 0.01;
TLorentzVector TracktoTLV( Track* track ) {
    TLV track_TLV = track->P4();
    double track_mass =0.;
    if      ( abs(track->PID)==11)  track_mass = 0.000511169;
    else if ( abs(track->PID)==13)  track_mass = 0.105652;
    else if ( abs(track->PID)==211) track_mass = 0.139526;
    else if ( abs(track->PID)==321) track_mass = 0.493652;
    TLorentzVector newTrack;
    newTrack.SetPtEtaPhiM( track_TLV.Pt(), track_TLV.Eta(), track_TLV.Phi(), track_mass);
    return newTrack;
}



class SecVtx  { 
public :
  SecVtx() {
    nDau_=0;
  }
  SecVtx(Track* trackA, Track* trackB, TClonesArray* branchParticle=nullptr) {

    float track1_mass = 0, track2_mass = 0;
    
    auto trackA_TLV = TracktoTLV(trackA);
    auto trackB_TLV = TracktoTLV(trackB);
    pos_ = trackA_TLV + trackB_TLV;

    dau_pid_[0] = trackA->PID;
    dau_pid_[1] = trackB->PID;
    charge_ = (trackA->Charge + trackB->Charge);
    daus_[0] = trackA_TLV;
    daus_[1] = trackB_TLV;

    if ( branchParticle != nullptr) {
      auto trackA_GEN = (GenParticle*)trackA->Particle.GetObject();
      GenParticle* trackAmother = (GenParticle*)branchParticle->At(trackA_GEN->M1);
      auto trackB_GEN = (GenParticle*)trackB->Particle.GetObject();
      GenParticle* trackBmother = (GenParticle*)branchParticle->At(trackB_GEN->M1);
      vx_[0] = trackAmother->X;
      vx_[1] = trackBmother->X;
      vy_[0] = trackAmother->Y;
      vy_[1] = trackBmother->Y;
      vz_[0] = trackAmother->Z;
      vz_[1] = trackBmother->Z;
    }
    nDau_ = 2;
    isFromTop_= 0;
  }
  SecVtx(SecVtx* secvtxA, Track* trackB, TClonesArray* branchParticle=nullptr) {
    float track2_mass = 0;
    auto trackB_TLV = TracktoTLV(trackB);
    pos_ = secvtxA->P4() + trackB_TLV;
    charge_ = trackB->PID/abs(trackB->PID);
    dau_pid_[0] = secvtxA->dau_pid(0);
    dau_pid_[1] = secvtxA->dau_pid(1);
    dau_pid_[2] = trackB->PID;
    daus_[0] = secvtxA->dau(0);
    daus_[1] = secvtxA->dau(1);
    daus_[2] = trackB_TLV;

    if ( branchParticle != nullptr) {
      vx_[0] = secvtxA->vx(0);
      vx_[1] = secvtxA->vx(1);
      vy_[0] = secvtxA->vy(0);
      vy_[1] = secvtxA->vy(1);
      vz_[0] = secvtxA->vz(0);
      vz_[1] = secvtxA->vz(1);
      auto trackB_GEN = (GenParticle*)trackB->Particle.GetObject();
      GenParticle* trackBmother = (GenParticle*)branchParticle->At(trackB_GEN->M1);
      vx_[2] = trackBmother->X;
      vy_[2] = trackBmother->Y;
      vz_[2] = trackBmother->Z;
    }
    nDau_ = 3;
    isFromTop_= 0;
  }
  TLV P4() { return pos_; } 
  Float_t mass() { return pos_.M(); }
  Float_t pt() { return pos_.Pt(); }
  Float_t eta() { return pos_.Eta(); }
  Float_t phi() { return pos_.Phi(); }
  Int_t charge() { return charge_; }
  Int_t dau_pid(Int_t idx) { return dau_pid_[idx]; }
  Int_t softLep() { return nLep_; }
  TLV dau(Int_t idx) { 
    if ( idx>3) { std::cout<<"Wrong idx. It is too high."<<std::endl; return TLorentzVector(); }
    //if ( idx > nDau_-1 ) { std::cout<<"Wrong idx. It is null ptr."<<std::endl; return TLorentzVector(); } 
    return daus_[idx]; 
  }
  Float_t vx() {
    if ( nDau_==2) {
      float x = (vx_[0]+vx_[1])/2.0;
      return x;
    }
    else if (nDau_==3) {
      float x = (vx_[0]+vx_[1]+vx_[2])/3.0;
      return x;
    }
    else return 0.0f;  
  }
  Int_t pid() {
    return pid_;
  } 
  Float_t vy() {
    if ( nDau_==2) {
      float y = (vy_[0]+vy_[1])/2.0;
      return y;
    }
    else if (nDau_==3) {
      float y = (vy_[0]+vy_[1]+vy_[2])/3.0;
      return y;
    }
    else return 0.0f;  
  } 
  Float_t vz() {
    if ( nDau_==2) {
      float z = (vz_[0]+vz_[1])/2.0;
      return z;
    }
    else if (nDau_==3) {
      float z = (vz_[0]+vz_[1]+vz_[2])/3.0;
      return z;
    }
    else return 0.0f;  
  } 
  Float_t vx(Int_t idx) {
    if ( idx< nDau_) { return vx_[idx]; }
    else { std::cout<<"Wrong vertex position. It is null ptr"<<std::endl; return 0.f; }
  }
  Float_t vy(Int_t idx) {
    if ( idx< nDau_) { return vy_[idx]; }
    else { std::cout<<"Wrong vertex position. It is null ptr"<<std::endl; return 0.f; }
  }
  Float_t vz(Int_t idx) {
    if ( idx< nDau_) { return vz_[idx]; }
    else { std::cout<<"Wrong vertex position. It is null ptr"<<std::endl; return 0.f; }
  }
  Float_t L3D() {
      return TMath::Sqrt(vx()*vx()+vy()*vy()+vz()*vz());
  }
  Float_t LXY() {
      return TMath::Sqrt(vx()*vx()+vy()*vy());
  }

  Float_t vxd() {
    if ( nDau_ ==2 ) {
      float vxd_temp = vx_[0]-vx_[1];
      return abs( vxd_temp);
    }
    else if ( nDau_ ==3 ) {
      float vxd_temp1 = vx_[0]-vx_[1];
      float vxd_temp2 = vx_[0]-vx_[2];
      float vxd_temp3 = vx_[1]-vx_[2];
      float vxd__ = (abs(vxd_temp1) + abs(vxd_temp2) + abs(vxd_temp3)) /3.;
      return vxd__;
    }
    else return 0; 
  }
  Float_t vyd() {
    if ( nDau_ ==2 ) {
      float vyd_temp = vy_[0]-vy_[1];
      return abs(vyd_temp) ;
    }
    else if ( nDau_ ==3 ) {
      float vyd_temp1 = vy_[0]-vy_[1];
      float vyd_temp2 = vy_[0]-vy_[2];
      float vyd_temp3 = vy_[1]-vy_[2];
      float vyd__ = (abs(vyd_temp1) + abs(vyd_temp2) + abs(vyd_temp3)) /3.;
      return vyd__;
    }
    else return 0; 
  }
  Float_t vzd() {
    if ( nDau_ ==2 ) {
      float vzd_temp = vz_[0]-vz_[1];
      return TMath::Sqrt( vzd_temp*vzd_temp) ;
    }
    else if ( nDau_ ==3 ) {
      float vzd_temp1 = vz_[0]-vz_[1];
      float vzd_temp2 = vz_[0]-vz_[2];
      float vzd_temp3 = vz_[1]-vz_[2];
      float vzd__ = (abs(vzd_temp1) + abs(vzd_temp2) + abs(vzd_temp3)) /3.;
      return vzd__;
    }
    else return 0; 
  }
  void setDRPT(DecayChannel dc, TClonesArray* branchParticle, int pdgId) {
    pid_ = pdgId;
    int mcTruth = dc.SearchParticle(branchParticle, pdgId, pos_ );
    if ( mcTruth == -1 ) return;
    int isFromTop =0;
    int isFromTopIdx = dc.isFromTop( branchParticle, mcTruth );
    if ( isFromTopIdx != -1 ) isFromTop_ = ((GenParticle*)branchParticle->At( isFromTopIdx))->PID;

    auto gen = (GenParticle*)branchParticle->At( mcTruth );
    dRTrue_ = gen->P4().DeltaR( pos_ );
    delPtTrue_ = abs(gen->PT-pos_.Pt())/gen->PT;

  }
  void setSoftLepton(int nLep) {nLep_ = nLep;}


  float DR() { return dRTrue_; }
  float DPT() { return delPtTrue_; }

  Int_t isFromTop() { return isFromTop_;}
  Int_t nDau() { return nDau_; }
  void setD0mass( float d0mass) {
    diffmass_ = pos_.M() - d0mass;
    d0mass_ = d0mass;
  }
  Float_t d0mass() { return d0mass_;}
  Float_t diffmass() { return diffmass_;}
private :
  Int_t charge_, pid_;
  TLV pos_;
  TLV daus_[3];
  Int_t nDau_;
  Int_t dau_pid_[3];
  Float_t vx_[3];
  Float_t vy_[3];
  Float_t vz_[3];
  Float_t diffmass_, d0mass_;
  Int_t isFromTop_;
  Float_t delPtTrue_;
  Float_t dRTrue_;
  Int_t nLep_;
};    

std::pair< TLorentzVector, TLorentzVector> lsvPairing( QLV* lep1 , QLV* lep2, SecVtx* sv) {
  int lep1_pid = lep1->first;
  int lep2_pid = lep2->first;
  TLorentzVector lep1_TLV = common::LVtoTLV(lep1->second);
  TLorentzVector lep2_TLV = common::LVtoTLV(lep2->second);

  TLorentzVector lsv1, lsv2;
  lsv1 = lep1_TLV+sv->P4();
  lsv2 = lep2_TLV+sv->P4();

  if ( sv->isFromTop() == 6) {
    // top to positron
    if ( lep1_pid <0 ) return make_pair(lsv1, lsv2);
    else               return make_pair(lsv2, lsv1);
  }
  else if ( sv->isFromTop() == -6) {
    // atop to electron
    if ( lep1_pid >0 ) return make_pair(lsv1, lsv2);
    else               return make_pair(lsv2, lsv1);
  }
  // if sv can not find top!, lsv1 is a lower inv mass set. 
  else {
    if (  lsv1.M() < lsv2.M() ) return make_pair(lsv1, lsv2); 
    else                        return make_pair(lsv2, lsv1); 
  }
}


class saveData {
public :
  //Int_t dataType;
  Int_t lep_charge[2], lep_pid[2];
  Float_t lep_pt[2], lep_eta[2], lep_phi[2], lep_mass[2];

  Int_t sv_pid, sv_charge, sv_isFromTop;
  Float_t sv_pt, sv_eta, sv_phi, sv_mass, sv_diffmass, sv_d0mass;
  Float_t sv_vx, sv_vy,  sv_vz;
  Float_t sv_dx, sv_dy,  sv_dz;
  Float_t sv_L3D, sv_LXY;
  Int_t sv_softlep;
  Float_t sv_dRTrue, sv_delPtTrue;
  

  Int_t sv_dau_pid[3];
  Float_t sv_dau_pt[3], sv_dau_eta[3], sv_dau_phi[3], sv_dau_mass[3];
  Int_t lsv_pid[2];
  Float_t lsv_pt[2], lsv_eta[2], lsv_phi[2], lsv_mass[2],lsv_diffmass[2], lsv_d0mass[2];
 
  void reset() {
    //dataType=0;
    sv_pid = 0 ;
    sv_charge =0; 
    sv_isFromTop=-9; 
    sv_pt = -9.0; 
    sv_eta = -9.0; 
    sv_phi = -9.0; 
    sv_mass = -9.0;
    sv_diffmass = -9.0;
    sv_d0mass = -9.0;
    sv_vx = 0.f;
    sv_vy = 0.f;
    sv_vz = 0.f;
    sv_dx = 0.f;
    sv_dy = 0.f;
    sv_dz = 0.f;
    sv_LXY= 0.f;
    sv_L3D= 0.f;
    sv_softlep=0;
    sv_dRTrue = 999.f;
    sv_delPtTrue = 1.f;
 
    for( int idx = 0; idx< 3 ; idx++) { 
      sv_dau_pid[idx] =0;  
      sv_dau_pt[idx] = -9.0; 
      sv_dau_eta[idx] = -9.0; 
      sv_dau_phi[idx] = -9.0; 
      sv_dau_mass[idx] = -9.0; 
    }
    for( int idx = 0; idx< 2 ; idx++) { 
 
      lep_charge[idx] =0;  
      lep_pid[idx] =0 ;
      lep_pt[idx] = -9.0; 
      lep_eta[idx] = -9.0; 
      lep_phi[idx] = -9.0; 
      lep_mass[idx] = -9.0; 

      
      lsv_pid[idx] = 0;
      lsv_pt[idx] = -9.0;
      lsv_eta[idx] = -9.0;
      lsv_phi[idx] = -9.0;
      lsv_mass[idx] = -9.0;
      lsv_diffmass[idx] = 0.0;
      lsv_d0mass[idx] = -9.0;
 
    }
  }
 void fillLSV(QLV* lep1_QLV, QLV* lep2_QLV, SecVtx* svx){
      std::pair<TLorentzVector,TLorentzVector> pair = lsvPairing(lep1_QLV, lep2_QLV, svx);
      TLorentzVector lsvx1 = pair.first;
      TLorentzVector lsvx2 = pair.second;

      lsv_pid[0]  = svx->pid();
      lsv_pt[0]   = lsvx1.Pt();
      lsv_eta[0]  = lsvx1.Eta();
      lsv_phi[0]  = lsvx1.Phi();
      lsv_mass[0] = lsvx1.M();
      lsv_d0mass[0] = svx->d0mass();
      lsv_diffmass[0] = svx->diffmass();

      lsv_pid[1]  = svx->pid();
      lsv_pt[1] = lsvx2.Pt();
      lsv_eta[1] = lsvx2.Eta();
      lsv_phi[1] = lsvx2.Phi();
      lsv_mass[1] = lsvx2.M();
      lsv_d0mass[1] = svx->d0mass();
      lsv_diffmass[1] = svx->diffmass();
  }
  void init(QLV* lep1_QLV, QLV* lep2_QLV, SecVtx* svx) {
    reset();
    LV* lep1 = &(lep1_QLV->second);
    LV* lep2 = &(lep2_QLV->second);
    int lep1_pid = lep1_QLV->first;
    int lep2_pid = lep2_QLV->first;
    sv_pid = svx->pid();
    sv_charge = svx->charge();
    sv_isFromTop = svx->isFromTop();
    sv_pt = svx->pt();
    sv_eta = svx->eta();
    sv_phi = svx->phi();
    sv_mass = svx->mass();
    sv_vx = svx->vx();
    sv_vy = svx->vy();
    sv_vz = svx->vz();
    sv_dx = svx->vxd();
    sv_dy = svx->vyd();
    sv_dz = svx->vzd();
    sv_LXY = svx->LXY();
    std::cout<<"X : "<<svx->vx()<<" LXY : "<<sv_LXY<<std::endl;
    sv_L3D = svx->L3D();
    sv_softlep = svx->softLep();
    sv_dRTrue = svx->DR();
    sv_delPtTrue = svx->DPT();

    for ( int idx =0 ; idx<svx->nDau() ; idx++) {
      sv_dau_pid[idx]   = svx->dau_pid(idx);
      std::cout<<"Daughter's pid : "<<sv_dau_pid[idx]<<std::endl;
      sv_dau_pt[idx]   = svx->dau(idx).Pt();
      sv_dau_eta[idx]  = svx->dau(idx).Eta();
      sv_dau_phi[idx]  = svx->dau(idx).Phi();
      sv_dau_mass[idx] = svx->dau(idx).M();
    }

    lep_charge[0] = -lep1_pid/abs(lep1_pid);
    lep_pid[0] = lep1_pid;
    lep_pt[0] = lep1->pt();
    lep_eta[0] = lep1->eta();
    lep_phi[0] = lep1->phi();
    lep_mass[0] = lep1->mass();
    std::cout<<"lep_mass[0] : "<<lep_mass[0]<<std::endl;

    lep_charge[1] = -lep2_pid/abs(lep2_pid);
    lep_pid[1] = lep2_pid;
    lep_pt[1] = lep2->pt();
    lep_eta[1] = lep2->eta();
    lep_phi[1] = lep2->phi();
    lep_mass[1] = lep2->mass();
    
    fillLSV(&lep1_QLV[0], &lep2_QLV[1], svx);
  }
};

class OutFileClass {
private :
  TFile* file_;
  std::map<std::string, TTree*> tree_;
  std::map<std::string, TH1*> th1_;
  std::map<std::string, TH2*> th2_;
public :
  OutFileClass(TFile* file) {
    SetFile(file);
  }
  saveData data;
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
    else { std::cout<<"Can not find tree!"<<std::endl; return nullptr; }
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


//------------------------------------------------------------------------------


std::vector<QLV> selectedLepton( TClonesArray* muons, TClonesArray* electrons) {
  std::vector<QLV> selLeptons;
  for( int i= 0 ; i< muons->GetEntriesFast() ; ++i) {
    Muon* muon = (Muon*) muons->At(i);
    if ( muon->PT<20 || abs(muon->Eta)>2.4) continue;
    // pid == 13, Q = -1
    int pid = muon->Charge* -13;
    TLorentzVector muonP4 = muon->P4();
    TLorentzVector newMuonP4 = TLorentzVector();
    newMuonP4.SetPtEtaPhiM( muonP4.Pt(), muonP4.Eta(), muonP4.Phi(), 0.105652);
    selLeptons.push_back( make_pair(pid, common::TLVtoLV( newMuonP4 ) ));
  }
  for( int i= 0 ; i< electrons->GetEntriesFast() ; ++i) {
    Electron* electron = (Electron*) electrons->At(i);
    if ( electron->PT<20 || abs(electron->Eta)>2.4) continue;
    // pid : 11 = Q:-1 * -11 
    int pid = electron->Charge* -11;
    TLorentzVector elecP4 = electron->P4();
    TLorentzVector newElecP4 = TLorentzVector();
    newElecP4.SetPtEtaPhiM( elecP4.Pt(), elecP4.Eta(), elecP4.Phi(), 0.000511169 );
    selLeptons.push_back( make_pair(pid, common::TLVtoLV( newElecP4) ));
  }
  std::sort( selLeptons.begin(), selLeptons.end(), []( QLV a, QLV b) { return a.second.pt() > b.second.pt(); });
  return selLeptons; 
}
std::vector<QJET> selectedJet( TClonesArray* jets ) {
  std::vector<QJET> selJets;
  for(int i = 0; i < jets->GetEntriesFast(); ++i)
  {
    Jet* jet = (Jet*) jets->At(i);
    if ( jet->PT < 30 || abs(jet->Eta)>2.4 ) continue;
    selJets.push_back( make_pair(jet->Charge, jet) );
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
void AnalyseEvents(ExRootTreeReader *treeReader, OutFileClass& ofc)
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
  
  
  TClonesArray *branchJet = treeReader->UseBranch("Jet");

  Long64_t allEntries = treeReader->GetEntries();

  D(cout << "** Chain contains " << allEntries << " events" << endl;)

  GenParticle *particle;

  //Int_t i, j, entry;
  int percent=0;
  TH1* nevt = ofc.GetTH1("nEvent");
  Int_t nJpsi= 0 , nD0=0, nDstar=0;
  for(int entry = 0; entry < allEntries; ++entry)
  {
    //Reset tree contents at event start!
    ofc.data.reset();
    treeReader->ReadEntry(entry);
    // Loop over all jets in event
    if (entry % (allEntries/100)==0 ) std::cout<<"Event : "<<entry<<"("<<percent++<<"%)"<<std::endl;
    nevt->Fill(1);

    // Event Selection using GenParticles.
    DecayChannel dc(branchParticle);
    if ( dc.channel() ==0 ) nevt->Fill(2); // Only choose dilepton channel except tau decay.
    //std::cout<<"Pass GEN EventSelector cut"<<std::endl;
    //nevt->Fill(2);
    int dcChannel = dc.channel();

    // Least POG Numbers, // Simple event Selection.
    //D( std::cout<<"This is dilepton events!"<<std::endl; )
    std::vector<QLV> selLepton = selectedLepton( branchMuon, branchElectron);
    if ( selLepton.size()<2  || selLepton[0].first*selLepton[1].first > 0 || (selLepton[0].second+selLepton[1].second).M()<20  ) continue;
    //std::cout<<"Pass Lepton cut"<<std::endl;
    nevt->Fill(3);
    if ( ((selLepton[0].second+selLepton[1].second).mass()-91.)<15) continue;
    nevt->Fill(4);
    std::vector<QJET> selJet = selectedJet(branchJet);
    if ( selJet.size()<2 ) continue;
    //std::cout<<"Pass Jet cut"<<std::endl;
    nevt->Fill(5);
    MissingET* met = (MissingET*) branchMET->At(0);
    if ( met->MET<40 ) continue;
    //std::cout<<"Pass MET cut"<<std::endl;
    nevt->Fill(6);
    ofc.GetTH1("nJet")->Fill( selJet.size() );

    if ( nBJet(branchParticle, selJet) <1) continue;
    nevt->Fill(7);

  
    vector<Track*> tracks;
    SecVtx *jpsi = nullptr, *d0 = nullptr, *dstar = nullptr;
    bool flag_jpsi = false, flag_d0 = false, flag_dstar = false;


    for( unsigned int jet_idx = 0 ; jet_idx < selJet.size() ; ++jet_idx) {
      Jet* jet = selJet[jet_idx].second;
      tracks.clear();
      for( unsigned int dau_idx =0 ; dau_idx < jet->Constituents.GetEntriesFast(); ++dau_idx ) {
        // temporary!!!!
        TObject* object = jet->Constituents.At(dau_idx);
        // Check if the constituent is accessible
        if(object == 0) continue;
        if (object->IsA() == Track::Class())
        {
          Track* track = (Track*) object;
          if ( track->PT < 1 ) continue;
          tracks.push_back(track);
        }
      }
      std::sort( tracks.begin(), tracks.end(), []( Track* a, Track* b) { return a->PT > b->PT; });
      std::cout<<"track idx : "<<tracks.size()<<" ";
      for( auto track : tracks) {
        std::cout<< track->PT << "( "<<track->PID<<")  ";
      }
      std::cout<<std::endl;
      /*
      int nLep = 0;
      for( int i= 0 ; i<tracks.size() ; i++) {
        if ( abs(tracks[i]->PID)==11 || abs(tracks[i]->PID)==13) nLep++;
      }
      */

      if ( tracks.size() <2 ) continue;
      for( unsigned int firstTrack = 0 ; firstTrack< tracks.size()-1 ; firstTrack++) {
        for( unsigned int secondTrack = firstTrack+1 ; secondTrack< tracks.size() ; secondTrack++) {
          int pid1 = tracks[firstTrack]->PID;
          int pid2 = tracks[secondTrack]->PID;
          Int_t pidMul = tracks[firstTrack]->PID*tracks[secondTrack]->PID;
          if ( pidMul >0 )  continue;
          // For J/Psi, 
          if ( pidMul == -121 || pidMul == -169 ) {
            TLorentzVector firstTrackLV = TracktoTLV(tracks[firstTrack]);
            TLorentzVector secondTrackLV = TracktoTLV(tracks[secondTrack]);

            TLorentzVector jpsiCand = firstTrackLV+secondTrackLV;
            if ( jpsiCand.M()>2 && jpsiCand.M()<4) {
              jpsi = new SecVtx( tracks[firstTrack], tracks[secondTrack], branchParticle);
              jpsi->setDRPT( dc, branchParticle, 443);
              ofc.data.init(&selLepton[0], &selLepton[1], jpsi);
              ofc.GetTree("tree")->Fill(); 
              flag_jpsi = true;
              nJpsi++;
              ofc.GetTH1("jpsi_mass")->Fill( jpsiCand.M());
              std::cout<<"jpsi"<<std::endl;
            }
          }
        }
      }
      for( unsigned int firstTrack = 0 ; firstTrack< tracks.size()-1 ; firstTrack++) {
        for( unsigned int secondTrack = firstTrack+1 ; secondTrack< tracks.size() ; secondTrack++) {
          int pid1 = tracks[firstTrack]->PID;
          int pid2 = tracks[secondTrack]->PID;
          Int_t pidMul = tracks[firstTrack]->PID*tracks[secondTrack]->PID;
          if ( pidMul >0 )  continue;

          // For D0,
          if ( pidMul == -67731 ) {
            bool pionplus = false;
            // If D0-> pion- kaon+,
            if ( pid1 == 211 || pid2 == 211   )  { pionplus = true;  }
            // If D0bar-> pion+ kaon-,
            if ( pid1 == -211 || pid2 == -211 )  { pionplus = false;  }
            // Third Track for soft lepton( e-, mu- ) ,
            int nsoftlep = 0;
            for( unsigned int softleptonTrack = 0 ; softleptonTrack < tracks.size() ; softleptonTrack++) {
              if ( firstTrack == softleptonTrack ) continue;
              if ( secondTrack == softleptonTrack ) continue;

              if ( pionplus && (tracks[softleptonTrack]->PID == 11 || tracks[softleptonTrack]->PID ==13) )    nsoftlep++;
              if ( !pionplus && (tracks[softleptonTrack]->PID == -11 || tracks[softleptonTrack]->PID ==-13) ) nsoftlep++;
            }
            if ( nsoftlep ==0 ) continue; 
                 
            auto d0Cand = TracktoTLV(tracks[firstTrack])+TracktoTLV(tracks[secondTrack]);
            if ( d0Cand.M() < 1 || d0Cand.M()>3 ) continue;
            d0 = new SecVtx( tracks[firstTrack], tracks[secondTrack], branchParticle);
            if ( pionplus) d0->setDRPT( dc, branchParticle, 421); 
            else           d0->setDRPT( dc, branchParticle, -421); 
            d0->setSoftLepton(nsoftlep); 
            ofc.data.init(&selLepton[0], &selLepton[1], d0);
            ofc.GetTree("tree")->Fill(); 
            flag_d0 = true;
            nD0++;
            ofc.GetTH1("d0_mass")->Fill( d0Cand.M());
            std::cout<<"D0"<<std::endl;

            if ( abs( d0Cand.M() - 1.864 ) > d0_massWindow*2 ) continue;
            std::cout<<"dstar cand"<<std::endl; 
            for( unsigned int pionTrack = 0 ; pionTrack < tracks.size() ; pionTrack++) {
              if ( firstTrack == pionTrack ) continue;
              if ( secondTrack == pionTrack ) continue;
              if ( tracks[pionTrack]->PID * tracks[firstTrack]->PID == 44521 || tracks[pionTrack]->PID*tracks[secondTrack]->PID == 44521 ) {
                
                auto dstarCand = d0Cand+ TracktoTLV(tracks[pionTrack]);
                float diffMass = dstarCand.M() - d0Cand.M(); 
                if ( diffMass > 0.135 && diffMass<0.170) {
                  dstar = new SecVtx( d0, tracks[pionTrack], branchParticle);
                  dstar->setDRPT( dc, branchParticle, 413*dstar->charge() ); 
                  dstar->setSoftLepton(nsoftlep);
                  dstar->setD0mass( d0Cand.M()); 
                  ofc.GetTH1("dstar_mass")->Fill( dstarCand.M());
                  ofc.GetTH1("dstar_diffmass")->Fill( dstarCand.M()- d0Cand.M());
                  flag_dstar = true;
                  std::cout<<"D*"<<std::endl;
                  nDstar++;
                  ofc.data.init(&selLepton[0], &selLepton[1], dstar);
                  ofc.GetTree("tree")->Fill(); 
                }
              } 
            }
          }
        }
      }
    }
    int cat = 0; 
    if ( flag_jpsi ) {
      cat = 2;
    }
    else if ( flag_dstar) {
      cat = 3; 
    }
    else if ( flag_d0) {
      cat = 4; 
    }
    else cat = 1; 
    ofc.GetTH1("sv_category")->Fill(cat);
  }        
  std::cout<<TString::Format("Filling! nJpsi: %d , nD0 : %d, nDstar : %d\n",nJpsi,nD0, nDstar)<<std::endl;
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
  TTree* tree = new TTree("tree", "tree");
  tree->Branch(treeName.c_str(), &ofc.data, treeTypes );
  ofc.AddTree(tree);
}

void addBranch(OutFileClass& ofc, std::string treeName, void* address, const char* treeTypes) {
  ofc.GetFile()->cd();
  TTree* tree = (TTree*)ofc.GetFile()->Get("tree");
  tree->Branch(treeName.c_str(), address, treeTypes);
}
void BookingHist(OutFileClass& ofc)
{
  ofc.GetFile()->cd();
  TH1F* h1 = new TH1F("nEvent","Number of Events",7,1,8);
  auto xaxis = h1->GetXaxis();
  xaxis->SetBinLabel(1,"Total");
  xaxis->SetBinLabel(2,"GenEvt");
  xaxis->SetBinLabel(3,"nlep2");
  xaxis->SetBinLabel(4,"ZVeto");
  xaxis->SetBinLabel(5,"nJet2");
  xaxis->SetBinLabel(6,"MET40");
  xaxis->SetBinLabel(7,"nbJet1");
  ofc.AddTH1(h1);

  TH1F* h2 = new TH1F("nJet","Number of Jets",10,0,10);
  ofc.AddTH1(h2);

  TH1F* h3 = new TH1F("jpsi_mass","Invariant mass of J/#Psi", 100,2,4);
  ofc.AddTH1(h3);
  TH1F* h4 = new TH1F("d0_mass","Invariant mass of D0 Candidate", 100,1.6,2.2);
  ofc.AddTH1(h4);
  TH1F* h5 = new TH1F("dstar_diffmass","Invariant diffmass of D*-D0", 100,0.135,0.17);
  ofc.AddTH1(h5);
  TH1F* h6 = new TH1F("dstar_mass","Invariant mass of D*", 100,1.6,2.2);
  ofc.AddTH1(h6);
  TH1F* h7 = new TH1F("sv_category","Categorization of charmed meson", 4,1,5);
  auto xaxis2 = h7->GetXaxis();
  xaxis2->SetBinLabel(1,"Other");
  xaxis2->SetBinLabel(2,"J/#psi");
  xaxis2->SetBinLabel(3,"D^{*}");
  xaxis2->SetBinLabel(4,"D^{0}");
  ofc.AddTH1(h7);


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
  TString brC = TString("");
  // Common
  //brC += "dataType/I:";
  // Lepton
  brC += "lep_charge[2]/I:lep_pid[2]/I:lep_pt[2]/F:lep_eta[2]/F:lep_phi[2]/F:lep_mass[2]/F:";
  // SV
  brC += "sv_pid/I:sv_charge/I:sv_isFromTop/I:sv_pt/F:sv_eta/F:sv_phi/F:sv_mass/F:sv_diffmass/F:sv_d0mass/F:sv_vx/F:sv_vy/F:sv_vz/F:sv_dx/F:sv_dy/F:sv_dz/F:sv_L3D/F:sv_LXY/F:sv_softlep/I:sv_dRTrue/F:sv_delPtTrue/F:sv_dau_pid[3]/F:sv_dau_pt[3]/F:sv_dau_eta[3]/F:sv_dau_phi[3]/F:sv_dau_mass[3]/F";
  // lSV
  TString brC2 = TString("");
  brC2 += "lsv_pid[2]/I:lsv_pt[2]/F:lsv_eta[2]/F:lsv_phi[2]/F:lsv_mass[2]/F:lsv_diffmass[2]/F:lsv_d0mass[2]/F";

  std::cout<<brC<<std::endl;
  BookingTree(ofc, "SVTree", brC.Data());
  addBranch( ofc,"LSVTree", &ofc.data.lsv_pid[0], brC2.Data());                      
  AnalyseEvents(treeReader,ofc);

  Write(ofc);
  cout << "** Exiting..." << endl;

  file->Close();
  //delete result;
  delete treeReader;
  delete chain;

}
