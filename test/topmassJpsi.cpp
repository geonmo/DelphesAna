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
  SecVtx(Track* trackA, Track* trackB) {

    float track1_mass = 0, track2_mass = 0;
    
    auto trackA_TLV = TracktoTLV(trackA);
    auto trackB_TLV = TracktoTLV(trackB);
    pos_ = trackA_TLV + trackB_TLV;

    dau_pid_[0] = trackA->PID;
    dau_pid_[1] = trackB->PID;
    charge_ = (trackA->Charge + trackB->Charge);
    daus_[0] = trackA_TLV;
    daus_[1] = trackB_TLV;
    vx_[0] = trackA->X;
    vx_[1] = trackB->X;
    vy_[0] = trackA->Y;
    vy_[1] = trackB->Y;
    vz_[0] = trackA->Z;
    vz_[1] = trackB->Z;
    nDau_ = 2;
  }
  SecVtx(SecVtx* secvtxA, Track* trackB) {
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
    vx_[0] = secvtxA->vx(0);
    vx_[1] = secvtxA->vx(1);
    vx_[2] = trackB->X;
    vy_[0] = secvtxA->vy(0);
    vy_[1] = secvtxA->vy(1);
    vy_[2] = trackB->Y;
    vz_[0] = secvtxA->vz(0);
    vz_[1] = secvtxA->vz(1);
    vz_[2] = trackB->Z;
    nDau_ = 3;
  }
  TLV P4() { return pos_; } 
  Float_t mass() { return pos_.M(); }
  Float_t pt() { return pos_.Pt(); }
  Float_t eta() { return pos_.Eta(); }
  Float_t phi() { return pos_.Phi(); }
  Int_t charge() { return charge_; }
  Int_t dau_pid(Int_t idx) { return dau_pid_[idx]; }
  TLV dau(Int_t idx) { 
    if ( idx>3) { std::cout<<"Wrong idx. It is too high."<<std::endl; return TLorentzVector(); }
    //if ( idx > nDau_-1 ) { std::cout<<"Wrong idx. It is null ptr."<<std::endl; return TLorentzVector(); } 
    return daus_[idx]; 
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
  Float_t vxd() {
    if ( nDau_ ==2 ) {
      float vxd_temp = vx_[0]-vx_[1];
      return TMath::Sqrt( vxd_temp*vxd_temp) ;
    }
    else if ( nDau_ ==3 ) {
      float vxd_temp1 = vx_[0]-vx_[1];
      float vxd_temp2 = vx_[0]-vx_[2];
      float vxd_temp3 = vx_[1]-vx_[2];
      float vxd__ = (TMath::Sqrt(vxd_temp1*vxd_temp1)+TMath::Sqrt(vxd_temp2*vxd_temp2)+TMath::Sqrt(vxd_temp3*vxd_temp3)) /3.;
      return vxd__;
    }
    else return 0; 
  }
  Float_t vyd() {
    if ( nDau_ ==2 ) {
      float vyd_temp = vy_[0]-vy_[1];
      return TMath::Sqrt( vyd_temp*vyd_temp) ;
    }
    else if ( nDau_ ==3 ) {
      float vyd_temp1 = vy_[0]-vy_[1];
      float vyd_temp2 = vy_[0]-vy_[2];
      float vyd_temp3 = vy_[1]-vy_[2];
      float vyd__ = (TMath::Sqrt(vyd_temp1*vyd_temp1)+TMath::Sqrt(vyd_temp2*vyd_temp2)+TMath::Sqrt(vyd_temp3*vyd_temp3)) /3.;
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
      float vzd__ = (TMath::Sqrt(vzd_temp1*vzd_temp1)+TMath::Sqrt(vzd_temp2*vzd_temp2)+TMath::Sqrt(vzd_temp3*vzd_temp3)) /3.;
      return vzd__;
    }
    else return 0; 
  }
  void setDRPT(DecayChannel dc, TClonesArray* branchParticle, int pdgId) {
    int mcTruth = dc.SearchParticle(branchParticle, pdgId, pos_ );
    if ( mcTruth == -1 ) return;
    auto gen = (GenParticle*)branchParticle->At( mcTruth );
    dRTrue_ = gen->P4().DeltaR( pos_ );
    delPtTrue_ = abs(gen->PT-pos_.Pt())/gen->PT;

    int isFromTop =0;
    int isFromTopIdx = dc.isFromTop( branchParticle, mcTruth );
    if ( isFromTopIdx != -1 ) isFromTop_ = ((GenParticle*)branchParticle->At( isFromTopIdx))->PID;
  }



//float dRTrue, float delPtTrue) { dRTrue_ = dRTrue; delPtTrue_ = delPtTrue_; }
  float DR() { return dRTrue_; }
  float DPT() { return delPtTrue_; }

  Int_t isFromTop() { return isFromTop_;}
  Int_t nDau() { return nDau_; }
private :
  Int_t charge_;
  TLV pos_;
  TLV daus_[3];
  Int_t nDau_;
  Int_t dau_pid_[3];
  Float_t vx_[3];
  Float_t vy_[3];
  Float_t vz_[3];
  Int_t isFromTop_;
  Float_t delPtTrue_;
  Float_t dRTrue_;
};    

std::pair< TLorentzVector, TLorentzVector> lsvPairing( QLV lep1 , QLV lep2, SecVtx* sv) {
  int lep1_pid = lep1.first;
  int lep2_pid = lep2.first;
  TLorentzVector lep1_TLV = common::LVtoTLV(lep1.second);
  TLorentzVector lep2_TLV = common::LVtoTLV(lep2.second);

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

struct LSV {
  float pt[2];
  float eta[2];
  float phi[2];
  float mass[2];
  float vx[2];
  float vy[2];
  float vz[2];
  float vxd[2];
  float vyd[2];
  float vzd[2];
};


class saveData {
public :
  Int_t dataType;
  Int_t lep_charge[2], lep_pid[2];
  Float_t lep_pt[2], lep_eta[2], lep_phi[2], lep_mass[2];

  Int_t jpsi_charge, jpsi_isFromTop;
  Float_t jpsi_pt, jpsi_eta, jpsi_phi, jpsi_mass;
  Float_t jpsi_dx, jpsi_dy, jpsi_dz;
  Float_t jpsi_dRTrue, jpsi_delPtTrue;
  

  Int_t jpsi_dau_pid[2];
  Float_t jpsi_dau_pt[2], jpsi_dau_eta[2], jpsi_dau_phi[2], jpsi_dau_mass[2];

  Int_t d0_charge, d0_isFromTop;
  Float_t d0_pt, d0_eta, d0_phi, d0_mass;
  Float_t d0_dx, d0_dy, d0_dz;
  Float_t d0_dRTrue, d0_delPtTrue;

  Int_t d0_dau_pid[2];
  Float_t d0_dau_pt[2], d0_dau_eta[2], d0_dau_phi[2], d0_dau_mass[2];

  Int_t dstar_charge, dstar_isFromTop;
  Float_t dstar_pt, dstar_eta, dstar_phi, dstar_mass, dstar_diffmass;
  Float_t dstar_dx, dstar_dy, dstar_dz;
  Float_t dstar_dRTrue, dstar_delPtTrue;

  Int_t dstar_dau_pid[3];
  Float_t dstar_dau_pt[3], dstar_dau_eta[3], dstar_dau_phi[3], dstar_dau_mass[3];


  Float_t ljpsi_pt[2], ljpsi_eta[2], ljpsi_phi[2], ljpsi_mass[2];
  Float_t ljpsi_vx[2], ljpsi_vy[2],  ljpsi_vz[2];
  Float_t ljpsi_dx[2], ljpsi_dy[2],  ljpsi_dz[2];
  Float_t ld0_pt[2], ld0_eta[2], ld0_phi[2], ld0_mass[2];
  Float_t ld0_vx[2], ld0_vy[2],  ld0_vz[2];
  Float_t ld0_dx[2], ld0_dy[2],  ld0_dz[2];
  Float_t ldstar_pt[2], ldstar_eta[2], ldstar_phi[2], ldstar_mass[2];
  Float_t ldstar_vx[2], ldstar_vy[2],  ldstar_vz[2];
  Float_t ldstar_dx[2], ldstar_dy[2],  ldstar_dz[2];
 


  void reset() {
    dataType=0;
    jpsi_charge =0; 
    jpsi_isFromTop=-9; 
    jpsi_pt = -9.0; 
    jpsi_eta = -9.0; 
    jpsi_phi = -9.0; 
    jpsi_mass = -9.0;
    jpsi_dx = 0.f;
    jpsi_dy = 0.f;
    jpsi_dz = 0.f;
    jpsi_dRTrue = 999.f;
    jpsi_delPtTrue = 1.f;
 
    d0_charge =0;  
    d0_isFromTop=-9; 
    d0_pt = -9.0; 
    d0_eta = -9.0; 
    d0_phi = -9.0; 
    d0_mass = -9.0; 
    d0_dx = 0.f;
    d0_dy = 0.f;
    d0_dz = 0.f;
    d0_dRTrue = 999.f;
    d0_delPtTrue = 1.f;

    dstar_charge =0;  
    dstar_isFromTop=-9; 
    dstar_pt = -9.0; 
    dstar_eta = -9.0; 
    dstar_phi = -9.0; 
    dstar_mass = -9.0; 
    dstar_diffmass = -9.0; 
    dstar_dx = 0.f;
    dstar_dy = 0.f;
    dstar_dz = 0.f;
    dstar_dRTrue = 999.f;
    dstar_delPtTrue = 1.f;
    for( int idx = 0; idx< 3 ; idx++) { 
      dstar_dau_pid[idx] =0;  
      dstar_dau_pt[idx] = -9.0; 
      dstar_dau_eta[idx] = -9.0; 
      dstar_dau_phi[idx] = -9.0; 
      dstar_dau_mass[idx] = -9.0; 
    }
    for( int idx = 0; idx< 2 ; idx++) { 
      jpsi_dau_pid[idx] =0;  
      jpsi_dau_pt[idx] = -9.0; 
      jpsi_dau_eta[idx] = -9.0; 
      jpsi_dau_phi[idx] = -9.0; 
      jpsi_dau_mass[idx] = -9.0;
 
      d0_dau_pid[idx] =0;  
      d0_dau_pt[idx] = -9.0; 
      d0_dau_eta[idx] = -9.0; 
      d0_dau_phi[idx] = -9.0; 
      d0_dau_mass[idx] = -9.0;
 
      lep_charge[idx] =0;  
      lep_pid[idx] =0 ;
      lep_pt[idx] = -9.0; 
      lep_eta[idx] = -9.0; 
      lep_phi[idx] = -9.0; 
      lep_mass[idx] = -9.0; 

      ljpsi_pt[idx] = -9.0;
      ljpsi_eta[idx] = -9.0;
      ljpsi_phi[idx] = -9.0;
      ljpsi_mass[idx] = -9.0;
      ljpsi_vx[idx] = 0.f;
      ljpsi_vy[idx] = 0.f;
      ljpsi_vz[idx] = 0.f;
      ljpsi_dx[idx] = 0.f;
      ljpsi_dy[idx] = 0.f;
      ljpsi_dz[idx] = 0.f;
 
      ld0_pt[idx] = -9.0;
      ld0_eta[idx] = -9.0;
      ld0_phi[idx] = -9.0;
      ld0_mass[idx] = -9.0;
      ld0_vx[idx] = 0.f;
      ld0_vy[idx] = 0.f;
      ld0_vz[idx] = 0.f;
      ld0_dx[idx] = 0.f;
      ld0_dy[idx] = 0.f;
      ld0_dz[idx] = 0.f;
 
      ldstar_pt[idx] = -9.0;
      ldstar_eta[idx] = -9.0;
      ldstar_phi[idx] = -9.0;
      ldstar_mass[idx] = -9.0;

      ldstar_vx[idx] = 0.f;
      ldstar_vy[idx] = 0.f;
      ldstar_vz[idx] = 0.f;
      ldstar_dx[idx] = 0.f;
      ldstar_dy[idx] = 0.f;
      ldstar_dz[idx] = 0.f;
 
    }
  }
  void init(Int_t lep1_pid, Int_t lep2_pid, LV* lep1, LV* lep2, SecVtx* jpsi, SecVtx* d0, SecVtx* dstar) {
    reset();
    if ( jpsi != nullptr) {
      jpsi_charge = jpsi->charge();
      jpsi_isFromTop = jpsi->isFromTop();
      jpsi_pt = jpsi->pt();
      jpsi_eta = jpsi->eta();
      jpsi_phi = jpsi->phi();
      jpsi_mass = jpsi->mass();
      jpsi_dx = jpsi->vxd();
      jpsi_dy = jpsi->vyd();
      jpsi_dz = jpsi->vzd();
      jpsi_dRTrue = jpsi->DR();
      jpsi_delPtTrue = jpsi->DPT();

      for ( int idx =0 ; idx<2 ; idx++) {
        jpsi_dau_pid[idx]   = jpsi->dau_pid(idx);
        jpsi_dau_pt[idx]   = jpsi->dau(idx).Pt();
        jpsi_dau_eta[idx]  = jpsi->dau(idx).Eta();
        jpsi_dau_phi[idx]  = jpsi->dau(idx).Phi();
        jpsi_dau_mass[idx] = jpsi->dau(idx).M();

        ljpsi_vx[idx] = jpsi->vx(0);
        ljpsi_vy[idx] = jpsi->vy(0);
        ljpsi_vz[idx] = jpsi->vz(0);
        ljpsi_dx[idx] = jpsi->vxd();
        ljpsi_dy[idx] = jpsi->vyd();
        ljpsi_dz[idx] = jpsi->vzd();
      }
    }

    if ( d0 != nullptr) {
      d0_charge = d0->charge();
      d0_isFromTop = d0->isFromTop();
      d0_pt = d0->pt();
      d0_eta = d0->eta();
      d0_phi = d0->phi();
      d0_mass = d0->mass();
      d0_dx = d0->vxd();
      d0_dy = d0->vyd();
      d0_dz = d0->vzd();
      d0_dRTrue = d0->DR();
      d0_delPtTrue = d0->DPT();
      for ( int idx =0 ; idx<2 ; idx++) {
        d0_dau_pid[idx]  = d0->dau_pid(idx);
        d0_dau_pt[idx]   = d0->dau(idx).Pt();
        d0_dau_eta[idx]  = d0->dau(idx).Eta();
        d0_dau_phi[idx]  = d0->dau(idx).Phi();
        d0_dau_mass[idx] = d0->dau(idx).M();

        ld0_vx[idx] = d0->vx(0);
        ld0_vy[idx] = d0->vy(0);
        ld0_vz[idx] = d0->vz(0);
        ld0_dx[idx] = d0->vxd();
        ld0_dy[idx] = d0->vyd();
        ld0_dz[idx] = d0->vzd();
      }
    }

    if ( dstar != nullptr) { 
      dstar_charge = dstar->charge();
      dstar_isFromTop = dstar->isFromTop();
      dstar_pt = dstar->pt();
      dstar_eta = dstar->eta();
      dstar_phi = dstar->phi();
      dstar_mass = dstar->mass();
      dstar_diffmass = dstar->mass() - d0->mass();
      dstar_dx = dstar->vxd();
      dstar_dy = dstar->vyd();
      dstar_dz = dstar->vzd();
      dstar_dRTrue = dstar->DR();
      dstar_delPtTrue = dstar->DPT();
      for ( int idx =0 ; idx<3 ; idx++) {
        dstar_dau_pid[idx]  = dstar->dau_pid(idx);
        dstar_dau_pt[idx]   = dstar->dau(idx).Pt();
        dstar_dau_eta[idx]  = dstar->dau(idx).Eta();
        dstar_dau_phi[idx]  = dstar->dau(idx).Phi();
        dstar_dau_mass[idx] = dstar->dau(idx).M();
      }
      for ( int idx =0 ; idx<2 ; idx++) {
        ldstar_vx[idx] = dstar->vx(0);
        ldstar_vy[idx] = dstar->vy(0);
        ldstar_vz[idx] = dstar->vz(0);
        ldstar_dx[idx] = dstar->vxd();
        ldstar_dy[idx] = dstar->vyd();
        ldstar_dz[idx] = dstar->vzd();
      }
    }

    lep_charge[0] = -lep1_pid/abs(lep1_pid);
    lep_pid[0] = lep1_pid;
    lep_pt[0] = lep1->pt();
    lep_eta[0] = lep1->eta();
    lep_phi[0] = lep1->phi();
    lep_mass[0] = lep1->mass();

    lep_charge[1] = -lep2_pid/abs(lep2_pid);
    lep_pid[1] = lep2_pid;
    lep_pt[1] = lep2->pt();
    lep_eta[1] = lep2->eta();
    lep_phi[1] = lep2->phi();
    lep_mass[1] = lep2->mass();
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
    if ( muon->PT<20 || abs(muon->Eta)>2.5) continue;
    selLeptons.push_back( make_pair(muon->Charge*-13, common::TLVtoLV(muon->P4()) ));
  }
  for( int i= 0 ; i< electrons->GetEntriesFast() ; ++i) {
    Electron* electron = (Electron*) electrons->At(i);
    if ( electron->PT<20 || abs(electron->Eta)>2.5) continue;
    selLeptons.push_back( make_pair(electron->Charge*-11, common::TLVtoLV(electron->P4()) ));
  }
  return selLeptons; 
}
std::vector<QJET> selectedJet( TClonesArray* jets ) {
  std::vector<QJET> selJets;
  for(int i = 0; i < jets->GetEntriesFast(); ++i)
  {
    Jet* jet = (Jet*) jets->At(i);
    if ( jet->PT < 30 || abs(jet->Eta)>2.5 ) continue;
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

    // Least POG Numbers, // Simple event Selection.
    //D( std::cout<<"This is dilepton events!"<<std::endl; )
    std::vector<QLV> selLepton = selectedLepton( branchMuon, branchElectron);
    if ( selLepton.size()<2  || selLepton[0].first*selLepton[1].first > 0  ) continue;
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
    SecVtx *jpsi = nullptr, *d0 = nullptr,*dstar = nullptr;
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
      if ( tracks.size() <=1 ) continue;
      for( unsigned int firstTrack = 0 ; firstTrack< tracks.size()-1 ; firstTrack++) {
        //int isFromTop=0, isFromTopIdx=-1;
        //int motherIdx = ((GenParticle*)(tracks[firstTrack]->Particle.GetObject()))->M1;
        //if ( motherIdx ==-1 ) isFromTop = 0; 
        //else { 
        //  isFromTopIdx = dc.isFromTop( branchParticle, motherIdx );
        //  if ( isFromTopIdx != -1 ) isFromTop = ((GenParticle*)branchParticle->At( isFromTopIdx))->PID;
        //  else isFromTop=0;
        //}
        for( unsigned int secondTrack = firstTrack+1 ; secondTrack< tracks.size() ; secondTrack++) {
          Int_t pidMul = tracks[firstTrack]->PID*tracks[secondTrack]->PID;
          if ( pidMul >0 )  continue;
          // For J/Psi, 
          if ( !flag_jpsi && (pidMul == -121 || pidMul == -169 )) {
            TLorentzVector firstTrackLV = TracktoTLV(tracks[firstTrack]);
            TLorentzVector secondTrackLV = TracktoTLV(tracks[secondTrack]);

            TLorentzVector jpsiCand = firstTrackLV+secondTrackLV;
            if ( jpsiCand.M() > 2 && jpsiCand.M()<4 ) { 
              ofc.GetTH1("jpsi_mass")->Fill( jpsiCand.M());
              jpsi = new SecVtx( tracks[firstTrack], tracks[secondTrack]);
              jpsi->setDRPT( dc, branchParticle, 443); 
              //jpsi->setFromTop(isFromTop);
              flag_jpsi = true;
              nJpsi++;
              std::cout<<"jpsi"<<std::endl;
            }
          }
          // For D0,
          if ( pidMul == -67731) {
            TLorentzVector d0Cand = TracktoTLV(tracks[firstTrack])+TracktoTLV(tracks[secondTrack]);
            if ( d0Cand.M()>1.6 && d0Cand.M()<2.2) {
              for( unsigned int pionTrack = 0 ; pionTrack < tracks.size() ; pionTrack++) {
                if ( firstTrack == pionTrack ) continue;
                if ( secondTrack == pionTrack ) continue;
                if ( !flag_dstar && (tracks[pionTrack]->PID * tracks[firstTrack]->PID == 44521 || tracks[pionTrack]->PID*tracks[secondTrack]->PID == 44521 )) {
                  auto dstarCand = d0Cand+ TracktoTLV(tracks[pionTrack]);
                  float diffMass = dstarCand.M() - d0Cand.M(); 
                  if ( diffMass> 0.135 && diffMass<0.170 ) {
                    if ( !flag_d0 ) nD0++;
                    d0 = new SecVtx( tracks[firstTrack], tracks[secondTrack]);
                    d0->setDRPT( dc, branchParticle, 421); 
                    flag_d0 = true;

                    dstar = new SecVtx( d0, tracks[pionTrack]);
                    dstar->setDRPT( dc, branchParticle, 421*dstar->charge() ); 
                    ofc.GetTH1("dstar_mass")->Fill( dstarCand.M());
                    ofc.GetTH1("dstar_diffmass")->Fill( dstarCand.M()- d0Cand.M());
                    flag_dstar = true;
                    std::cout<<"D*"<<std::endl;
                    nDstar++;
                    // if D* is found, D0 also be kept.
                    ofc.GetTH1("d0_mass")->Fill( d0Cand.M());
                    std::cout<<"D0"<<std::endl;
                  }
                }
              }
              // If D0 is found first, it will be kept.
              if ( (!flag_d0)  ) {
                ofc.GetTH1("d0_mass")->Fill( d0Cand.M());
                std::cout<<"D0"<<std::endl;
                nD0++;
                d0 = new SecVtx( tracks[firstTrack], tracks[secondTrack]);
                d0->setDRPT( dc, branchParticle, 421); 
                //d0->setFromTop(isFromTop);
                flag_d0 = true;
              }
            }
          }
        }
      }
    }
    ofc.data.init(selLepton[0].first, selLepton[1].first, &(selLepton[0].second), &(selLepton[1].second), jpsi, d0, dstar);
    if ( jpsi != nullptr) {

      std::pair<TLorentzVector,TLorentzVector> pair = lsvPairing(selLepton[0], selLepton[1], jpsi);
      TLorentzVector ljpsi1 = pair.first;
      TLorentzVector ljpsi2 = pair.second;

      ofc.data.ljpsi_pt[0] = ljpsi1.Pt();
      ofc.data.ljpsi_eta[0] = ljpsi1.Eta();
      ofc.data.ljpsi_phi[0] = ljpsi1.Phi();
      ofc.data.ljpsi_mass[0] = ljpsi1.M();

      ofc.data.ljpsi_pt[1] = ljpsi2.Pt();
      ofc.data.ljpsi_eta[1] = ljpsi2.Eta();
      ofc.data.ljpsi_phi[1] = ljpsi2.Phi();
      ofc.data.ljpsi_mass[1] = ljpsi2.M();
    }
    if ( d0 != nullptr) {
      std::pair<TLorentzVector,TLorentzVector> pair = lsvPairing(selLepton[0], selLepton[1], d0);
      TLorentzVector ld01 = pair.first;
      TLorentzVector ld02 = pair.second;

      ofc.data.ld0_pt[0] = ld01.Pt();
      ofc.data.ld0_eta[0] = ld01.Eta();
      ofc.data.ld0_phi[0] = ld01.Phi();
      ofc.data.ld0_mass[0] = ld01.M();

      ofc.data.ld0_pt[1] = ld02.Pt();
      ofc.data.ld0_eta[1] = ld02.Eta();
      ofc.data.ld0_phi[1] = ld02.Phi();
      ofc.data.ld0_mass[1] = ld02.M();
    }
    if ( dstar != nullptr) {
      std::pair<TLorentzVector,TLorentzVector> pair = lsvPairing(selLepton[0], selLepton[1], dstar);
      TLorentzVector ldstar1 = pair.first;
      TLorentzVector ldstar2 = pair.second;

      ofc.data.ldstar_pt[0] = ldstar1.Pt();
      ofc.data.ldstar_eta[0] = ldstar1.Eta();
      ofc.data.ldstar_phi[0] = ldstar1.Phi();
      ofc.data.ldstar_mass[0] = ldstar1.M();

      ofc.data.ldstar_pt[1] = ldstar2.Pt();
      ofc.data.ldstar_eta[1] = ldstar2.Eta();
      ofc.data.ldstar_phi[1] = ldstar2.Phi();
      ofc.data.ldstar_mass[1] = ldstar2.M();

    }
    ofc.data.dataType= dc.channel();
    ofc.GetTree("tree")->Fill();
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
  TH1F* h6 = new TH1F("vertexDiff","Resolution of vertex", 1000,-1e-5,1e-5);
  ofc.AddTH1(h6);
  TH2F* h7 = new TH2F("vertexDiff2D","Resolution of vertex", 1000,-1e-5,1e-5,1000,-1e-5,1e-5);
  ofc.AddTH1(h7);
  TH1F* h8 = new TH1F("dstar_mass","Invariant mass of D*", 100,1.6,2.2);
  ofc.AddTH1(h8);


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
  brC += "dataType/I:";
  // Lepton
  brC += "lep_charge[2]/I:lep_pid[2]/I:lep_pt[2]/F:lep_eta[2]/F:lep_phi[2]/F:lep_mass[2]/F:jpsi_charge/I:jpsi_isFromTop/I:";
  // Jpsi
  brC += "jpsi_pt/F:jpsi_eta/F:jpsi_phi/F:jpsi_mass/F:jpsi_dx/F:jpsi_dy/F:jpsi_dz/F:jpsi_dRTrue/F:jpsi_delPtTrue/F:jpsi_dau_pid[2]/I:jpsi_dau_pt[2]/F:jpsi_dau_eta[2]/F:jpsi_dau_phi[2]/F:jpsi_dau_mass[2]/F:";
  // d0
  brC += "d0_charge/I:d0_isFromTop/I:d0_pt/F:d0_eta/F:d0_phi/F:d0_mass/F:d0_dx/F:d0_dy/F:d0_dz/F:d0_dRTrue/F:d0_delPtTrue/F:d0_dau_pid[2]/I:d0_dau_pt[2]/F:d0_dau_eta[2]/F:d0_dau_phi[2]/F:d0_dau_mass[2]/F:";
  // dstar
  brC += "dstar_charge/I:dstar_isFromTop/I:dstar_pt/F:dstar_eta/F:dstar_phi/F:dstar_mass/F:dstar_diffmass/F:dstar_dx/F:dstar_dy/F:dstar_dz/F:dstar_dRTrue/F:dstar_delPtTrue/F:dstar_dau_pid[3]/F:dstar_dau_pt[3]/F:dstar_dau_eta[3]/F:dstar_dau_phi[3]/F:dstar_dau_mass[3]/F";
  // lSV
  TString brC2 = TString("");
  brC2 += "ljpsi_pt[2]/F:ljpsi_eta[2]/F:ljpsi_phi[2]/F:ljpsi_mass[2]/F:ljpsi_vx[2]/F:ljpsi_vy[2]/F:ljpsi_vz[2]/F:ljpsi_dx[2]/F:ljpsi_dy[2]/F:ljpsi_dz[2]/F:";
  brC2 += "ld0_pt[2]/F:ld0_eta[2]/F:ld0_phi[2]/F:ld0_mass[2]/F:ld0_vx[2]/F:ld0_vy[2]/F:ld0_vz[2]/F:ld0_dx[2]/F:ld0_dy[2]/F:ld0_dz[2]/F:";
  brC2 += "ldstar_pt[2]/F:ldstar_eta[2]/F:ldstar_phi[2]/F:ldstar_mass[2]/F:ldstar_vx[2]/F:ldstar_vy[2]/F:ldstar_vz[2]/F:ldstar_dx[2]/F:ldstar_dy[2]/F:ldstar_dz[2]/F";

  std::cout<<brC<<std::endl;
  BookingTree(ofc, "SVTree", brC.Data());
  addBranch( ofc,"LSVTree", &ofc.data.ljpsi_pt[0], brC2.Data());                      
  //BookingTree(ofc,"JpsiTree", "dataType/I:lep_charge[2]/I:lep_pid[2]/I:lep_pt[2]/F:lep_eta[2]/F:lep_phi[2]/F:lep_mass[2]/F:jpsi_charge/I:jpsi_isFromTop/I:jpsi_pt/F:jpsi_eta/F:jpsi_phi/F:jpsi_mass/F:jpsi_dx/F:jpsi_dy/F:jpsi_dz/F:jpsi_dRTrue/F:jpsi_delPtTrue/F:jpsi_dau_pid[2]/I:jpsi_dau_pt[2]/F:jpsi_dau_eta[2]/F:jpsi_dau_phi[2]/F:jpsi_dau_mass[2]/F:d0_charge/I:d0_isFromTop/I:d0_pt/F:d0_eta/F:d0_phi/F:d0_mass/F:d0_dx/F:d0_dy/F:d0_dz/F:d0_dRTrue/F:d0_delPtTrue/F:d0_dau_pid[2]/I:d0_dau_pt[2]/F:d0_dau_eta[2]/F:d0_dau_phi[2]/F:d0_dau_mass[2]/F:dstar_charge/I:dstar_isFromTop/I:dstar_pt/F:dstar_eta/F:dstar_phi/F:dstar_mass/F:dstar_diffmass/F:dstar_dx/F:dstar_dy/F:dstar_dz/F:dstar_dRTrue/F:dstar_delPtTrue/F:dstar_dau_pid[3]/F:dstar_dau_pt[3]/F:dstar_dau_eta[3]/F:dstar_dau_phi[3]/F:dstar_dau_mass[3]/F");
  AnalyseEvents(treeReader,ofc);

  Write(ofc);
  cout << "** Exiting..." << endl;

  file->Close();
  //delete result;
  delete treeReader;
  delete chain;

}
