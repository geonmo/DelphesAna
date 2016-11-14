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

#include <OutFileClass.h>
#include <EventSelector.h>

using namespace std;

class ExRootTreeReader;
class saveData;

typedef std::pair<int, LV> QLV;
typedef std::pair<int, Jet*> QJET;
typedef TLorentzVector TLV;


float jpsi_massWindow = 0.1;
float d0_massWindow = 0.15;
float dstar_massWindow = 0.05;

float muonTrackPTMargin = 4;
float elecTrackPTMargin = 1;
float hadronPTMargin = 1;

float track_vertex_margin = 0.01f;



TLorentzVector TracktoTLV( Track* track, int assumePID=0 ) {
  TLV track_TLV = track->P4();
  double track_mass =0.;
  if      ( abs(track->PID)==11)  track_mass = 0.000511169;
  else if      ( abs(track->PID)==13)  track_mass = 0.105658;
  else {
    if ( abs( assumePID)==211 ) track_mass = 0.139526;
    else if ( abs( assumePID)==321) track_mass = 0.493652;
  }
  TLorentzVector newTrack;
  newTrack.SetPtEtaPhiM( track_TLV.Pt(), track_TLV.Eta(), track_TLV.Phi(), track_mass);
  return newTrack;
}


class SecVtx  { 
  public :
    SecVtx() {
      reset();
      nDau_=0;
    }
    SecVtx(Track* trackA, Track* trackB, TClonesArray* branchParticle=nullptr, float dau1_V3D=1000.0f, float dau2_V3D=1000.0f){
      reset();
      float track1_mass = 0, track2_mass = 0;

      auto trackA_TLV = TracktoTLV(trackA, 321);
      auto trackB_TLV = TracktoTLV(trackB, 211);
      pos_ = trackA_TLV + trackB_TLV;
      dau_pid_[0] = trackA->PID;
      dau_pid_[1] = trackB->PID;

      dau_V3D_[0] = dau1_V3D; 
      dau_V3D_[1] = dau2_V3D; 
      dau_V3D_[2] = 999.f;

      charge_ = (trackA->Charge + trackB->Charge);
      daus_[0] = trackA_TLV;
      daus_[1] = trackB_TLV;

      vx_[0] = trackA->X;
      vy_[0] = trackA->Y;
      vz_[0] = trackA->Z;
      vx_[1] = trackB->X;
      vy_[1] = trackB->Y;
      vz_[1] = trackB->Z;
      nDau_ = 2;
      isFromTop_= 0;
    }

    SecVtx(Track* trackA, Track* trackB, Track* trackC, TClonesArray* branchParticle=nullptr, float dau1_V3D=1000.0f, float dau2_V3D=1000.0f, float dau3_V3D=1000.0f) {
      reset();
      auto trackA_TLV = TracktoTLV(trackA, 321);
      auto trackB_TLV = TracktoTLV(trackB, 211);
      auto trackC_TLV = TracktoTLV(trackC, 211);
      pos_ = trackA_TLV+trackB_TLV+trackC_TLV;
      charge_ = trackC->Charge;
      dau_pid_[0] = trackA->PID;
      dau_pid_[1] = trackB->PID;
      dau_pid_[2] = trackC->PID;
      dau_V3D_[0] = dau1_V3D; 
      dau_V3D_[1] = dau2_V3D; 
      dau_V3D_[2] = dau3_V3D;
      daus_[0] = trackA_TLV;
      daus_[1] = trackB_TLV;
      daus_[2] = trackC_TLV;

      vx_[0] = trackA->X;
      vy_[0] = trackA->Y;
      vz_[0] = trackA->Z;
      vx_[1] = trackB->X;
      vy_[1] = trackB->Y;
      vz_[1] = trackB->Z;
      vx_[2] = trackC->X;
      vy_[2] = trackC->Y;
      vz_[2] = trackC->Z;
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
    void setIsFromTop( DecayChannel& dc, TClonesArray* branchParticle, int mcTruthIdx) {
      int isFromTop =0;
      int isFromTopIdx = dc.isFromTop( branchParticle, mcTruthIdx );
      if ( isFromTopIdx != -1 ) isFromTop_ = ((GenParticle*)branchParticle->At( isFromTopIdx))->PID;
      else isFromTop_ = isFromTop;
    }
    void setDRPT(DecayChannel& dc, TClonesArray* branchParticle, int pdgId) {
      pid_ = pdgId;
      int mcTruthIdx = dc.SearchParticle(branchParticle, pdgId, pos_ );
      if ( mcTruthIdx == -1 ) {
        dRTrue_ = 999.;
        delPtTrue_ = 1.0f;
        return ;
      }
      auto gen = (GenParticle*)branchParticle->At( mcTruthIdx );
      dRTrue_ = gen->P4().DeltaR( pos_ );
      delPtTrue_ = abs(gen->PT-pos_.Pt())/gen->PT;
      setIsFromTop( dc, branchParticle, mcTruthIdx);

    }
    void setDRPT(DecayChannel& dc, TClonesArray* branchParticle, int pid,  Track* trackA, Track* trackB=nullptr, Track* trackC=nullptr) {
      pid_ = pid;
      GenParticle *motherA=nullptr;
      bool isFound = true;
      int motherAIdx =  ((GenParticle*) trackA->Particle.GetObject())->M1;
      motherA =  (GenParticle*) branchParticle->At(motherAIdx);
      if ( motherA->PID != pid ) isFound = false;

      if ( trackB != nullptr) {
        int motherBIdx = ((GenParticle*) trackB->Particle.GetObject())->M1;
        if ( motherBIdx != motherAIdx ) isFound = false;
      }
      if ( trackC != nullptr) {
        int motherCIdx = ((GenParticle*) trackC->Particle.GetObject())->M1;
        if ( motherCIdx != motherAIdx ) isFound = false;
      }

      if ( isFound ) {
        dRTrue_= motherA->P4().DeltaR(pos_);
        delPtTrue_ = abs( motherA->PT- pos_.Pt())/motherA->PT;
      }
      else {
        dRTrue_ = 999.;
        delPtTrue_ = 1.0f;
      }
      setIsFromTop(dc, branchParticle, motherAIdx);
    }
    void setSoftLepton(int nLep) {nLep_ = nLep;}
    TLorentzVector softLeptonTLV() { 
      return softlepTLV_;
    }
    void AddSoftLepton( Track* softlepTrack) {
      softlepTLV_= TracktoTLV( softlepTrack);
      softlepPID_ = softlepTrack->PID;
      softlepPT_  = softlepTrack->PT;
    }
    Int_t softLep() { return nLep_; }
    Int_t softLepPID() { return softlepPID_;}
    Int_t softLepPT() { return softlepPT_;}

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
    Float_t dau_V3D(int idx) { return dau_V3D_[idx];}

    void reset() {
      charge_=0;
      pid_=0;
      pos_ = TLorentzVector();
      softlepTLV_ = TLorentzVector();
      nDau_ = 0;
      diffmass_=-9; 
      d0mass_=-9;
      isFromTop_=0 ;
      delPtTrue_=1;
      dRTrue_ = 999.f;
      nLep_=0;
      softlepPID_=0;
      softlepPT_ = 9.f;
      for(int i=0 ; i< 3 ; i++) {
        daus_[i]=TLorentzVector();
        dau_pid_[i]=0;
        dau_V3D_[i]=999.f;
        vx_[i]=0.f;
        vy_[i]=0.f;
        vz_[i]=0.f;
    }
  }
  std::pair< TLV, TLV> lsvPairing( QLV* lep1 , QLV* lep2) {
    int lep1_pid = lep1->first;
    int lep2_pid = lep2->first;
    LV lep1_LV = lep1->second;
    LV lep2_LV = lep2->second;

    TLorentzVector lsv1, lsv2;
    lsv1 = common::LVtoTLV(lep1_LV+common::TLVtoLV(pos_));
    lsv2 = common::LVtoTLV(lep2_LV+common::TLVtoLV(pos_));

    if (  lsv1.M() < lsv2.M() )     return make_pair(lsv1, lsv2); 
    else                        return make_pair(lsv2, lsv1); 
  }
      
  private :
    Int_t charge_, pid_;
    TLV pos_, softlepTLV_;
    Int_t nDau_;
    Float_t diffmass_, d0mass_;
    Int_t isFromTop_;
    Float_t delPtTrue_;
    Float_t dRTrue_;
    Int_t nLep_;
    Int_t softlepPID_;
    Float_t softlepPT_;

    TLV daus_[3];
    Int_t dau_pid_[3];
    Float_t dau_V3D_[3];
    Float_t vx_[3];
    Float_t vy_[3];
    Float_t vz_[3];
};    



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
    Int_t sv_softlepPID;
    Float_t sv_softlepPT;
    Float_t sv_dRTrue, sv_delPtTrue;


    Int_t sv_dau_pid[3];
    Float_t sv_dau_pt[3], sv_dau_eta[3], sv_dau_phi[3], sv_dau_mass[3], sv_dau_v3d[3];
    Int_t lsv_pid[2];
    Float_t lsv_pt[2], lsv_eta[2], lsv_phi[2], lsv_mass[2],lsv_diffmass[2], lsv_d0mass[2], llsv_mass[2];

    void reset() {
      sv_pid = 0 ;
      sv_charge =0; 
      sv_isFromTop=-9; 
      sv_pt = -9.0f; 
      sv_eta = -9.0f; 
      sv_phi = -9.0f; 
      sv_mass = -9.0f;
      sv_diffmass = -9.0f;
      sv_d0mass = -9.0f;
      sv_vx = 0.f;
      sv_vy = 0.f;
      sv_vz = 0.f;
      sv_dx = 0.f;
      sv_dy = 0.f;
      sv_dz = 0.f;
      sv_LXY= 0.f;
      sv_L3D= 0.f;
      sv_softlep=0;
      sv_softlepPID=0;
      sv_softlepPT=-9.f;
      sv_dRTrue = 999.f;
      sv_delPtTrue = 1.f;

      for( int idx = 0; idx< 3 ; idx++) { 
        sv_dau_pid[idx] =0;  
        sv_dau_pt[idx] = -9.0f; 
        sv_dau_eta[idx] = -9.0f; 
        sv_dau_phi[idx] = -9.0f; 
        sv_dau_mass[idx] = -9.0f; 
        sv_dau_v3d[idx] = 999.0f; 
      }
      for( int idx = 0; idx< 2 ; idx++) { 

        lep_charge[idx] =0;  
        lep_pid[idx] =0 ;
        lep_pt[idx] = -9.0f; 
        lep_eta[idx] = -9.0f; 
        lep_phi[idx] = -9.0f; 
        lep_mass[idx] = -9.0f; 


        lsv_pid[idx] = 0;
        lsv_pt[idx] = -9.0f;
        lsv_eta[idx] = -9.0f;
        lsv_phi[idx] = -9.0f;
        lsv_mass[idx] = -9.0f;
        lsv_diffmass[idx] = 0.0f;
        lsv_d0mass[idx] = -9.0f;
        llsv_mass[idx] = -9.0f;

      }
    }
    void fillLSV(QLV* lep1_QLV, QLV* lep2_QLV, SecVtx* svx){
      std::pair<TLorentzVector,TLorentzVector> pair = svx->lsvPairing(lep1_QLV, lep2_QLV);
      TLorentzVector lsvx1  = pair.first;
      TLorentzVector llsvx1 = pair.first  + svx->softLeptonTLV();
      TLorentzVector lsvx2  = pair.second + svx->softLeptonTLV();
      TLorentzVector llsvx2 = pair.second + svx->softLeptonTLV();

      lsv_pid[0]  = svx->pid();
      lsv_pt[0]   = lsvx1.Pt();
      lsv_eta[0]  = lsvx1.Eta();
      lsv_phi[0]  = lsvx1.Phi();
      lsv_mass[0] = lsvx1.M();
      lsv_d0mass[0] = svx->d0mass();
      lsv_diffmass[0] = svx->diffmass();
      llsv_mass[0] = llsvx1.M();

      lsv_pid[1]  = svx->pid();
      lsv_pt[1] = lsvx2.Pt();
      lsv_eta[1] = lsvx2.Eta();
      lsv_phi[1] = lsvx2.Phi();
      lsv_mass[1] = lsvx2.M();
      lsv_d0mass[1] = svx->d0mass();
      lsv_diffmass[1] = svx->diffmass();
      llsv_mass[1] = llsvx2.M();
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
      sv_diffmass = svx->diffmass();
      sv_d0mass = svx->d0mass();
      sv_vx = svx->vx();
      sv_vy = svx->vy();
      sv_vz = svx->vz();
      sv_dx = svx->vxd();
      sv_dy = svx->vyd();
      sv_dz = svx->vzd();
      sv_LXY = svx->LXY();
      sv_L3D = svx->L3D();
      sv_softlep = svx->softLep();
      sv_softlepPID = svx->softLepPID();
      sv_softlepPT = svx->softLepPT();
      sv_dRTrue = svx->DR();
      sv_delPtTrue = svx->DPT();

      for ( int idx =0 ; idx<svx->nDau() ; idx++) {
        sv_dau_pid[idx]   = svx->dau_pid(idx);
        sv_dau_pt[idx]   = svx->dau(idx).Pt();
        sv_dau_eta[idx]  = svx->dau(idx).Eta();
        sv_dau_phi[idx]  = svx->dau(idx).Phi();
        sv_dau_mass[idx] = svx->dau(idx).M();
        sv_dau_v3d[idx]  = svx->dau_V3D(idx);
      }

      lep_charge[0] = (int)-lep1_pid/abs(lep1_pid);
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


//------------------------------------------------------------------------------

class OutFileClassCustom :public OutFileClass{
public :
  OutFileClassCustom(TFile* file):OutFileClass(file) { }
  saveData data;  
};

void recoDstar(DecayChannel& dc, TClonesArray* branchParticle, Track* kaonTrack, Track* pionTrack, std::vector<Track*>& hadronTracks, SecVtx*& dstar, bool& flag_dstar, float& best_dstar_diffmass, int nsoftlep, Track* softlepTrack) { 
  TLorentzVector d0Cand = TracktoTLV(kaonTrack,321)+TracktoTLV(pionTrack,211);
  if ( abs( d0Cand.M() - 1.864 ) > d0_massWindow ) return ;
  for( unsigned int thirdTrack = 0 ; thirdTrack < hadronTracks.size() ; thirdTrack++) {
    Track* pion2Track = hadronTracks[thirdTrack];
    if ( kaonTrack == pion2Track ) continue;
    if ( pionTrack == pion2Track ) continue;
    if ( pionTrack->Charge == pion2Track->Charge ) {
      auto dstarCand = d0Cand+ TracktoTLV(pion2Track,211);
      float diffMass = dstarCand.M() - d0Cand.M(); 
      if ( diffMass > 0.135 && diffMass<0.170 && abs(0.145-diffMass)<abs(0.145-best_dstar_diffmass)) {
        float track1_v3D = dc.VertexDistance( branchParticle, kaonTrack,  pionTrack );
        float track2_v3D = dc.VertexDistance( branchParticle, kaonTrack,  pion2Track);
        float track3_v3D = dc.VertexDistance( branchParticle, pionTrack,  pion2Track);

        /*
        if ( track1_v3D > track_vertex_margin ) continue;
        if ( track2_v3D > track_vertex_margin ) continue;
        if ( track3_v3D > track_vertex_margin ) continue;
        */

        best_dstar_diffmass = diffMass;

        dstar = new SecVtx( kaonTrack, pionTrack, pion2Track, branchParticle, track1_v3D, track2_v3D, track3_v3D );
        dstar->setD0mass( d0Cand.M()); 
        dstar->setDRPT( dc, branchParticle, 413*dstar->charge(), pion2Track ); 
        dstar->setSoftLepton(nsoftlep);
        flag_dstar = true;
        if ( softlepTrack != nullptr) dstar->AddSoftLepton(softlepTrack);
      }
    } 
  }
}







//------------------------------------------------------------------------------
void AnalyseEvents(ExRootTreeReader *treeReader, OutFileClassCustom& ofc)
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
    // Event selection step1
    if ( selLepton.size()<2  || selLepton[0].first*selLepton[1].first > 0 ) continue;
    if ( (selLepton[0].second+selLepton[1].second).M()<20  ) continue;
    //std::cout<<"Pass Lepton cut"<<std::endl;
    nevt->Fill(3);
    int lepPIDmul = selLepton[0].first* selLepton[1].first;
    if ( (lepPIDmul == -121 || lepPIDmul == -169) && ((selLepton[0].second+selLepton[1].second).mass()-91.<15)  ) continue;
    nevt->Fill(4);
    std::vector<QJET> selJet = selectedJet(branchJet);
    if ( selJet.size()<2 ) continue;
    nevt->Fill(5);
    MissingET* met = (MissingET*) branchMET->At(0);
    if ( met->MET<40 ) continue;
    nevt->Fill(6);
    ofc.GetTH1("nJet")->Fill( selJet.size() );

    if ( nBJet(branchParticle, selJet) <1) continue;
    nevt->Fill(7);


    vector<Track*> leptonTracks, hadronTracks;
    SecVtx *jpsi[2] = {nullptr,nullptr};
    SecVtx *d0[2]   = {nullptr,nullptr};
    SecVtx *dstar[2] = {nullptr,nullptr};
    bool flag_jpsi[2] = {false,false};
    bool flag_d0[2] = {false,false};
    bool flag_dstar[2] = {false, false};

    if ( selJet.size() > 2 ) selJet.resize(2);
    for( unsigned int jet_idx = 0 ; jet_idx < selJet.size() ; ++jet_idx) {
      Jet* jet = selJet[jet_idx].second;
      leptonTracks.clear();
      hadronTracks.clear();
      for( unsigned int dau_idx =0 ; dau_idx < jet->Constituents.GetEntriesFast(); ++dau_idx ) {
        // temporary!!!!
        TObject* object = jet->Constituents.At(dau_idx);
        // Check if the constituent is accessible
        if(object == 0 || object->IsA() != Track::Class() ) continue;
        Track* track = (Track*) object;
        // Neutral Particles are ruled out.
        if ( track->Charge ==0 ) continue;

        // Fill Lepton Tracks
        if ( abs(track->PID) ==13      && track->PT > muonTrackPTMargin ) leptonTracks.push_back(track);
        else if ( abs( track->PID)==11 && track->PT > elecTrackPTMargin ) leptonTracks.push_back(track);
        // Fill hadron Tracks.
        else if ( track->PT > hadronPTMargin ) hadronTracks.push_back(track);
      }
      std::sort( leptonTracks.begin(), leptonTracks.end(), []( Track* a, Track* b) { return a->PT > b->PT; });
      std::sort( hadronTracks.begin(), hadronTracks.end(), []( Track* a, Track* b) { return a->PT > b->PT; });

      float best_jpsi_mass = -9.f, best_d0_mass = -9.f, best_dstar_diffmass = -9.f;

      // For J/Psi, 
      if ( leptonTracks.size() >=2 ) {
        for( unsigned int firstTrack = 0 ; firstTrack< leptonTracks.size()-1 ; firstTrack++) {
          for( unsigned int secondTrack = firstTrack+1 ; secondTrack< leptonTracks.size() ; secondTrack++) {
            int pid1 = leptonTracks[firstTrack]->PID;
            int pid2 = leptonTracks[secondTrack]->PID;

            Int_t pidMul = leptonTracks[firstTrack]->PID*leptonTracks[secondTrack]->PID;
            if ( pidMul >0 )  continue;
            // Only muon pair.
            if (  pidMul == -169 ) {
              TLorentzVector firstTrackLV = TracktoTLV(leptonTracks[firstTrack]);
              TLorentzVector secondTrackLV = TracktoTLV(leptonTracks[secondTrack]);

              TLorentzVector jpsiCand = firstTrackLV+secondTrackLV;
              if ( jpsiCand.M()>2 && jpsiCand.M()<4  &&  abs(3.096-jpsiCand.M())<abs(3.096-best_jpsi_mass) ) {
                best_jpsi_mass= jpsiCand.M();
                float track1_v3D = dc.VertexDistance( branchParticle, leptonTracks[firstTrack], leptonTracks[secondTrack]);
                jpsi[jet_idx] = new SecVtx( leptonTracks[firstTrack], leptonTracks[secondTrack], branchParticle, track1_v3D);
                jpsi[jet_idx]->setDRPT( dc, branchParticle, 443, leptonTracks[firstTrack], leptonTracks[secondTrack]);
                flag_jpsi[jet_idx] = true;
              }
            }
          }
        }
      }
      // For D0, 
      Track* softlepTrack = nullptr;
      int nsoftlep = 0;
      // D0 statement,
      if ( hadronTracks.size() <2 ) continue;
      int isD0 = 0 ;    // 0 : no softlepton, -1 : D0bar , +1 : D0

      // soft lepton( e-, mu- ) ,
      // If D0 + l^-    -> kaon- pion+ l-
      // If D0bar + l^+ -> kaon+ pion- l^+,
      if( leptonTracks.size()>0  ) { 
        softlepTrack = leptonTracks[0]; 
        if ( softlepTrack->Charge <0 ) { nsoftlep= 1 ; isD0 = 1; }
        else                           { nsoftlep= -1; isD0 = -1; }
      }

      for( unsigned int firstTrack = 0 ; firstTrack< hadronTracks.size()-1 ; firstTrack++) {
        for( unsigned int secondTrack = firstTrack+1 ; secondTrack< hadronTracks.size() ; secondTrack++) {
           
          int track1Charge = hadronTracks[firstTrack]->Charge;  
          int track2Charge = hadronTracks[secondTrack]->Charge;
          if ( track1Charge*track2Charge >0 ) continue;

          int recoPID = 421*isD0;

          Track *kaonTrack=nullptr, *pionTrack=nullptr, *pion2Track=nullptr;
          if ( (isD0==1 && track1Charge <0) || (isD0 == -1 && track1Charge>0) ) { kaonTrack = hadronTracks[firstTrack]; pionTrack = hadronTracks[secondTrack]; }
          else if ( (isD0==1 && track1Charge>0) || (isD0 ==-1 && track1Charge<0) ) { kaonTrack = hadronTracks[secondTrack]; pionTrack = hadronTracks[firstTrack]; } 
          if ( isD0 != 0 ) {
            TLorentzVector d0Cand = TracktoTLV(kaonTrack,321)+TracktoTLV(pionTrack,211);
            if ( d0Cand.M() > 1 && d0Cand.M()<3 && (abs(1.864 - d0Cand.M()) < abs(1.864-best_d0_mass))  ) {
              float track1_v3D = dc.VertexDistance( branchParticle, kaonTrack, pionTrack );
              if ( track1_v3D > track_vertex_margin ) continue;
              best_d0_mass = d0Cand.M();
              d0[jet_idx] = new SecVtx( kaonTrack, pionTrack, branchParticle, track1_v3D );
              d0[jet_idx]->setDRPT( dc, branchParticle, recoPID, kaonTrack, pionTrack); 
              d0[jet_idx]->setSoftLepton(nsoftlep);
              d0[jet_idx]->AddSoftLepton(softlepTrack);
              d0[jet_idx]->setD0mass( d0Cand.M() ); 
              flag_d0[jet_idx] = true;

              // For D*,
              kaonTrack = hadronTracks[firstTrack]; pionTrack = hadronTracks[secondTrack];
              recoDstar(dc, branchParticle, kaonTrack, pionTrack, hadronTracks, dstar[jet_idx], flag_dstar[jet_idx], best_dstar_diffmass , nsoftlep, softlepTrack);
            }
          }
        }
      }
    }
    // For each hemisphere, 
    for(int i=0 ; i<2 ; i++) { 
      int cat = 0; 
      if ( flag_jpsi[i] ) {
        ofc.GetTH1("jpsi_mass")->Fill( jpsi[i]->mass() );
        ofc.data.init(&selLepton[0], &selLepton[1], jpsi[i]);
        if ( abs( 3.096-jpsi[i]->mass())< jpsi_massWindow) { 
          nJpsi++;
          cat = 2;
          ofc.data.fillLSV(&selLepton[0], &selLepton[1], jpsi[i]);
          ofc.GetTH1("sv_category")->Fill(cat); 
        }
        ofc.GetTree("tree")->Fill(); 
      }
      else if ( flag_dstar[i]) {
        ofc.GetTH1("dstar_mass")->Fill( dstar[i]->mass() );
        ofc.GetTH1("dstar_diffmass")->Fill( dstar[i]->diffmass() );
        ofc.data.init(&selLepton[0], &selLepton[1], dstar[i]);
        if ( abs( 0.145-dstar[i]->diffmass())< dstar_massWindow ) {
          nDstar++;
          cat = 3;
          ofc.data.fillLSV(&selLepton[0], &selLepton[1], dstar[i]);
          ofc.GetTH1("sv_category")->Fill(cat); 
        } 
        ofc.GetTree("tree")->Fill(); 
      }
      else if ( flag_d0[i]) {
        ofc.GetTH1("d0_mass")->Fill( d0[i]->mass() );
        ofc.data.init(&selLepton[0], &selLepton[1], d0[i]);
        if ( abs( 1.864-d0[i]->mass())< d0_massWindow ) {
          nD0++;
          cat = 4;
          ofc.data.fillLSV(&selLepton[0], &selLepton[1], d0[i]);
          ofc.GetTH1("sv_category")->Fill(cat); 
        } 
        ofc.GetTree("tree")->Fill();
      }
      if ( !flag_jpsi[i]&&!flag_d0[i]&& !flag_dstar[i]) {
        cat = 1; 
        ofc.GetTH1("sv_category")->Fill(cat); 
      }
    }
  }        
  std::cout<<TString::Format("Filling! nJpsi: %d , nD0 : %d, nDstar : %d\n",nJpsi,nD0, nDstar)<<std::endl;
}

//------------------------------------------------------------------------------


//------------------------------------------------------------------------------


void BookingHist(OutFileClassCustom& ofc)
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
  const char* inputFile ;
  const char* outputFile;

  TChain *chain = new TChain("Delphes");
  if ( argc ==1 ) {
    for( int i=10 ; i<30 ; i++) {
      TString input = TString::Format("./Delphes_%d.root",i).Data();
      //inputFile = std::string("Delphes_")+TString::Format("Delphes_%d.root",i).Data();
      //sprintf(inputFile, "Delphes_%d.root",i);
      chain->Add( input.Data()); 
    } 
    outputFile = "result_output.root";
  }
  else if ( argc !=3 ) {
    std::cout<<"Argument is wrong. Please, run \"JetCharge Delphes_input.root result_output.root\""<<std::endl;
    return -1;
  }
  else {
    inputFile = argv[1];
    outputFile = argv[2];
    chain->Add(inputFile);
  }

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  TFile* file = TFile::Open(outputFile,"RECREATE"); 
  auto ofc = OutFileClassCustom(file);
  BookingHist(ofc);
  TString brC = TString("");
  // Common
  //brC += "dataType/I:";
  // Lepton
  brC += "lep_charge[2]/I:lep_pid[2]/I:lep_pt[2]/F:lep_eta[2]/F:lep_phi[2]/F:lep_mass[2]/F:";
  // SV
  brC += "sv_pid/I:sv_charge/I:sv_isFromTop/I:sv_pt/F:sv_eta/F:sv_phi/F:sv_mass/F:sv_diffmass/F:sv_d0mass/F:sv_vx/F:sv_vy/F:sv_vz/F:sv_dx/F:sv_dy/F:sv_dz/F:sv_L3D/F:sv_LXY/F:";
  brC += "sv_nsoftlep/I:sv_softlepPID/I:sv_softlepPT/F:sv_dRTrue/F:sv_delPtTrue/F:sv_dau_pid[3]/I:sv_dau_pt[3]/F:sv_dau_eta[3]/F:sv_dau_phi[3]/F:sv_dau_mass[3]/F:sv_dau_v3d[3]/F";
  // lSV
  TString brC2 = TString("");
  brC2 += "lsv_pid[2]/I:lsv_pt[2]/F:lsv_eta[2]/F:lsv_phi[2]/F:lsv_mass[2]/F:lsv_diffmass[2]/F:lsv_d0mass[2]/F:llsv_mass[2]/F";

  std::cout<<brC<<std::endl;
  ofc.BookingTree("SVTree", &ofc.data, brC.Data());
  ofc.addBranch( "LSVTree", &ofc.data.lsv_pid[0], brC2.Data());                      
  AnalyseEvents(treeReader,ofc);

  ofc.WriteAll();
  cout << "** Exiting..." << endl;

  file->Close();
  delete treeReader;
  delete chain;

}
