#ifndef NDEBUG
#define D(x) x 
#else
#define D(x)  
#endif

#include<iostream>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#include <TTree.h>
#include "DecayChannel.h"
using namespace std;

class ExRootTreeReader;

DecayChannel::DecayChannel( TClonesArray* genParticles){
  channel_ = channelSelection(genParticles ) ;   

}

int DecayChannel::isFromTop(TClonesArray* genParticles, int baseIdx ) {
  if ( baseIdx <0 ) return -1;
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 6 ) return baseIdx;

  // First, M1
  if ( base->M1 >=0  ) {
    int nextResult = isFromTop( genParticles, base->M1 );
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}
int DecayChannel::isFromB(TClonesArray* genParticles, int baseIdx ) {
  if ( baseIdx <0 ) return -1;
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 5 ) return baseIdx;

  // First, M1
  if ( base->M1 >= 0 ) {
    int nextResult = isFromTop( genParticles, base->M1 );
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}

int DecayChannel::FindJetParton(TClonesArray* genParticles, int baseIdx) {
  if ( baseIdx <0 ) return -1;
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID < 6 ) return baseIdx;

  // First, M1
  if ( base->M1 !=-1 ) {
    int nextResult = FindJetParton( genParticles, base->M1);
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;

}

int DecayChannel::FindJetParton(TClonesArray* genParticles, Jet* jet) {
  int value = -1;
  std::vector<GenParticle*> particles ;
  std::vector<Track*> tracks ;
  for ( unsigned int i=0 ; i< jet->Constituents.GetEntriesFast() ; i++ ) {
    TObject* object = jet->Constituents.At(i) ; 
    if ( object ==0 ) continue;
    if ( object->IsA() == GenParticle::Class()) {
      GenParticle* particle = (GenParticle*)object;
      if ( particle->Charge==0 ) continue;
      particles.push_back(particle);
    }
    else if ( object->IsA() == Track::Class()) {
      Track* track = (Track*) object;
      tracks.push_back(track);
    }
  } 
  GenParticle* particle ; 
  if ( tracks.size() >0 ) { 
    std::sort( tracks.begin(), tracks.end(), []( Track* a, Track* b) { return a->PT > b->PT; });
    particle = (GenParticle*)tracks[0]->Particle.GetObject();
    value = FindJetParton(genParticles,  particle->M1);
  }
  if ( particles.size()>0 ) {
    std::sort( particles.begin(), particles.end(), []( GenParticle* a, GenParticle* b) { return a->PT > b->PT; });
    particle = particles[0];
    value = FindJetParton(genParticles,  particle->M1); 
  }
  return value;
}


int DecayChannel::SearchParticle(TClonesArray* genParticles, int pid, TLorentzVector cand) {
  float dR = 999.f;
  float delPt = 1.0f;
  int particle_idx = -1;
  bool found = false;
  for( int i=0 ; i< genParticles->GetEntriesFast(); i++) {
    GenParticle* gen = (GenParticle*)genParticles->At(i);
    if ( gen->PID != pid) continue;
    int top_idx = isFromTop(genParticles, i);
    int b_idx = isFromB(genParticles, i);
    if ( top_idx ==-1 || b_idx == -1 ) continue;

    float dR_current = gen->P4().DeltaR( cand );
    float delPt_current = abs(gen->PT-cand.Pt())/gen->PT;
    if ( dR_current<dR && delPt_current< delPt ) { dR = dR_current; delPt = delPt_current; particle_idx = i; }
  } 
  return particle_idx;
}
GenParticle* DecayChannel::SearchParticleRef(TClonesArray* genParticles, int pid, TLorentzVector cand) {
  float dR = 999.f;
  float delPt = 1.0f;
  int particle_idx = -1;
  bool found = false;
  for( int i=0 ; i< genParticles->GetEntriesFast(); i++) {
    GenParticle* gen = (GenParticle*)genParticles->At(i);
    if ( gen->PID != pid) continue;
    int top_idx = isFromTop(genParticles, i);
    int b_idx = isFromB(genParticles, i);
    if ( top_idx ==-1 || b_idx == -1 ) continue;

    float dR_current = gen->P4().DeltaR( cand );
    float delPt_current = abs(gen->PT-cand.Pt())/gen->PT;
    if ( dR_current<dR && delPt_current< delPt ) { dR = dR_current; delPt = delPt_current; particle_idx = i; }
  }
  if ( particle_idx != -1) return (GenParticle*) genParticles->At(particle_idx);
  else return nullptr;
}
float DecayChannel::VertexDistance(TClonesArray* genParticles, GenParticle* gen, Track* cand) {

  if ( gen == nullptr || cand == nullptr) return 999.f;
  TVector3 genV3( gen->X, gen->Y, gen->Z);
  TVector3 candV3(cand->X, cand->Y, cand->Z);
  GenParticle* candMother=nullptr;
  if ( candV3.Mag() ==0 ) {
    int candMotherIdx = ((GenParticle*)cand->Particle.GetObject())->M1;
    if ( candMotherIdx >= 0 ) candMother = (GenParticle*)genParticles->At( candMotherIdx );
    candV3.SetXYZ( candMother->X, candMother->Y, candMother->Z);
  }
 
  D(cout << "    Track pt: " <<    cand->PT << ", eta: " <<    cand->Eta << ", phi: " <<    cand->Phi << " , mass: "<<cand->P4().M()<<" ,PID: "<<   cand->PID<<endl;)
  D(cout << "    Gen   pt: " << gen->PT << ", eta: " << gen->Eta << ", phi: " << gen->Phi << " , mass: "<<gen->Mass<<" , PID: "<<gen->PID<<endl;)
  D(if ( candMotherIdx>= 0 ) cout << "    Mother  pt: " << candMother->PT << ", eta: " << candMother->Eta << ", phi: " << candMother->Phi << " , mass: "<<candMother->Mass<<" , PID: "<<candMother->PID<<endl;)
  return genV3.DeltaR(candV3);
}
float DecayChannel::VertexDistance(TClonesArray* genParticles, Track* track1, Track* track2) {
  TVector3 track1V3( track1->X, track1->Y, track1->Z);
  TVector3 track2V3( track2->X, track2->Y, track2->Z);
  GenParticle *track1Mother=nullptr, *track2Mother=nullptr;

  if ( track1V3.Mag() ==0 ) {
    int track1MotherIdx = ((GenParticle*)track1->Particle.GetObject())->M1;
    if ( track1MotherIdx >= 0 ) { 
      track1Mother = (GenParticle*)genParticles->At( track1MotherIdx );
      track1V3.SetXYZ( track1Mother->X, track1Mother->Y, track1Mother->Z);
    }
  }
  if ( track2V3.Mag() ==0 ) {
    int track2MotherIdx = ((GenParticle*)track2->Particle.GetObject())->M1;
    if ( track2MotherIdx >= 0 ) { 
      track2Mother = (GenParticle*)genParticles->At( track2MotherIdx );
      track2V3.SetXYZ( track2Mother->X, track2Mother->Y, track2Mother->Z);
    }
  }
 
  return track1V3.DeltaR(track2V3);
}

int DecayChannel::FindWboson(TClonesArray* genParticles, int baseIdx )
{
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 24 ) return baseIdx;

  // First, D1
  if ( base->D1 !=-1 ) {
    int nextResult = FindWboson( genParticles, base->D1);
    if ( nextResult !=-1 ) return nextResult;
  }
  // Next, D2
  if ( base->D2 !=-1 ) {
    int nextResult = FindWboson( genParticles, base->D2);
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}
int DecayChannel::FindLepton(TClonesArray* genParticles, int baseIdx )
{
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  if ( absPID == 11 || absPID == 13 || absPID == 15 ) return baseIdx;

  // First, D1
  if ( base->D1 !=-1 ) {
    int nextResult = FindLepton( genParticles, base->D1);
    if ( nextResult !=-1 ) return nextResult;
  }
  // Next, D2
  if ( base->D2 !=-1 ) {
    int nextResult = FindLepton( genParticles, base->D2);
    if ( nextResult !=-1 ) return nextResult;
  }
  return -1;
}
int DecayChannel::channelSelection( TClonesArray* genParticles  ) {
  // searcing Top or anti top
  int nGenParticle = genParticles->GetEntriesFast();
  D( std::cout<<"Num of genParticles : "<<nGenParticle<<std::endl; )
    int top_idx =-1, antitop_idx=-1;

  for (int  i =0 ; i < nGenParticle ; ++i) {
    GenParticle* genParticle = (GenParticle*) genParticles->At(i);
    D( std::cout<<"genParticle PID : "<<genParticle->PID<<std::endl; )
      if ( top_idx==-1 && genParticle->PID == 6) top_idx= i;
    if ( antitop_idx==-1 && genParticle->PID == -6 ) antitop_idx= i;
    if ( top_idx !=-1 && antitop_idx != -1 ) break;
  }

  if ( top_idx == -1 || antitop_idx == -1 ) {
    D( std::cout<<"Can not find top quark for this events. \nPlease, check your data which include the partons(status ==3 or 13 for pythia6,8)"<<std::endl;)
    D( std::cout<<"Skip decay channel selection."<<std::endl;)
    return 0; 
  }

  int WbosonIdx1 = FindWboson( genParticles,     top_idx);
  int WbosonIdx2 = FindWboson( genParticles, antitop_idx);

  bool     topToLepton = false;
  bool antitopToLepton = false;

  int lep1Idx = FindLepton( genParticles, WbosonIdx1);
  int lep2Idx = FindLepton( genParticles, WbosonIdx2);
  int lep1=-1, lep2=-1;
  float lep1PT=0.0, lep2PT=0.0;
  if ( lep1Idx != -1) { 
    GenParticle* lep1Part = ((GenParticle*) genParticles->At( lep1Idx) );
    lep1 = lep1Part->PID;
    lep1PT = lep1Part->PT;
  }
  if ( lep2Idx != -1) {
    GenParticle* lep2Part = ((GenParticle*) genParticles->At( lep2Idx) );
    lep2 = lep2Part->PID;
    lep2PT = lep2Part->PT;
  }

  D(std::cout<<std::endl;)
    int trackingIdx;
  if ( lep1Idx != -1 ) {
    trackingIdx = lep1Idx;
    int upper_count = 0;
    while( 1 ) {
      GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
      D(std::cout<< base->PID;)
        if ( abs(base->PID) == 24 ) {break;}
      D(std::cout<<">>";)
        trackingIdx = base->M1;
      upper_count++;
    }
    if ( upper_count >1 ) { 
      D(std::cout<<"oops"<<std::endl;) 
        lep1 = -1; 
      lep1Idx = -1; 
    }
  } 
  D(std::cout<<std::endl;)
    if ( lep2Idx != -1 ) {
      trackingIdx = lep2Idx;
      int upper_count = 0;
      while( 1 ) {
        GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
        D(std::cout<< base->PID;)
          if ( abs(base->PID) == 24 ) {break;}
        D(std::cout<<">>";)
          trackingIdx = base->M1;
        upper_count++;
      }
      if ( upper_count >1 ) { 
        D(std::cout<<"oops"<<std::endl; )
          lep2 = -1; 
        lep2Idx = -1; 
      }
    } 
  D(std::cout<<std::endl;)

    int mulValue = lep1*lep2;

  if ( lep1Idx != -1 ) {
    D(std::cout<<"Lep1 PT : "<<lep1PT<<std::endl; )
  }
  if ( lep2Idx != -1 ) {
    D(std::cout<<"Lep2 PT : "<<lep2PT<<std::endl; )
  }

  if ( abs(mulValue) > 100 ) { 
    if ( abs(mulValue) %15 ==0 ) {
      D(std::cout<<"Dilepton tau"<<std::endl;)
        return 3;
    }
    D(std::cout<<"Dilepton"<<std::endl;)
      D(std::cout<<"lep1 :"<<lep1<<"  lep2 : "<<lep2<<std::endl;)
      return 0;
  }
  else if ( abs(mulValue) > 10 ) {
    if ( abs(mulValue) %15 ==0 ) {
      D(std::cout<<"Semilepton tau"<<std::endl;)
        return 4;
    }
    D(std::cout<<"Semilepton"<<std::endl;)
      return 1;
  }
  else {
    D(std::cout<<"Hardronic"<<std::endl;)
      return 2; 
  }
  return -1;

}


