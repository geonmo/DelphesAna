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

int DecayChannel::FindJetParton(TClonesArray* genParticles, int baseIdx) {
  if ( baseIdx <0 ) return -1;
  GenParticle* base = (GenParticle*)genParticles->At(baseIdx);
  int absPID = abs(base->PID);
  std::cout<<"abs PID : "<<absPID<<std::endl;
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
  GenParticle* particle ;
  for ( unsigned int i=0 ; i< jet->Constituents.GetEntriesFast() ; i++ ) {
    TObject* object = jet->Constituents.At(i) ; 
    if ( object ==0 ) continue;
    if ( object->IsA() == GenParticle::Class()) {
      particle = (GenParticle*)object;
      return value;
    }
    else if ( object->IsA() == Track::Class()) {
      Track* track = (Track*) object;
      particle = (GenParticle*)track->Particle.GetObject();
    }
    /*
    else if ( object->IsA() == Tower::Class()) {
      Tower* tower = (Tower*) object;
      particle = (GenParticle*)tower->Particles.At(0);
    }
    */
    else continue; 
    value = FindJetParton(genParticles,  particle->M1);
    return value;
  } 
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


