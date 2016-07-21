#ifndef NDEBUG
#define D(x) 
#else
#define D(x) x
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
  channel_ = channelSelection( genParticles) ;    
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
int DecayChannel::channelSelection( TClonesArray* genParticles) {
  // searcing Top or anti top
  int nGenParticle = genParticles->GetEntriesFast();

  int top_idx =0, antitop_idx=0;

  for (int  i =0 ; i < genParticles->GetEntriesFast() ; ++i) {
    GenParticle* genParticle = (GenParticle*) genParticles->At(i);
    if ( top_idx==0 && genParticle->PID == 6) top_idx= i;
    if ( antitop_idx==0 && genParticle->PID == -6 ) antitop_idx= i;
    if ( top_idx !=0 && antitop_idx != 0 ) break;
  }

  int WbosonIdx1 = FindWboson(genParticles,     top_idx);
  int WbosonIdx2 = FindWboson(genParticles, antitop_idx);

  bool     topToLepton = false;
  bool antitopToLepton = false;

  int lep1Idx = FindLepton( genParticles, WbosonIdx1);
  int lep2Idx = FindLepton( genParticles, WbosonIdx2);
  int lep1=-1, lep2=-1;
  if ( lep1Idx != -1) lep1 = ((GenParticle*) genParticles->At( lep1Idx) )->PID;
  if ( lep2Idx != -1) lep2 = ((GenParticle*) genParticles->At( lep2Idx) )->PID;

  D(std::cout<<std::endl;)
    int trackingIdx;
  if ( lep1Idx != -1 ) {
    trackingIdx = lep1Idx;
    int upper_count = 0;
    while( 1 ) {
      GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
      D(std::cout<< base->PID;)
        if ( abs(base->PID) == 24 ) break;
      D(std::cout<<">>";)
        trackingIdx = base->M1;
      upper_count++;
    }
    if ( upper_count >1 ) { 
      D(std::cout<<"oops"<<std::endl;) 
        lep1 = -1; 
    }
  } 
  D(std::cout<<std::endl;)
    if ( lep2Idx != -1 ) {
      trackingIdx = lep2Idx;
      int upper_count = 0;
      while( 1 ) {
        GenParticle* base = (GenParticle*) genParticles->At( trackingIdx);
        D(std::cout<< base->PID;)
          if ( abs(base->PID) == 24 ) break;
        D(std::cout<<">>";)
          trackingIdx = base->M1;
        upper_count++;
      }
      if ( upper_count >1 ) { 
        D(std::cout<<"oops"<<std::endl; )
          lep2 = -1; 
      }
    } 
  D(std::cout<<std::endl;)

    int mulValue = lep1*lep2;

  if ( abs(mulValue) > 100 ) { 
    if ( abs(mulValue) %15 ==0 ) {
      D(std::cout<<"Dilepton tau"<<std::endl;)
        return 4;
    }
    D(std::cout<<"Dilepton"<<std::endl;)
      D(std::cout<<"lep1 :"<<lep1<<"  lep2 : "<<lep2<<std::endl;)
      return 1;
  }
  else if ( abs(mulValue) > 10 ) {
    if ( abs(mulValue) %15 ==0 ) {
      D(std::cout<<"Semilepton tau"<<std::endl;)
        return 5;
    }
    D(std::cout<<"Semilepton"<<std::endl;)
      return 2;
  }
  else {
    D(std::cout<<"Hardronic"<<std::endl;)
      return 3; 
  }

}


