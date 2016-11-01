#ifndef __Geonmo_DecayChannel__
#define __Geonmo_DecayChannel__

#include <TClonesArray.h>
using namespace std;


class DecayChannel {

  private :
    int channel_;
    TClonesArray* genParticles_;
  public :
    DecayChannel( TClonesArray* genParticles);
    static int FindJetParton(TClonesArray*,Jet*);
    static int FindJetParton(TClonesArray*,int);
    int isFromTop(TClonesArray* genParticles, int);
    int isFromB(TClonesArray* genParticles, int);
    int SearchParticle(TClonesArray* genParticles, int pid, TLorentzVector cand);
    GenParticle* SearchParticleRef(TClonesArray* genParticles, int pid, TLorentzVector cand);
    GenParticle* SearchParticleRef(TClonesArray* genParticles, int pid, Track* cand) { return SearchParticleRef(genParticles, pid, cand->P4());}
    float VertexDistance(TClonesArray* , GenParticle* gen, Track* cand);
    float VertexDistance(TClonesArray* , Track* cand1, Track* cand2);
    int FindWboson(TClonesArray* , int baseIdx );
    int FindLepton(TClonesArray* , int baseIdx );
    int channelSelection(TClonesArray*  ); 
    int channel() { return channel_;}
};

#endif
