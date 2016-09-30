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
    int FindWboson(TClonesArray* , int baseIdx );
    int FindLepton(TClonesArray* , int baseIdx );
    int channelSelection(TClonesArray*  ); 
    int channel() { return channel_;}
};

#endif
