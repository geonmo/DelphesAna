#include <TClonesArray.h>
using namespace std;


class DecayChannel {

  private :
    int channel_;
  public :
    DecayChannel( TClonesArray* genParticles);
    int FindWboson(TClonesArray* genParticles, int baseIdx );
    int FindLepton(TClonesArray* genParticles, int baseIdx );
    int channelSelection( TClonesArray* genParticles); 
};

