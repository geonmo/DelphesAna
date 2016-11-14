#ifndef __DELPHESANA__EVENTSELECTOR__
#define __DELPHESANA__EVENTSELECTOR__
#include "classes/DelphesClasses.h"
#include <TMath.h>
#include <TLorentzVector.h>

std::vector<std::pair<int, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > selectedLepton( TClonesArray* muons, TClonesArray* electrons) {
  std::vector<std::pair<int, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > > selLeptons;
  for( int i= 0 ; i< muons->GetEntriesFast() ; ++i) {
    Muon* muon = (Muon*) muons->At(i);
    if ( muon->PT<20 || abs(muon->Eta)>2.4) continue;
    // pid == 13, Q = -1
    int pid = muon->Charge* -13;
    TLorentzVector muonP4 = muon->P4();
    TLorentzVector newMuonP4 = TLorentzVector();
    newMuonP4.SetPtEtaPhiM( muonP4.Pt(), muonP4.Eta(), muonP4.Phi(), 0.105652);
    selLeptons.push_back( make_pair(pid, common::TLVtoLV ( newMuonP4 ) ));
  }
  for( int i= 0 ; i< electrons->GetEntriesFast() ; ++i) {
    Electron* electron = (Electron*) electrons->At(i);
    if ( electron->PT<20 || abs(electron->Eta)>2.4) continue;
    // pid : 11 = Q:-1 * -11 
    int pid = electron->Charge* -11;
    TLorentzVector elecP4 = electron->P4();
    TLorentzVector newElecP4 = TLorentzVector();
    newElecP4.SetPtEtaPhiM( elecP4.Pt(), elecP4.Eta(), elecP4.Phi(), 0.000511169 );
    selLeptons.push_back( make_pair(pid, common::TLVtoLV ( newElecP4) ));
  }
  std::sort( selLeptons.begin(), selLeptons.end(), []( std::pair<int, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > a, std::pair<int, ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > > b) { return a.second.pt() > b.second.pt(); });
  return selLeptons; 
}
std::vector<std::pair<int, Jet*>> selectedJet( TClonesArray* jets ) {
  std::vector<std::pair<int, Jet*>> selJets;
  for(int i = 0; i < jets->GetEntriesFast(); ++i)
  {
    Jet* jet = (Jet*) jets->At(i);
    if ( jet->PT < 30 || abs(jet->Eta)>2.4 ) continue;
    selJets.push_back( make_pair(jet->Charge, jet) );
  }
  std::sort( selJets.begin(), selJets.end(), []( std::pair<int, Jet*> a, std::pair<int, Jet*> b) { return a.second->PT > b.second->PT; });
  return selJets; 
}

int nBJet( TClonesArray* branchParticle, std::vector<std::pair<int, Jet*>> jets) {
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
#endif
