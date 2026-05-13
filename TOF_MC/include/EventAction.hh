#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4Event;

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event) override;
    
    virtual void EndOfEventAction(const G4Event* event) override;

    void SetTimeBar1(G4double t){fTimeBar1 = t;}

    G4double GetTimeBar1(){return fTimeBar1;}

    void SetTimeBar2(G4double t){fTimeBar2 = t;}

    G4double GetTimeBar2(){return fTimeBar2;}

    void SetTimeBarSmall(G4double t){fTimeBarSmall = t;}

    G4double GetTimeBarSmall(){return fTimeBarSmall;}

    void SetTimeVeto(G4double t){fTimeVeto = t;}

    G4double GetTimeVeto(){return fTimeVeto;}

    void SetHitBarPos(G4ThreeVector pos)  { fHitBarPos = pos; }

    void SetHitSurfPos(G4ThreeVector pos) { fHitSurfPos = pos; }

private:

    G4double fTimeBar1;

    G4double fTimeBar2;

    G4double fTimeBarSmall;

    G4double fTimeVeto;

    G4ThreeVector fHitBarPos;
    
    G4ThreeVector fHitSurfPos;

    G4double fEdepPb = 0.;

    G4double fTrackLenPb = 0.;
};

#endif