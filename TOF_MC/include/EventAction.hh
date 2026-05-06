#ifndef EventAction_h
#define EventAction_h

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;

class EventAction : public G4UserEventAction
{
public:
    EventAction();
    virtual ~EventAction();

    virtual void BeginOfEventAction(const G4Event* event) override;
    virtual void EndOfEventAction(const G4Event* event) override;

    void SetTimeBar1(G4double t)      { fTimeBar1 = t; }
    void SetTimeBar2(G4double t)      { fTimeBar2 = t; }
    void SetTimeBarSmall(G4double t)  { fTimeBarSmall = t; }

private:
    G4double fTimeBar1;
    G4double fTimeBar2;
    G4double fTimeBarSmall;
};

#endif