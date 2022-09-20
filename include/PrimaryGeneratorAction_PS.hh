//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//

#ifndef PrimaryGeneratorAction_PS_hh
#define PrimaryGeneratorAction_PS_hh 1

using namespace std;

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4UImessenger.hh"

class PrimaryGeneratorAction_PS : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction_PS();
  virtual ~PrimaryGeneratorAction_PS();

  virtual void GeneratePrimaries(G4Event *);

  void SetPS(G4String name); //peakE in keV

private:
  G4ParticleGun *fPrimary;
  G4UImessenger *messenger;
  G4int nThRotation, PSnum;
  vector<G4double> data;
};

#include "G4UIcmdWithAString.hh"
class PrimaryMessenger_PS : public G4UImessenger
{
public:
  PrimaryMessenger_PS(PrimaryGeneratorAction_PS *_primary)
      : G4UImessenger(), fPrimary(_primary)
  {
    fBeamDir = new G4UIdirectory("/beam/");
    fPSCmd = new G4UIcmdWithAString("/beam/PS", this);
  };
  virtual ~PrimaryMessenger_PS()
  {
    delete fPSCmd;
    delete fBeamDir;
  };
  virtual void SetNewValue(G4UIcommand *command, G4String newValue)
  {
    if (command == fPSCmd)
    {
      fPrimary->SetPS(newValue);
    }
  };

private:
  PrimaryGeneratorAction_PS *fPrimary;
  G4UIdirectory *fBeamDir;
  G4UIcmdWithAString *fPSCmd;
};

#endif
