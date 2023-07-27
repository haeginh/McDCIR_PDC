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

#include "PrimaryGeneratorAction_PS.hh"
#include "G4ParticleTable.hh"
#include <fstream>

PrimaryGeneratorAction_PS::PrimaryGeneratorAction_PS()
	: G4VUserPrimaryGeneratorAction(), nThRotation(0)
{
	fPrimary = new G4ParticleGun();
	fPrimary->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));

	messenger = new PrimaryMessenger_PS(this);
}

PrimaryGeneratorAction_PS::~PrimaryGeneratorAction_PS()
{
	delete fPrimary;
	delete messenger;
}

void PrimaryGeneratorAction_PS::GeneratePrimaries(G4Event *anEvent)
{
	if(anEvent->GetEventID()>(nThRotation+1)*PSnum)
	{
		nThRotation++;
	}
	G4int start = (anEvent->GetEventID()-nThRotation*PSnum)*7;
	fPrimary->SetParticleEnergy(data[start]);
	fPrimary->SetParticlePosition(G4ThreeVector(data[start+1], data[start+2], data[start+3]));
	fPrimary->SetParticleMomentumDirection(G4ThreeVector(data[start+4], data[start+5], data[start+6]));
	fPrimary->GeneratePrimaryVertex(anEvent);
}

void PrimaryGeneratorAction_PS::SetPS(G4String name)
{
	data.clear();
	std::ifstream ifs(name, std::ios::binary);
	G4int num(0);
	ifs.read((char*) &num, sizeof(G4int));
	data.resize(num*7);
	ifs.read((char*) &data[0], num*7*sizeof(G4double));
	G4cout<<"read phase-space file for "<<num<<" particles"<<G4endl;
	ifs.close();
}
