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


#ifndef Run_h
#define Run_h 1

#include "G4Run.hh"
#include "G4Event.hh"
#include "G4THitsMap.hh"
#include "G4SDManager.hh"
#include "TETModelImport.hh"

using namespace std;

typedef map<G4int, std::pair<G4double, G4double>> HITSMAP;

class Run : public G4Run
{
public:
	Run(TETModelImport* _tetData_doctor);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

	HITSMAP* GetdoctorEdepMap()  { return &doctor_edepMap; }

private:
	void StackHits(G4HCofThisEvent* HCE, G4int fCollID_DRF, G4int fCollID_eDep ,G4bool doseOrganized, 
				   map<G4int, G4double> rbmFactor, map<G4int, G4double> bsFactor, map<G4int, vector<G4int>> organ2dose, 
				   HITSMAP &edepMap);


private:
	// TETModelImport* tetData_patient;
	// G4int fCollID_patient_eDep;
	// G4int fCollID_patient_DRF;
	// map<G4int, vector<G4int>>  patient_organ2dose;
	// map<G4int, G4double>  patient_rbmFactor;
	// map<G4int, G4double>  patient_bsFactor;
	// G4bool patient_doseOrganized;
	// HITSMAP patient_edepMap;

	TETModelImport* tetData_doctor;
	G4int fCollID_doctor_eDep;
	G4int fCollID_doctor_DRF;
	map<G4int, vector<G4int>>  doctor_organ2dose;
	map<G4int, G4double>  doctor_rbmFactor;
	map<G4int, G4double>  doctor_bsFactor;
	G4bool doctor_doseOrganized;
	HITSMAP doctor_edepMap;
};

#endif
