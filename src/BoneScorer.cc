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

#include "BoneScorer.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4RunManager.hh"
//using namespace std;
///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring cell charge.
//   The Cell Charge is defined by  a sum of charge inside the cell
//  which calculates the deposited charge in the cell.
///////////////////////////////////////////////////////////////////////////////

//Unit: Gy
BoneScorer::BoneScorer(G4String name,ParallelPhantom* _phantom)
  :G4VPrimitiveScorer(name), phantom(_phantom), HCID(-1), EvtMap(0)
{
	// energyBin={0.01 ,0.015 ,0.02 ,0.03 ,0.04 ,0.05 ,0.06 ,0.08
	// 			,0.1 ,0.15 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.8
	// 			,1 ,1.5 ,2 ,3 ,4 ,5 ,6 ,8 ,10};
	// RBMratio = phantom->GetRBMratio();
	// BSratio = phantom->GetBSratio();
	gamma = G4ParticleTable::GetParticleTable()->FindParticle("gamma");
}

BoneScorer::~BoneScorer()
{
	EvtMap->clear();
}

G4int BoneScorer::GetIndex(G4Step* aStep)
{
	G4int copyNo = aStep->GetPreStepPoint()->GetTouchable()->GetCopyNumber();
	return phantom->GetOrganID(copyNo);
}

G4bool BoneScorer::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (aStep -> GetTrack() -> GetDefinition()!=gamma) return FALSE;

	G4double stepLength = aStep->GetStepLength();
	if (stepLength==0.) return FALSE;

	G4int index = GetIndex(aStep);
	G4double rbmR, bsR;
	if(!phantom->GetBoneMassRatio(index, rbmR, bsR)) return FALSE;
	G4double edep = aStep->GetTotalEnergyDeposit();
	// G4double massInv = 1./(phantom->GetVolume(index)*phantom->GetMateiral(index)->GetDensity());
	// G4double rbmDose(edep * rbmR), bsDose(edep * bsR);
	// EvtMap->add(2, rbmDose); //rbm simple mass R
	// EvtMap->add(3, bsDose); // bs simple mass R

	G4double CellFlux = stepLength / phantom->GetVolume(index);
	G4double energy=aStep->GetPreStepPoint()->GetKineticEnergy();
	
	// G4double RBMdose = GetRBMdose(energy, CellFlux, index);
	// G4double BSdose = GetBSdose(energy, CellFlux, index);
	G4double rbmDrfDose, bsDrfDose;
	tie(rbmDrfDose, bsDrfDose) = GetDRFdose(energy, CellFlux, index);

	if(std::isnan(rbmDrfDose)) rbmDrfDose=0;
	if(std::isnan(bsDrfDose)) bsDrfDose=0;

	rbmDrfDose *= rbmR;
	bsDrfDose *= bsR;

	EvtMap->add(0, rbmDrfDose);
	EvtMap->add(1, bsDrfDose);

//	G4int copyNo = aStep->GetTrack()->GetOriginTouchable()->Get;
//	G4cout<<PhantomData->GetMaterialIndex(copyNo)<<"\t"<<index<<G4endl;
//	if(PhantomData->GetMaterialIndex(copyNo)!=index) {
//		EvtMap->add(-2, RBMdose);
//		EvtMap->add(-1, BSdose);
//	}

	return TRUE;
}

void BoneScorer::Initialize(G4HCofThisEvent* HCE)
{
	EvtMap = new G4THitsMap<G4double>(detector->GetName(), GetName());
    if(HCID < 0) HCID = GetCollectionID(0);
    HCE->AddHitsCollection(HCID,EvtMap);
}

void BoneScorer::EndOfEvent(G4HCofThisEvent*)
{;}

void BoneScorer::clear()
{
	EvtMap->clear();
}

std::pair<G4double, G4double> BoneScorer::GetDRFdose(G4double energy, G4double cellFlux, G4int organID)
{
	G4double rbmDRF, bsDRF;
	phantom->GetDRF(energy, organID, rbmDRF, bsDRF);
    return std::make_pair(rbmDRF*1e6*cellFlux*gray, bsDRF*1e6*cellFlux*gray);
}

