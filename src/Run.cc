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

#include "Run.hh"
#include "RunAction.hh"
#include "G4RunManager.hh"
#include "G4VUserDetectorConstruction.hh"

Run::Run(RunAction* runAction)
:G4Run(), fCollID(-1), fCollID_drf(-1)
{
	parallel = static_cast<ParallelPhantom*>(
		G4RunManager::GetRunManager()->GetUserDetectorConstruction()->GetParallelWorld(0));
	G4int num = parallel->GetNumberOfPhantoms();
	edep = Eigen::MatrixXd::Zero(num, 33);
	edep2 = Eigen::MatrixXd::Zero(num, 33);
	currentEdep = Eigen::MatrixXd::Zero(num, 33);
	// eff = Eigen::VectorXd::Zero(num);
	// eff2 = Eigen::VectorXd::Zero(num);
	doseMap = runAction->GetDoseMap();
	effDoseVec = runAction->GetEffDoseVec();
}

Run::~Run()
{
}

void Run::RecordEvent(const G4Event* event)
{
	if(fCollID<0)
	{
		fCollID	= G4SDManager::GetSDMpointer()->GetCollectionID("phantom/edep");
		fCollID_drf	= G4SDManager::GetSDMpointer()->GetCollectionID("phantom/drf");
	}
	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));

	// sum up the energy deposition and the square of it
	currentEdep.setZero();
	for (auto itr : *evtMap->GetMap()) {
		auto idx = parallel->GetIdx(itr.first);
		for(auto d:doseMap[idx.second])
			currentEdep(d, idx.first) += *itr.second;
	}

	G4THitsMap<G4double>* evtMap_drf =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_drf));
	for (auto itr : *evtMap_drf)
	{
		G4int phantomID = std::floor(itr.first*0.1);
		G4int idx(0);
		if(itr.first-phantomID*10>0) idx = 10;
		currentEdep(idx, phantomID) += *itr.second;
	}

	currentEdep = currentEdep.array() * massInv.array();
	VectorXd currentEff = currentEdep * effDoseVec;
	currentEdep.rightCols(1) = currentEff;
	edep += currentEdep;
	edep2 += currentEdep.cwiseAbs2();
}

void Run::Merge(const G4Run* run)
{
	// merge the data from each thread
	auto localRun = static_cast<const Run*>(run);

	edep += localRun->edep;
	edep2 += localRun->edep2;
	G4Run::Merge(run);
}

