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

Run::Run(std::vector<ParallelPhantom*> _phantoms)
:G4Run(), phantoms(_phantoms)
{
	maxNum = 5;
	if(std::getenv("DCIR_MAX_NUM_OF_PEOPLE")) maxNum = std::atoi(std::getenv("DCIR_MAX_NUM_OF_PEOPLE"));
	for(G4int i=0;i<maxNum;i++)
	{
		collIDs.push_back(G4SDManager::GetSDMpointer()->GetCollectionID("phantom"+std::to_string(i)+"_MFD/edep"));
		collIDs_bone.push_back(G4SDManager::GetSDMpointer()->GetCollectionID("phantom"+std::to_string(i)+"_MFD/bone"));
		doseVals.push_back(Eigen::VectorXd::Zero(phantoms[i]->phantomData->organ2dose.rows()));
		doseVals2.push_back(Eigen::VectorXd::Zero(phantoms[i]->phantomData->organ2dose.rows()));
		boneDose.push_back(Eigen::ArrayXd::Zero(4));
		boneDose2.push_back(Eigen::ArrayXd::Zero(4));
		skinDoses.push_back(Eigen::VectorXd::Zero(phantoms[i]->GetNumOfSkinFacet()));
		skinDoses2.push_back(Eigen::VectorXd::Zero(phantoms[i]->GetNumOfSkinFacet()));
		effVal.push_back(0); effVal2.push_back(0);
		boneVec.push_back(phantoms[i]->phantomData->GetBoneIdxVec());
	}
}

Run::~Run()
{
	collIDs.clear();
	// organDose.clear();
	// boneDose.clear();
}

void Run::RecordEvent(const G4Event* event)
{
	// Hits collections
	//
	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
	if(!HCE) return;

	for(int i=0;i<maxNum;i++)
	{
		G4THitsMap<G4double>* evtMap =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(collIDs[i]));
		G4THitsMap<G4double>* evtMap_bone =
			static_cast<G4THitsMap<G4double>*>(HCE->GetHC(collIDs_bone[i]));
		Eigen::VectorXd edep = Eigen::VectorXd::Zero(phantoms[i]->phantomData->organ2dose.cols());
		Eigen::ArrayXd edep_bone = Eigen::VectorXd::Zero(4);
	
		for (auto itr : *evtMap->GetMap())
		{	
			if(itr.first>0) edep(itr.first) = *itr.second;
			else
			{
				skinDoses[i](itr.first) += *itr.second;
				skinDoses2[i](itr.first) += (*itr.second)*(*itr.second);
			}
		}
		for (auto itr : *evtMap_bone->GetMap()) 
			edep_bone(itr.first) = *itr.second;

		for(G4int b:boneVec[i]){
			G4double rbmR, bsR, massInv(phantoms[i]->GetBoneMassInv(b));
			phantoms[i]->GetBoneMassRatio(b, rbmR, bsR);
			edep_bone(2) += edep(b) * rbmR * massInv;
			edep_bone(3) += edep(b) * bsR * massInv;
		}
		
		edep = phantoms[i]->phantomData->organ2dose * edep;
		edep(0) = edep_bone(0);
		edep(1) = edep_bone(1);

		doseVals[i] += edep;
		doseVals2[i] += (edep.array()*edep.array()).matrix();

		boneDose[i] += edep_bone;
		boneDose2[i] += edep_bone*edep_bone;
		
		G4double eff = phantoms[i]->phantomData->effCal * edep;
		effVal[i] += eff;
		effVal2[i] += eff*eff;
	}
}

void Run::Merge(const G4Run* run)
{
	// merge the data from each thread
	auto localRun = static_cast<const Run*>(run);

	for(int i=0;i<maxNum;i++){
		doseVals[i]+=localRun->doseVals[i];
		doseVals2[i]+=localRun->doseVals2[i];
		effVal[i] += localRun->effVal[i];
		effVal2[i] += localRun->effVal2[i];
		boneDose[i] += localRun->boneDose[i];
		boneDose2[i] += localRun->boneDose2[i];
	}

	G4Run::Merge(run);
}

