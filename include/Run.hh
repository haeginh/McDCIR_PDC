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
#include "ParallelPhantom.hh"
#include <Eigen/Core>

typedef std::map<G4int ,std::map<G4int, std::pair<G4double, G4double>>> EDEPMAP;

class Run : public G4Run
{
public:
	Run(std::vector<ParallelPhantom*>);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

    Eigen::VectorXd* GetDoseVal(G4int idx) {return &(doseVals[idx]);}
    Eigen::VectorXd* GetSkinDoseVal(G4int idx) {return &(skinDoses[idx]);}
    Eigen::ArrayXd* GetBoneDoseVal(G4int idx) {return &(boneDose[idx]);}
    Eigen::VectorXd* GetDoseVal2(G4int idx) {return &(doseVals2[idx]);}
    Eigen::VectorXd* GetSkinDoseVal2(G4int idx) {return &(skinDoses2[idx]);}
    Eigen::ArrayXd* GetBoneDoseVal2(G4int idx) {return &(boneDose2[idx]);}
    G4double GetEffDose(G4int idx) {return effVal[idx];}
    G4double GetEffDose2(G4int idx) {return effVal2[idx];}

private:
	std::vector<ParallelPhantom*> phantoms;
	G4int maxNum;
	// EDEPMAP organDose, boneDose;
	std::vector<Eigen::VectorXd> doseVals, doseVals2, skinDoses, skinDoses2;
	std::vector<Eigen::ArrayXd> boneDose, boneDose2;
	std::vector<G4int> effVal, effVal2;
	std::vector<G4int> collIDs, collIDs_bone;
	std::vector<std::vector<G4int>> boneVec;
};

#endif
 