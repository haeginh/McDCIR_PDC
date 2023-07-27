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

class RunAction;
class Run : public G4Run
{
public:
	Run(RunAction*);
	virtual ~Run();

	virtual void RecordEvent(const G4Event*);
    virtual void Merge(const G4Run*);

	void DoseRead(G4String);
    Eigen::MatrixXd GetEdepMat() {return edep;}
    Eigen::MatrixXd GetEdep2Mat() {return edep2;}
	void SetMassInv(const MatrixXd &mass) {
		massInv = (mass.array()==0).cast<double>().matrix() + mass;
		massInv = massInv.cwiseInverse();
		massInv.conservativeResize(massInv.rows(), massInv.cols()+1);
		massInv.rightCols(1) = Eigen::VectorXd::Ones(massInv.rows());
	}

private:
	Eigen::MatrixXd edep, edep2, currentEdep, massInv;
	std::map<G4int, std::vector<G4int>>            doseMap; //suborgan:dose
	Eigen::VectorXd                                effDoseVec;
	G4int fCollID, fCollID_drf;
	ParallelPhantom* parallel;
};

#endif
