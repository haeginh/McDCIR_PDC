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
// TETRunAction.cc
// \file   MRCP_GEANT4/External/src/TETRunAction.cc
// \author Haegin Han
//

#include "G4Timer.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include <iostream>
#include "RunAction.hh"
#include "ParallelPhantom.hh"
#include "DetectorConstruction.hh"

RunAction::RunAction(G4String _output, G4Timer* _init)
: fRun(0), numOfEvent(0), runID(0), outputFile(_output), initTimer(_init), runTimer(0)
{
	DoseRead("dose.txt");
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(this);
	return fRun;
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(G4int(numOfEvent*0.1));
	
	auto parallel = static_cast<ParallelPhantom*>(
		G4RunManager::GetRunManager()->GetUserDetectorConstruction()->GetParallelWorld(0));
	mass = Eigen::MatrixXd::Zero(parallel->GetNumberOfPhantoms(), 32);
	for(G4int i=0;i<mass.rows();i++)
	{
		auto massMap = parallel->GetMassMap(i);
		for(auto iter:massMap)
		{
			for(G4int d:doseMap[iter.first])
				mass(i, d) += iter.second;
		}
	}
	fRun->SetMassInv(mass);
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;

	// get the run ID
	runID = aRun->GetRunID();

	// Print the run result by G4cout and std::ofstream
	//
	// print by G4cout
	// PrintResult(G4cout);

	// print by std::ofstream
	std::ofstream ofs(outputFile + to_string(runID) + ".out");
	PrintResult(ofs);
	ofs.close();
}

void RunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	using namespace std;

	auto parallel = static_cast<ParallelPhantom*>(
		G4RunManager::GetRunManager()->GetUserDetectorConstruction()->GetParallelWorld(0));

	// out << G4endl
	//     << "=====================================================================" << G4endl
	//     << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	//     << "=====================================================================" << G4endl
	//     << "organ ID| "
	// 	<< setw(19) << "Organ Mass (g)"
	// 	<< setw(19) << "Dose (Gy/source)"
	// 	<< setw(19) << "Relative Error" << G4endl;

	out.precision(3);
	for(G4int i=0;i<32;i++) out<<organNameMap[i]<<"\t";
	out<<"eff"<<endl;
	Eigen::MatrixXd edep = fRun->GetEdepMat()/(double)numOfEvent;
	Eigen::MatrixXd err = ((fRun->GetEdep2Mat()/(double)numOfEvent) - edep.cwiseAbs2()).cwiseSqrt().array() / edep.array();
	out<<edep/gray<<endl<<endl<<err<<endl<<endl<<mass;
}

void RunAction::DoseRead(G4String file){
	std::ifstream ifs(file);
	if(!ifs.is_open())
		G4Exception("PhantomData::DataRead","",FatalErrorInArgument,
		G4String("      There is no " + file ).c_str());
	
	G4String line;
	effDoseVec = Eigen::VectorXd::Zero(32);
	while(getline(ifs, line))
	{
		if(line.empty()) continue;
		std::stringstream ss(line);
		G4int id, id1(-1);
		G4double eff, mass(0);
		G4String name;
		ss>>id>>name>>eff;
		if(name.empty()) continue;
		organNameMap[id] = name;
		effDoseVec[id] = eff;
		while(ss>>id1)
		{
			if(id1<0) continue;
			doseMap[id1].push_back(id);
			id1=-1;
		}
	}
}