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
#include <iostream>
#include "RunAction.hh"

RunAction::RunAction(std::vector<ParallelPhantom*> _phantoms, G4Timer* _init)
:phantoms(_phantoms), fRun(0), numOfEvent(0), initTimer(_init), runTimer(0), frameID(0)
{
	// if(isMaster) runTimer = new G4Timer();
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(phantoms);
	return fRun;
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(G4int(numOfEvent*0.1));
	// if(runTimer) runTimer->Start();
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	// print the result only in the Master
	if(!isMaster) return;
	// runTimer->Stop();
	
	// get the run ID
	// runID = aRun->GetRunID();

	// Print the run result by G4cout and std::ofstream
	//
	// print by G4cout
	// PrintResult(G4cout);

	// print by std::ofstream
	// std::ofstream ofs("./PDClog/" + to_string(frameID) + ".out");
	std::ofstream ofs(to_string(frameID) + ".out");
	PrintResult(ofs);
	ofs.close();
}

void RunAction::PrintResult(std::ostream &out)
{
	// Print run result
	//
	// using namespace std;

	G4int maxNum = atoi(std::getenv("DCIR_MAX_NUM_OF_PEOPLE"));
	for(G4int i=0;i<maxNum;i++)
	{
		out<<"Radiologist #"<<i<<endl;
		ArrayXd dose = ((*fRun->GetDoseVal(i)).array() / phantoms[i]->doseMass)/(double)numOfEvent;
		ArrayXd dose2 = ((*fRun->GetDoseVal2(i)).array() / phantoms[i]->doseMass.square())/(double)numOfEvent;
		ArrayXd error = ((dose2 - dose*dose)/(double)numOfEvent).cwiseSqrt()/dose;
		for(auto iter:phantoms[i]->phantomData->doseNameMap)
			out<<iter.first<<"\t"<<iter.second<<"\t"<<dose(iter.first)<<"\t"<<error(iter.first)<<endl;
		dose = (*fRun->GetBoneDoseVal(i)).array()/(double)numOfEvent;
		dose2 = (*fRun->GetBoneDoseVal2(i)).array() / (double)numOfEvent;
		error = ((dose2 - dose*dose)/(double)numOfEvent).cwiseSqrt()/dose;
		for(int j=0;j<4;j++) out<<dose(j)<<"\t"<<error(j)<<endl;
	}
		// fRun->GetDoseVal()
	// EDEPMAP edepMap = *fRun->GetEdepMap();

	// out << G4endl
	//     << "=====================================================================" << G4endl
	//     << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	//     << "=====================================================================" << G4endl
	//     << "organ ID| "
	// 	<< setw(19) << "Organ Mass (g)"
	// 	<< setw(19) << "Dose (Gy/source)"
	// 	<< setw(19) << "Relative Error" << G4endl;

	// out.precision(3);
	// auto massMap = tetData->GetMassMap();
	// for(auto itr : massMap){
	// 	G4double meanDose    = edepMap[itr.first].first  / itr.second / numOfEvent;
	// 	G4double squareDoese = edepMap[itr.first].second / (itr.second*itr.second);
	// 	G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
	// 	G4double relativeE   = sqrt(variance)/meanDose;

	// 	// out << setw(8)  << itr.first << "| "
	// 	// 	<< setw(19) << fixed      << itr.second/g;
	// 	// out	<< setw(19) << scientific << meanDose/(joule/kg);
	// 	// out	<< setw(19) << fixed      << relativeE << G4endl;
	// 	out << setw(8)  << itr.first 
	// 		<< setw(19) << fixed      << itr.second/g;
	// 	out	<< setw(19) << scientific << meanDose/(joule/kg)*1e12; //pGy
	// 	out	<< setw(19) << fixed      << relativeE << G4endl;
	// }
	// out << "=====================================================================" << G4endl << G4endl;
}
