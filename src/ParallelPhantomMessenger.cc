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
// ParallelMessenger.cc
// \author Haegin Han
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithAString.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "RunAction.hh"
#include "ParallelPhantomMessenger.hh"
#include "ParallelPhantom.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

ParallelPhantomMessenger::ParallelPhantomMessenger(ParallelPhantom* _phantom)
:G4UImessenger(), fPhantom(_phantom)
{
	fPhantomDir = new G4UIdirectory("/phantom/");
	fInitCmd = new G4UIcmdWithAString("/phantom/init", this);
	fDeformCmd = new G4UIcmdWithAString("/phantom/deform", this);
	//(id) (x,y,z,w) (x,y,z,w) .... (x,y,z)_root
	fSetFrameCmd = new G4UIcmdWithAnInteger("/phantom/setFrame", this);
}

ParallelPhantomMessenger::~ParallelPhantomMessenger() {
	delete fPhantomDir;
	delete fInitCmd;
	delete fDeformCmd; 
	delete fSetFrameCmd;
}

void ParallelPhantomMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
	if(command == fInitCmd){
		fPhantom->InitializePhantoms(newValue);
	}
	else if(command == fDeformCmd){
		std::stringstream ss(newValue);
		G4int id;
		ss>>id;
		G4double dat;
		std::vector<G4double> datVec;
		while(ss>>dat)
			datVec.push_back(dat);
		RotationList vQ;
		Vector3d root(datVec[datVec.size()-3], datVec[datVec.size()-2], datVec[datVec.size()-1]);
		for(int i=0;i<datVec.size()-3;i+=4)
			vQ.push_back(Quaterniond(datVec[i],datVec[i+1],datVec[i+2],datVec[i+3]));
		fPhantom->Deform(id, vQ, root*cm);
	}
	else if(command == fSetFrameCmd){
		fPhantom->SetNewFrame(fSetFrameCmd->GetNewIntValue(newValue));
	}
}

// void ParallelPhantomMessenger::ReadPostureData(G4String fileName)
// {
// 	ifstream ifs(fileName);
// 	if(!ifs.is_open())
// 	{
// 		cout<<"fileName is not open"<<endl;
// 	}
// 	vQ_vec.clear();
// 	roots.clear();

// 	G4String line;
// 	while(getline(ifs, line))
// 	{
// 		if(line.empty()) continue;
// 		stringstream ss(line);
// 		G4double x, y, z, w;
// 		ss>>x>>y>>z;
// 		roots.push_back(Vector3d(x,y,z)*cm);
// 		RotationList vQ;
// 		for(int i=0;i<22;i++)
// 		{
// 			ss>>w>>x>>y>>z;
// 			vQ.push_back(Quaterniond(w, x, y, z));
// 		}
// 		vQ_vec.push_back(vQ);
// 	}
// 	ifs.close();
// }
