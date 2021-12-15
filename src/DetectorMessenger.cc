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
// DetectorMessenger.cc
// \author Haegin Han
//

#include "G4UIdirectory.hh"
#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWithABool.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include "RunAction.hh"
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction *_det)
	: G4UImessenger(), fDet(_det), tablePivot(0)
{
	fMachineDir = new G4UIdirectory("/machine/");
	
	fIsoCenterCmd = new G4UIcmdWith3VectorAndUnit("/machine/isocenter", this);
	fIsoCenterCmd->AvailableForStates(G4State_Idle, G4State_PreInit);
	fIsoCenterCmd->SetParameterName("isoX", "isoY", "isoZ", false);
	
	fTableRefCmd = new G4UIcmdWith3VectorAndUnit("/machine/tableRef", this);
	fTableRefCmd->AvailableForStates(G4State_PreInit);
	
	fTableTransCmd = new G4UIcmdWith3VectorAndUnit("/machine/tableTrans", this);
	fTablePivotCmd = new G4UIcmdWithADoubleAndUnit("/mahcine/pivot", this);
	fDetCmd = new G4UIcmdWith3Vector("/machine/c-arm", this);
	fRemoveDetCmd = new G4UIcmdWithoutParameter("/machine/removeDet", this);
	fGlassTransCmd = new G4UIcmdWith3VectorAndUnit("/machine/glassTrans", this);
	fGlassRotCmd = new G4UIcmdWith3Vector("/machine/glassRot", this);
	fRemoveGlassCmd = new G4UIcmdWithoutParameter("/machine/removeGlass", this);
	fCurtainCmd=  new G4UIcmdWithABool("/machine/useCurtain", this);

	fIsoCenterCmd->SetDefaultUnit("cm");
	fTableRefCmd->SetDefaultUnit("cm");
	fTableTransCmd->SetDefaultUnit("cm");
	fTablePivotCmd->SetDefaultUnit("deg");
	fGlassTransCmd->SetDefaultUnit("cm");

	fDetCmd->SetParameterName("lao[deg]", "caud[deg]", "sid[cm]", true, true);
	fGlassRotCmd->SetParameterName("axisX*[deg]", "axisY*[deg]", "axisZ*[deg]", false);
}

DetectorMessenger::~DetectorMessenger()
{
	delete fMachineDir;
	delete fIsoCenterCmd;
	delete fTableRefCmd;
	delete fTableTransCmd; //trans
	delete fTablePivotCmd; //pivot
	delete fDetCmd;		   //primary, secondary, sid
	delete fRemoveDetCmd;
	delete fGlassTransCmd;
	delete fGlassRotCmd; // axis * angle(in deg)
	delete fRemoveGlassCmd;
	delete fCurtainCmd;
}

#include "G4UImanager.hh"
void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if (command == fIsoCenterCmd)
	{
		G4ThreeVector iso = fIsoCenterCmd->GetNew3VectorValue(newValue);
		fDet->SetIsoCenter(iso);
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->SetIsoCenter(iso);
	}
	else if (command == fTableRefCmd)
	{
		fDet->SetTableRefPos(fTableRefCmd->GetNew3VectorValue(newValue));
	}
	else if (command == fTableTransCmd)
	{
		tableTrans = fTableTransCmd->GetNew3VectorValue(newValue);
		fDet->SetTablePose(tableTrans, tablePivot);
	}
	else if (command == fTablePivotCmd)
	{
		tablePivot = fTablePivotCmd->GetNewDoubleValue(newValue);
		fDet->SetTablePose(tableTrans, tablePivot);
	}
	else if (command == fDetCmd)
	{
		G4ThreeVector det = fDetCmd->GetNew3VectorValue(newValue);
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->SetCarmDetPose(det.x() * deg, det.y() * deg, det.z() * cm);
	}
	else if (command == fRemoveDetCmd)
	{
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->RemoveCarmDet();
	}
	else if (command == fGlassTransCmd)
	{
		glassTrans = fGlassTransCmd->GetNew3VectorValue(newValue);
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->SetGlassPose(glassTrans, glassAxis, glassTheta);
	}
	else if (command == fGlassRotCmd)
	{
		G4ThreeVector rot = fGlassRotCmd->GetNew3VectorValue(newValue);
		glassTheta = rot.mag() * deg;
		glassAxis = rot.unit();
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->SetGlassPose(glassTrans, glassAxis, glassTheta);
	}
	else if (command == fRemoveGlassCmd)
	{
		// ((ParallelGlass*) fDet->GetParallelWorld(0))->RemoveGlass();
	}
	else if (command == fCurtainCmd)
	{
		fDet->UseCurtain(fCurtainCmd->GetNewBoolValue(newValue));
	}
	// G4UImanager::GetUIpointer()->ApplyCommand("/vis/geometry/restore");
}
