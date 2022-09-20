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
#include "G4UIcmdWithAnInteger.hh"
#include "G4RunManager.hh"
#include <sstream>
#include <vector>
#include <fstream>
#include "RunAction.hh"
#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4RunManager.hh"
#include "G4GeometryManager.hh"

DetectorMessenger::DetectorMessenger(DetectorConstruction *_det)
	: G4UImessenger(), fDet(_det), tablePivot(0)
{
	maxPeople = std::atoi(std::getenv("DCIR_MAX_NUM_OF_PEOPLE"));

	fMachineDir = new G4UIdirectory("/machine/");
	
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

	fTableRefCmd->SetDefaultUnit("cm");
	fTableTransCmd->SetDefaultUnit("cm");
	fTablePivotCmd->SetDefaultUnit("deg");
	fGlassTransCmd->SetDefaultUnit("cm");

	fDetCmd->SetParameterName("lao[deg]", "caud[deg]", "sid[cm]", true, true);
	fGlassRotCmd->SetParameterName("axisX*[deg]", "axisY*[deg]", "axisZ*[deg]", false);

	fRecordDir = new G4UIdirectory("/phantom/");
	fInstallPhantom = new G4UIcmdWithAnInteger("/phantom/install", this);
	fUninstallPhantom = new G4UIcmdWithAnInteger("/phantom/uninstall", this);

	fRecordDir = new G4UIdirectory("/record/");
	fReadRecordData = new G4UIcmdWithAString("/record/read", this);
	fSetFrameFromRecordData = new G4UIcmdWithAnInteger("/record/setFrame", this);
}

DetectorMessenger::~DetectorMessenger()
{
	delete fMachineDir;
	delete fTableRefCmd;
	delete fTableTransCmd; //trans
	delete fTablePivotCmd; //pivot
	delete fDetCmd;		   //primary, secondary, sid
	delete fRemoveDetCmd;
	delete fGlassTransCmd;
	delete fGlassRotCmd; // axis * angle(in deg)
	delete fRemoveGlassCmd;
	delete fCurtainCmd;

	delete fRecordDir;
	delete fReadRecordData;
	delete fSetFrameFromRecordData;
}

#include "G4UImanager.hh"
void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if (command == fTableRefCmd)
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
	else if (command == fInstallPhantom)
	{
		G4int parallelID = maxPeople-1-fSetFrameFromRecordData->GetNewIntValue(newValue);
		G4GeometryManager::GetInstance()->OpenGeometry();
		fDet->InstallPhantom(parallelID);
	}
	else if (command == fUninstallPhantom)
	{
		G4int parallelID = maxPeople-1-fSetFrameFromRecordData->GetNewIntValue(newValue);
		G4GeometryManager::GetInstance()->OpenGeometry();
		fDet->UninstallPhantom(parallelID);
		// G4GeometryManager::GetInstance()->CloseGeometry();
	}
	else if (command == fReadRecordData)
	{
		std::ifstream ifs(newValue);
		if(!ifs.is_open())
		{
			G4Exception("Messenger", "", JustWarning, G4String("There is no "+newValue).c_str());
			return;
		}
		recordData.clear();
		G4int dataSize = 18+maxPeople*145;
        std::vector<float> frameData(dataSize);
		while (ifs.read((char *)frameData.data(), dataSize * sizeof(G4float)))
            recordData.push_back(frameData);
		ifs.close();
		G4cout<<"Read "<<recordData.size()<<" frames"<<G4endl;
		return;
	}
	else if (command == fSetFrameFromRecordData)
	{
		G4int frameNo = fSetFrameFromRecordData->GetNewIntValue(newValue);
		if(frameNo>=recordData.size())
		{
			G4Exception("Messenger", "", JustWarning, G4String("Wrong frame number. Record data size is "+std::to_string(frameNo)).c_str());
			return;
		}
		DataSet data = ReadAFrameData(recordData[frameNo]);
		for(int i=0;i<maxPeople;i++) // i is for phanotm ID (parallel ID is reversed)
		{
			G4int parallelID = maxPeople-1-i;
			if(data.bodyMap.find(i)==data.bodyMap.end())
				fDet->UninstallPhantom(parallelID);
			else 
				fDet->InstallPhantom(parallelID, data.bodyMap[i].posture, data.bodyMap[i].jointC.row(0));
		}
	}
	G4RunManager::GetRunManager()->GeometryHasBeenModified(true);
	// G4UImanager::GetUIpointer()->ApplyCommand("/vis/geometry/restore");
}

DataSet DetectorMessenger::ReadAFrameData(std::vector<G4float> frameData)
{
    DataSet data;
    data.time = frameData[0];
    data.kVp = frameData[1];
    data.mA = frameData[2];
    data.dap = frameData[3];
    data.cArm(0) = frameData[4];
    data.cArm(1) = frameData[5];
    data.cArm(2) = frameData[6];
    data.bed(0) = frameData[7];
    data.bed(1) = frameData[8];
    data.bed(2) = frameData[9];
    data.glassChk = frameData[10];
    Affine3d glassAff = Affine3d::Identity();
    glassAff.rotate(Quaterniond(frameData[11], frameData[12], frameData[13], frameData[14]));
    glassAff.translate(Vector3d(frameData[15], frameData[16], frameData[17]));
    data.glass_aff = glassAff;
    int pos = 18;
    for (int i = 0; i < 5; i++)
    {
        int num = pos + i * 145;
        if (!frameData[num++])
            continue;
        for (int b = 0; b < data.bodyMap[i].posture.size(); b++)
        {
            data.bodyMap[i].posture[b].x() = frameData[num++];
            data.bodyMap[i].posture[b].y() = frameData[num++];
            data.bodyMap[i].posture[b].z() = frameData[num++];
            data.bodyMap[i].posture[b].w() = frameData[num++];
        }
        for (int c = 0; c < data.bodyMap[i].jointC.rows(); c++)
        {
            data.bodyMap[i].jointC(c, 0) = frameData[num++];
            data.bodyMap[i].jointC(c, 1) = frameData[num++];
            data.bodyMap[i].jointC(c, 2) = frameData[num++];
        }
    }
    return data;
}
