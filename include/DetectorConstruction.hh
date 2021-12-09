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
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

#include <cmath>

#include "globals.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "DetectorMessenger.hh"
#include "ParallelPhantom.hh"
#include "ParallelGlass.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	// pre init setting
	void SetIsoCenter(G4ThreeVector iso){
		carm_isocenter = iso;
		// pv_frame->SetTranslation(frame_default-carm_isocenter+table_trans);
		// ((ParallelGlass*)GetParallelWorld(0))->SetIsoCenter(iso);
		((ParallelPhantom*)GetParallelWorld(0))->SetIsoCenter(iso);
	}
	void SetTableRefPos(G4ThreeVector ref){table_ref_pos = ref;}

    // Operating Table
    void SetTablePose(G4ThreeVector _table_trans, G4double table_pivot_angle) 
	//_table_trans is the translation from when rotation center is at the center of the table
	{  				
		G4GeometryManager::GetInstance()->OpenGeometry(pv_frame);
		table_trans = _table_trans;
		G4ThreeVector frameOrigin = frame_default - carm_isocenter + table_trans; //before rotation
		if(table_pivot_angle==0) 
		{
			pv_frame->SetTranslation(frameOrigin);
			return;
		}

		frame_rotation_matrix->set(G4RotationMatrix::IDENTITY.axisAngle());
		frame_rotation_matrix->rotateZ(-table_pivot_angle);
		pv_frame->SetTranslation(frameOrigin + frame_rotation_matrix->inverse() * (table_rotation_center - frameOrigin));
		G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_frame);
	}

	void UseCurtain(G4bool use = true)
	{
		if(use) lv_curtain->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"));
		else lv_curtain->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
	}

	void SetPatientName(G4String _patient){patient = _patient;}

private:
    void SetupWorldGeometry();

	void ConstructOperatingTable();
	G4LogicalVolume* ConstructPatient(G4String _patient);
	// void ConstructCarmDet();
	// void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;
	G4ThreeVector carm_isocenter;

	// Operating Table
	G4VPhysicalVolume* pv_frame;
	G4LogicalVolume* lv_curtain;
	G4ThreeVector table_ref_pos; // right top corner
	G4ThreeVector table_rotation_center;//, table_default_trans; //program input (relative coordinate to ChArUco)
	G4RotationMatrix* frame_rotation_matrix; //inverse
	G4ThreeVector frame_default, table_trans;
	G4String patient;

	//frame
	G4double head_margin, curtain_margin;

	//messenger
	DetectorMessenger* messenger;
	
};


#endif

