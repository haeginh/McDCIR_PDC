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
#include <Eigen/Core>

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction(TETModelImport* tetData, G4UIExecutive* ui);
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
	//virtual void ConstructSDandField();

	// pre init setting
	void SetTableRefPos(G4ThreeVector ref){ // right top corner (rel. to iso.)
		table_center.setX(ref.x()+table_half_size.x());
		table_center.setY(ref.y()-table_half_size.y());
		table_center.setZ(ref.z()-table_half_size.z());
	}

    // Operating Table
    void SetTablePose(G4ThreeVector _table_trans, G4double table_pivot_angle) 
	{  				
		G4GeometryManager::GetInstance()->OpenGeometry(pv_world);
		
		table_rot->set(G4RotationMatrix::IDENTITY.axisAngle());
		table_rot->rotateZ(-table_pivot_angle);

		Eigen::Rotation2Dd rot(table_pivot_angle);
		Eigen::Translation2d trans0(table_rotation_center.x(), table_rotation_center.y());
		auto T = trans0*rot*trans0.inverse()*Eigen::Translation2d(_table_trans.x(), _table_trans.y());
		Eigen::Vector2d transT = (T*Eigen::Translation2d(table_center.x(), table_center.y())).translation();
		Eigen::Vector2d transC = (T*Eigen::Translation2d(curtain_center.x(), curtain_center.y())).translation();
		Eigen::Vector2d transP = (T*Eigen::Translation2d(phantom_center.x(), phantom_center.y())).translation();
		pv_table->SetTranslation(G4ThreeVector(transT(0), transT(1), table_center.z()+_table_trans.z()));
		pv_curtain->SetTranslation(G4ThreeVector(transC(0), transC(1), curtain_center.z()+_table_trans.z()));
		pv_phantom->SetTranslation(G4ThreeVector(transP(0), transP(1), phantom_center.z()+_table_trans.z()));

		G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_world);
		G4RunManager::GetRunManager()->GeometryHasBeenModified();
	}

	void UseCurtain(G4bool use = true)
	{
		if(use) lv_curtain->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"));
		else lv_curtain->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
	}

	void SetPatientName(G4String _patient){patient = _patient;}

private:
	void SetTableSize(G4ThreeVector size){ // full size
		table_half_size = size*0.5;
	}

    void SetupWorldGeometry();
	void ConstructOperatingTable();
	G4LogicalVolume* ConstructPatient(G4String _patient);
	void ConstructFlatPanelDetector();
	void ConstructPbCurtain();
	void ConstructPbGlass();
	void ConstructXrayTube();
	void VisualizePhantom();

	G4ThreeVector RowToG4Vec(const RowVector3d& row){return G4ThreeVector(row(0), row(1), row(2));}

private:
	// phantom
	TETModelImport* tetData;
	//messenger
	DetectorMessenger* messenger;
	
	// World
	G4LogicalVolume*   lv_world;
	G4VPhysicalVolume* pv_world;

	// Operating Table
	G4VPhysicalVolume *pv_table, *pv_phantom, *pv_curtain;
	G4LogicalVolume* lv_curtain;
	G4ThreeVector table_rotation_center, table_half_size;//relative coordinate to isocenter
	G4ThreeVector table_center, curtain_center, phantom_center; // default: w/o table trans.
	G4RotationMatrix* table_rot; //inverse

	// C-arm
	G4VPhysicalVolume* pv_FPD;
	G4VPhysicalVolume* pv_xrayTube;

	// Patient
	G4String patient;
	G4LogicalVolume*   lv_patient;
	G4VPhysicalVolume* pv_patient;
	G4LogicalVolume*   lv_pTet;

	// Glass
	G4VPhysicalVolume* pv_glass;

	//frame
	G4double head_margin, curtain_margin;

	
	
	
	
	// Phantom Visualize
	G4UIExecutive* ui;
	G4bool isInstall;
	map<G4int, G4Material*> materialMap;
	map<G4int, G4int>       numTetMap;
	map<G4int, G4LogicalVolume*> lvTessMap;
	map<G4int, G4Colour>    colourMap;
	map<G4int, G4PVPlacement*> pvTessMap;
	
};


#endif

