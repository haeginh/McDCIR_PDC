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
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	// pre init setting
	void SetTableRefPos(G4ThreeVector ref){ // right top corner (rel. to iso.)
		table_center.setX(ref.x()+table_half_size.x());
		table_center.setY(ref.y()-table_half_size.y());
		table_center.setZ(ref.z()-table_half_size.z());
	}

    // Operating Table
    void SetTablePose(G4ThreeVector _table_trans, G4double table_pivot_angle) 
	{  				
		table_rot->rotateZ(-table_pivot_angle);
		Eigen::Rotation2Dd rot(table_pivot_angle);
		Eigen::Translation2d trans0(table_rotation_center.x(), table_rotation_center.y());
		auto T = rot.inverse()*trans0*rot*trans0.inverse()*Eigen::Translation2d(_table_trans.x(), _table_trans.y());
		Eigen::Vector2d transT = (T*Eigen::Translation2d(table_center.x(), table_center.y())).translation();
		Eigen::Vector2d transC = (T*Eigen::Translation2d(curtain_center.x(), curtain_center.y())).translation();
		Eigen::Vector2d transP = (T*Eigen::Translation2d(phantom_center.x(), phantom_center.y())).translation();
		pv_table->SetTranslation(G4ThreeVector(transT(0), transT(1), table_center.z()+_table_trans.z()));
		pv_curtain->SetTranslation(G4ThreeVector(transC(0), transC(1), curtain_center.z()+_table_trans.z()));
		pv_phantom->SetTranslation(G4ThreeVector(transP(0), transP(1), phantom_center.z()+_table_trans.z()));


		// G4GeometryManager::GetInstance()->OpenGeometry(pv_frame);
		// table_trans = _table_trans;
		// G4ThreeVector frameOrigin = frame_default - carm_isocenter + table_trans; //before rotation
		// if(table_pivot_angle==0) 
		// {
		// 	pv_frame->SetTranslation(frameOrigin);
		// 	return;
		// }

		// frame_rotation_matrix->set(G4RotationMatrix::IDENTITY.axisAngle());
		// frame_rotation_matrix->rotateZ(-table_pivot_angle);
		// pv_frame->SetTranslation(frameOrigin + frame_rotation_matrix->inverse() * (table_rotation_center - frameOrigin));
		// G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_frame);
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
	// void ConstructCarmDet();
	// void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;

	// Operating Table
	G4VPhysicalVolume *pv_table, *pv_phantom, *pv_curtain;
	G4LogicalVolume* lv_curtain;
	G4ThreeVector table_rotation_center, table_half_size;//relative coordinate to isocenter
	G4ThreeVector table_center, curtain_center, phantom_center; // default: w/o table trans.
	G4RotationMatrix* table_rot; //inverse
	// G4RotationMatrix* frame_rotation_matrix; //inverse
	// G4ThreeVector table_trans; //frame_default, 
	G4String patient;

	//frame
	G4double head_margin, curtain_margin;

	//messenger
	DetectorMessenger* messenger;
	
};


#endif

