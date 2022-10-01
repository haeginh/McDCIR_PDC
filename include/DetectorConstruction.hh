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
#include "G4Tet.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "DetectorMessenger.hh"
#include "G4GeometryManager.hh"
#include <Eigen/Geometry>

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    virtual ~DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();

	// pre init setting
	void SetTableRefPos(G4ThreeVector ref){table_ref_pos = ref;}

    // Operating Table
    void SetTablePose(G4ThreeVector _table_trans, G4double table_pivot_angle) 
	{  				
		// G4GeometryManager::GetInstance()->OpenGeometry(worldPhysical);
		
		auto table_rot = pv_table->GetRotation();
		table_rot->set(G4RotationMatrix::IDENTITY.axisAngle());
		table_rot->rotateZ(-table_pivot_angle); //inverse

		Eigen::Rotation2Dd rot(table_pivot_angle);
		Eigen::Translation2d trans0(table_rotation_center.x(), table_rotation_center.y());
		auto T = trans0*rot*trans0.inverse()*Eigen::Translation2d(_table_trans.x(), _table_trans.y());
		Eigen::Vector2d transT = (T*Eigen::Translation2d(table_center0.x(), table_center0.y())).translation();
		Eigen::Vector2d transC = (T*Eigen::Translation2d(curtain_center0.x(), curtain_center0.y())).translation();
		Eigen::Vector2d transP = (T*Eigen::Translation2d(phantom_center0.x(), phantom_center0.y())).translation();
		G4cout<<table_center0<<" "<<_table_trans<<" "<<transT<<G4endl;
		pv_table->SetTranslation(G4ThreeVector(transT(0), transT(1), table_center0.z()+_table_trans.z()));
		pv_curtain->SetTranslation(G4ThreeVector(transC(0), transC(1), curtain_center0.z()+_table_trans.z()));
		pv_phantom->SetTranslation(G4ThreeVector(transP(0), transP(1), phantom_center0.z()+_table_trans.z()));
	}

	void UseCurtain(G4bool use = true)
	{
		if(use) pv_curtain->GetLogicalVolume()->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"));
		else pv_curtain->GetLogicalVolume()->SetMaterial(G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"));
	}

	void SetPatientName(G4String _patient){patient = _patient;}

private:
	void SetTableSize(G4ThreeVector size){ // full size
		table_size = size;
	}

    void SetupWorldGeometry();

	void ConstructOperatingTable();
	G4LogicalVolume* ConstructPatient(G4String _patient);
	// void ConstructCarmDet();
	// void ConstructPbGlass();

	G4LogicalVolume*   worldLogical;
	G4VPhysicalVolume* worldPhysical;

	G4double isoZ;
	// Operating Table
	G4ThreeVector table_ref_pos; // right top corner
	G4ThreeVector table_rotation_center;
	G4ThreeVector table_size, curtain_size;
	G4ThreeVector table_center0, phantom_center0, curtain_center0;
	G4VPhysicalVolume *pv_table, *pv_phantom, *pv_curtain;

	G4String patient;

	//frame
	G4double head_margin, curtain_margin;

	//messenger
	DetectorMessenger* messenger;

	//phantom
	std::vector<G4Tet*> tetVec;
	Eigen::ArrayXi idArray;
	std::map<G4int, G4Material*> matMap;
public:
	G4Tet* GetTet(G4int id){return tetVec[id];}
	// G4int GetOrganId(G4int id){return idArray(id);}
	G4Material* GetMaterial(G4int id){return matMap[idArray(id)];}
	// G4Material* GetTissueMaterial(){return matMap[12200];}
};


#endif

