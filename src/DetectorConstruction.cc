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

#include "DetectorConstruction.hh"
#include "TETParameterisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Tet.hh"
#include "Eigen/Core"
#include "G4Trd.hh"
#include "G4EllipticalTube.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4MultiFunctionalDetector.hh"
#include "TETPSEnergyDeposit.hh"
#include "G4SDManager.hh"
#include "TETDRFScorer.hh"

using namespace std;

DetectorConstruction::DetectorConstruction(TETModelImport* _tetData, G4UIExecutive* _ui)
:lv_world(0) ,pv_world(0), table_rotation_center(-50.*mm, -1625*mm, 0.), head_margin(5.*cm), curtain_margin(10.*cm), tetData(_tetData), ui(_ui)
,isInstall(false)
{
	SetTableSize(G4ThreeVector(55*cm, 290*cm, 1.43*mm));
	SetTableRefPos(G4ThreeVector(-327, 398, -175)*mm);

	messenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
	delete messenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	cout << "Detector Construct()" << endl;
	SetupWorldGeometry();
	ConstructOperatingTable();
	ConstructXrayTube();
	ConstructFlatPanelDetector();
	ConstructPbCurtain();
	ConstructPbGlass();

	if (tetData == nullptr)
		return pv_world;

	if(!ui) {
		ConstructPatient(patient);
	}
	else {
		materialMap = tetData->GetMaterialMap();
		numTetMap   = tetData->GetNumTetMap();
		colourMap   = tetData->GetColourMap();
		VisualizePhantom();
	}
	return pv_world;
}

void DetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10 * 10 * 5 m3)
    G4double worldHalfX = 5. * m;
    G4double worldHalfY = 5. * m;
    G4double worldHalfZ = 2.5 * m;

	G4VSolid* worldSolid = new G4Box("worldSolid", worldHalfX, worldHalfY, worldHalfZ);
	lv_world = new G4LogicalVolume(worldSolid,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),"lv_world");
	pv_world = new G4PVPlacement(0,G4ThreeVector(), lv_world,"pv_world", 0, false,0,false);
	lv_world->SetVisAttributes(G4VisAttributes::GetInvisible());
	G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	va_World->SetForceWireframe(true);
	lv_world->SetVisAttributes(va_World);

	// // Room
	// G4VisAttributes* vis_room  = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.2));
	// vis_room->SetForceAuxEdgeVisible(true);

	// // Define floor
	// G4double zMargin = 0. * cm;
	// G4VSolid* sol_floor = new G4Box("sol_floor", worldHalfX, worldHalfY, 0.5*(worldHalfZ-zMargin));
	// G4LogicalVolume* lv_floor = new G4LogicalVolume(sol_floor, G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE"), "lv_floor");
	// lv_floor->SetVisAttributes(vis_room);
	// new G4PVPlacement(0, G4ThreeVector(0.,0.,-((G4Box*)sol_floor)->GetZHalfLength()-zMargin), lv_floor, "pv_floor", lv_world, 0, 0 );

	// // Ceiling
	// G4double wall_halfThickness = 0.25 * m;
	// G4VSolid* sol_ceiling = new G4Box("sol_ceiling", worldHalfX, worldHalfY, wall_halfThickness);
	// G4LogicalVolume* lv_ceiling = new G4LogicalVolume(sol_ceiling, G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE"), "lv_ceiling");
	// lv_ceiling->SetVisAttributes(vis_room);
	// new G4PVPlacement(0, G4ThreeVector(0.,0.,worldHalfZ-((G4Box*)sol_ceiling)->GetZHalfLength()), lv_ceiling, "pv_ceiling", lv_world, 0, 0 );

	// // Wall
	// G4VSolid* sol_wall = new G4Box("sol_wall", wall_halfThickness, (worldHalfY-wall_halfThickness), worldHalfZ-((G4Box*)sol_floor)->GetZHalfLength()-((G4Box*)sol_ceiling)->GetZHalfLength());
	// G4LogicalVolume* lv_wall = new G4LogicalVolume(sol_wall, G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE"), "lv_wall");
	// lv_wall->SetVisAttributes(vis_room);
	// new G4PVPlacement(0, G4ThreeVector(-worldHalfX+wall_halfThickness, wall_halfThickness,((G4Box*)sol_wall)->GetZHalfLength()-zMargin), lv_wall, "pv_wall1", lv_world, 0, 0);
	// new G4PVPlacement(0, G4ThreeVector(+worldHalfX-wall_halfThickness, -wall_halfThickness,((G4Box*)sol_wall)->GetZHalfLength()-zMargin), lv_wall, "pv_wall2", lv_world, 0, 0);
	// new G4PVPlacement(new G4RotationMatrix(0,0,90*deg), G4ThreeVector(wall_halfThickness, +worldHalfY-wall_halfThickness,((G4Box*)sol_wall)->GetZHalfLength()-zMargin), lv_wall, "pv_wall3", lv_world, 0, 0);
	// new G4PVPlacement(new G4RotationMatrix(0,0,90*deg), G4ThreeVector(-wall_halfThickness, -worldHalfY+wall_halfThickness,((G4Box*)sol_wall)->GetZHalfLength()-zMargin), lv_wall, "pv_wall4", lv_world, 0, 0);
}


#include "G4GeometryManager.hh"
void DetectorConstruction::ConstructOperatingTable()
{
	// Operating table
	G4Box* table = new G4Box("sol_table", table_half_size.x(), table_half_size.y(), table_half_size.z());
	G4LogicalVolume* lv_table = new G4LogicalVolume(table, G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"), "lv_table");
	lv_table->SetVisAttributes( new G4VisAttributes(G4Colour(0.2,0.2,0.2)) );
	G4Box* mattress = new G4Box("sol_mattress", table_half_size.x(), table_half_size.y(), 7.*cm * 0.5);

	G4Material* viscoelastic_foam = new G4Material("viscoelastic_foam", 57*kg/m3, G4NistManager::Instance()->FindOrBuildMaterial("G4_C"), kStateSolid, NTP_Temperature, STP_Pressure);
	G4LogicalVolume* lv_mattress = new G4LogicalVolume(mattress, viscoelastic_foam, "lv_mattress");
	
	lv_mattress->SetVisAttributes( new G4VisAttributes(G4Colour(1.,1.,0.)) );
	
	
	pv_table = new G4PVPlacement(0, G4ThreeVector(0, -851.72236, -251.1489
									-(G4double)(table)->GetZHalfLength()
									-(G4double)(mattress)->GetZHalfLength()), 	
									lv_table, "pv_table", lv_world, false, 0);
	new G4PVPlacement(0, G4ThreeVector(0, -851.72236, -251.1489), lv_mattress, "pv_table", lv_world, false, 0);


	// table_rot = new G4RotationMatrix;

	// // Operating table
	// G4Box* table = new G4Box("sol_table", table_half_size.x(), table_half_size.y(), table_half_size.z());
	// G4LogicalVolume* lv_table = new G4LogicalVolume(table, G4NistManager::Instance()->FindOrBuildMaterial("G4_Al"), "lv_table");
	// pv_table = new G4PVPlacement(table_rot, table_center, lv_table, "pv_table", lv_world, false, 1);
	// lv_table->SetVisAttributes( new G4VisAttributes(G4Colour(1.,1.,0.)));

	// // phantom box
	// G4LogicalVolume* lv_phantomBox = ConstructPatient(patient);
	// lv_phantomBox->SetVisAttributes(G4VisAttributes::GetInvisible());
	// phantom_center = table_center + 
	//                         G4ThreeVector(0., table_half_size.y()-((G4Box*) lv_phantomBox->GetSolid())->GetYHalfLength(),
	// 						                  table_half_size.z() + ((G4Box*) lv_phantomBox->GetSolid())->GetZHalfLength());
	// pv_phantom = new G4PVPlacement(table_rot, phantom_center,lv_phantomBox, "phantom box", lv_world, false, 0);

	// // Pb Curtain - lead equivalence 0.5 mm Pb
	// G4ThreeVector curtainHalfSize(0.25*mm, 250*mm, 200*mm);
	// G4Box* curtain = new G4Box("sol_curtain", curtainHalfSize.x(), curtainHalfSize.y(), curtainHalfSize.z());
	// curtain_center.setX(table_center.x() - table_half_size.x() - curtainHalfSize.x());
	// curtain_center.setY(table_center.y() + table_half_size.y() - curtain_margin - curtainHalfSize.y());
	// curtain_center.setZ(table_center.z() -table_half_size.z() - curtainHalfSize.z());
	// lv_curtain = new G4LogicalVolume(curtain, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_curtain");
	// pv_curtain = new G4PVPlacement(table_rot, curtain_center, lv_curtain, "pv_curtain", lv_world, false, 1);
	// lv_curtain->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.,0.8)) );
}

G4LogicalVolume* DetectorConstruction::ConstructPatient(G4String _patient)
{
	G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
	G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin()) * 0.5 + G4ThreeVector(1.,1.,1.) * cm;
	G4VSolid* sol_patientBox = new G4Box("sol_patientBox", halfSize.x(), halfSize.y(), halfSize.z());
	G4LogicalVolume* lv_patientBox = new G4LogicalVolume(sol_patientBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lv_patientBox");
	G4VisAttributes* va_patientBox = new G4VisAttributes();
	va_patientBox->SetForceWireframe(true);
	lv_patientBox->SetVisAttributes(va_patientBox);

	pv_patient = new G4PVPlacement(0,center,lv_patientBox,"lv_patientBox",lv_world,false,0);

	G4VSolid* sol_tet = new G4Tet("sol_tet", 
								  G4ThreeVector(0),
								  G4ThreeVector(1,0,0),
								  G4ThreeVector(0,1,0),
								  G4ThreeVector(0,0,1));
	lv_pTet = new G4LogicalVolume(sol_tet, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "lv_pTet");
	lv_pTet->SetVisAttributes(G4VisAttributes::GetInvisible());
	new G4PVParameterised("pv_patient", lv_pTet, lv_patientBox, kUndefined, tetData->GetNumTetrahedron(), new TETParameterisation(tetData));

	return lv_patientBox;
}


void DetectorConstruction::ConstructPbGlass()
{
	// Pb Glass - lead equivalence 0.5 mm Pb
	//G4Box* sol_glass = new G4Box("sol_glass", 60*0.5*cm, 0.5*0.5*mm, 75*0.5*cm);
	G4ThreeVector glassHalfSize((60)*0.5*cm, (75)*0.5*cm, 0.5*0.5*mm);
	G4Box* sol_glass = new G4Box("sol_glass", glassHalfSize.x(), glassHalfSize.y(), glassHalfSize.z());
	G4EllipticalTube* sol_glass_hole = new G4EllipticalTube("sol_glass_hole", 35*cm, 25*cm, 0.5*mm);
	G4VSolid* subtract = new G4SubtractionSolid("sub_glass", sol_glass, sol_glass_hole, 0, G4ThreeVector(glassHalfSize.x(), glassHalfSize.y(), 0));

	G4LogicalVolume* lv_glass = new G4LogicalVolume(subtract, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_glass");
	lv_glass->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.,0.5)) );

	Quaterniond q(0.6876994, 0.724684, -0.0300256, -0.0316404);
	Vector3d xAxis = q.matrix()*Vector3d::UnitX(); xAxis.normalized();
	AngleAxisd aa(M_PI, xAxis);
	Matrix3d m = aa.matrix() * q.matrix();
	G4ThreeVector xPrime(m(0), m(3), m(6));
	G4ThreeVector yPrime(m(1), m(4), m(7));
	G4ThreeVector zPrime(m(2), m(5), m(8));
	G4RotationMatrix* rot = new G4RotationMatrix();
	rot->rotateAxes(xPrime, yPrime, zPrime);

	pv_glass = new G4PVPlacement(rot, G4ThreeVector(-250.31587, -300, 220), lv_glass, "pv_glass", lv_world, false, 0);
}

void DetectorConstruction::ConstructPbCurtain()
{
	// Pb Curtain - lead equivalence 0.5 mm Pb
	G4Box* sol_curtain = new G4Box("sol_curtain", 80*0.5*cm, 0.5*0.5*mm, (62+5)*0.5*cm);
	lv_curtain = new G4LogicalVolume(sol_curtain, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_curtain");
	lv_curtain->SetVisAttributes( new G4VisAttributes(G4Colour(0.,1.,0.,0.5)) );
	
	Quaterniond q(0.6779879, 0.6779879, -0.2008292, -0.2008292);
	Vector3d xAxis = q.matrix()*Vector3d::UnitX(); xAxis.normalized();
	AngleAxisd aa(0.5 * M_PI, xAxis);
	Matrix3d m = aa.matrix() * q.matrix();
	G4ThreeVector xPrime(m(0), m(3), m(6));
	G4ThreeVector yPrime(m(1), m(4), m(7));
	G4ThreeVector zPrime(m(2), m(5), m(8));
	G4RotationMatrix* rot = new G4RotationMatrix();
	rot->rotateAxes(xPrime, yPrime, zPrime);
	new G4PVPlacement(rot, G4ThreeVector(-550.04186, -360.5048, -670.87331), lv_curtain, "pv_curtain", lv_world, false, 0);
}

void DetectorConstruction::ConstructXrayTube()
{
	G4VSolid* sol_xrayTube1 = new G4Trd("sol_xrayTube1", 20*cm, 10*cm, 20*cm, 10*cm, 10*cm);
	G4VSolid* sol_xrayTube1_hole = new G4Trd("sol_xrayTube2", 20*cm-2*mm, 10*cm-2*mm, 20*cm-2*mm, 10*cm-2*mm, 10*cm);
	G4EllipticalTube* sol_xrayTube2 = new G4EllipticalTube("sol_xrayTube2", 12*cm, 40*cm, 20*cm);
	G4EllipticalTube* sol_xrayTube2_hole = new G4EllipticalTube("sol_xrayTube2_hole", 12*cm-2*mm, 40*cm-2*mm, 20*cm-2*mm);

	G4RotationMatrix* rot1 = new G4RotationMatrix();
	rot1->rotateY(90*deg);
	G4SubtractionSolid* sub_xrayTube1 = new G4SubtractionSolid("sub_xrayTube1", sol_xrayTube1, sol_xrayTube1_hole, 0, G4ThreeVector(0, 0, 0));
	G4SubtractionSolid* sub_xrayTube1_1 = new G4SubtractionSolid("sub_xrayTube1_1", sub_xrayTube1, sol_xrayTube2, rot1, G4ThreeVector(0, 0, -20*cm));

	G4RotationMatrix* rot2 = new G4RotationMatrix();
	rot2->rotateY(-90*deg);
	G4SubtractionSolid* sub_xrayTube2 = new G4SubtractionSolid("sub_xrayTube2", sol_xrayTube2, sol_xrayTube2_hole, 0, G4ThreeVector(0, 0, 0));
	G4SubtractionSolid* sub_xrayTube2_2 = new G4SubtractionSolid("sub_xrayTube2_2", sub_xrayTube2, sol_xrayTube1, rot1, G4ThreeVector(-20*cm, 0, 0));

	G4UnionSolid* union_xrayTube = new G4UnionSolid("union_xrayTube", sub_xrayTube1_1, sub_xrayTube2_2, rot2, G4ThreeVector(0,0,-20*cm));
	G4LogicalVolume* lv_xrayTube = new G4LogicalVolume(union_xrayTube,G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"),"lv_xrayTube");
	lv_xrayTube->SetVisAttributes(new G4VisAttributes(G4Colour(0, 1, 1, 0.5)));

	G4RotationMatrix* rot3 = new G4RotationMatrix();
	rot3->rotateZ(90*deg);
	pv_xrayTube = new G4PVPlacement(rot3, G4ThreeVector(0,0,-810*mm)+G4ThreeVector(0,0,18*cm), lv_xrayTube, "pv_xrayTube", lv_world, false, 0, false);
}

void DetectorConstruction::ConstructFlatPanelDetector()
{
	// C-arm Flat Panel Detector, SID=100 cm
	G4LogicalVolume* lv_FPD = new G4LogicalVolume(new G4Box("pv_fpd", 52*cm*0.5, 42*cm*0.5, 1*cm), G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_epd");
	lv_FPD->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0,0.5)));
	pv_FPD = new G4PVPlacement(0, G4ThreeVector(0,0,19*cm), lv_FPD, "pv_FPD", lv_world, false, 0);
}


void DetectorConstruction::VisualizePhantom()
{
	if (isInstall) {
		for (auto itr: materialMap) {
			G4int idx = itr.first;
			delete lvTessMap[idx];
			delete pvTessMap[idx];
		}
		delete pv_patient;
		lvTessMap.clear();
		pvTessMap.clear();
	}
	
	G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
	G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin()) * 0.5 + G4ThreeVector(1.,1.,1.) * cm;
	G4VSolid* sol_patientBox = new G4Box("sol_patientBox", halfSize.x(), halfSize.y(), halfSize.z());
	G4LogicalVolume* lv_patientBox = new G4LogicalVolume(sol_patientBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lv_patientBox");
	G4VisAttributes* va_patientBox = new G4VisAttributes();
	va_patientBox->SetForceWireframe(true);
	lv_patientBox->SetVisAttributes(va_patientBox);
	pv_patient = new G4PVPlacement(0,center,lv_patientBox,"lv_patientBox",lv_world,false,0);

	MatrixXd U = tetData->GetAnimator()->GetU();
	
	for (auto itr: materialMap) {
		G4int idx = itr.first;
		G4String organName = itr.second->GetName();
		if (numTetMap[idx] == 0) continue;
		else if (organName.find("spongiosa")   != string::npos) continue;
		else if (organName.find("medullary")   != string::npos) continue;
		else if (organName.find("Muscle")      != string::npos) continue;
		else if (organName.find("RST")         != string::npos) continue;
		else if (organName.find("arteries")    != string::npos) continue;
		else if (organName.find("veins")       != string::npos) continue;
		// else if (organName.find("insensitive") != string::npos) continue;
		else if (organName.find("(") != string::npos && organName.find("surface") == string::npos) continue;
		else {
			cout << itr.first << " " << itr.second->GetName() << " " << numTetMap[idx] << flush;
			vector<tuple<int,int,int>> facePick = tetData->GetOrganSurfaceVec(idx);
			G4TessellatedSolid* tess = new G4TessellatedSolid;
			for (size_t i=0; i<facePick.size(); i++) {
			tess->AddFacet(new G4TriangularFacet(RowToG4Vec(U.row(get<0>(facePick[i])))-center, 
											 	 RowToG4Vec(U.row(get<1>(facePick[i])))-center,
											 	 RowToG4Vec(U.row(get<2>(facePick[i])))-center, ABSOLUTE) );
			} tess->SetSolidClosed(true);

			lvTessMap[idx] = new G4LogicalVolume(tess, itr.second, organName);
			lvTessMap[idx]->SetVisAttributes(new G4VisAttributes(colourMap[idx]));
			pvTessMap[idx] = new G4PVPlacement(0, G4ThreeVector(), lvTessMap[idx], to_string(idx)+"_"+organName, lv_patientBox, false, 0);
			cout << "...done" << endl;
		}
	}
	isInstall = true;
}

// void DetectorConstruction::ConstructSDandField()
// {
// 	if(tetData == nullptr) return;
// 	// Define detector and scorer
// 	G4SDManager* SDman = G4SDManager::GetSDMpointer();

// 	// Detector
// 	G4MultiFunctionalDetector* mfd_patient = new G4MultiFunctionalDetector("MFD_patient");
// 	SDman->AddNewDetector( mfd_patient );

// 	// Scorer
// 	mfd_patient->RegisterPrimitive(new TETPSEnergyDeposit("eDep", tetData));
// 	mfd_patient->RegisterPrimitive(new TETDRFScorer("DRF", tetData));
	
// 	// Attach the detector to logical volume
// 	if(!ui) SetSensitiveDetector(lv_pTet, mfd_patient);
// }

