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
// \author: Haegin Han
//
#include "ParallelPhantom.hh"
#include "TETParameterisation.hh"

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "TETParameterisation.hh"
#include "TETPSEnergyDeposit.hh"
#include "ParallelPhantomMessenger.hh"
#include "G4UImanager.hh"
#include "TETDRFScorer.hh"

ParallelPhantom
::ParallelPhantom(G4String parallelWorldName, TETModelImport* _tetData, G4UIExecutive* _ui)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false), tetData(_tetData), ui(_ui), isInstall(false)
{
  messenger = new ParallelPhantomMessenger(this);
}

ParallelPhantom::~ParallelPhantom()
{
  delete messenger;
}

void ParallelPhantom::Construct()
{
  cout << "Parallel Construct()" << endl;
  if(fConstructed) return;
  fConstructed = true;

  //
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  lv_world = ghostWorld->GetLogicalVolume();

  if (!ui) {
    G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
    G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin()) * 0.5 + G4ThreeVector(1.,1.,1.) * cm;
    G4VSolid* sol_doctorBox = new G4Box("sol_doctorBox", halfSize.x(), halfSize.y(), halfSize.z());
    G4LogicalVolume* lv_doctorBox = new G4LogicalVolume(sol_doctorBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lv_doctorBox");
    G4VisAttributes* va_doctorBox = new G4VisAttributes();
    va_doctorBox->SetForceWireframe(true);
    lv_doctorBox->SetVisAttributes(va_doctorBox);
    pv_doctor = new G4PVPlacement(0,G4ThreeVector(),lv_doctorBox,"lv_doctorBox",lv_world,false,0);
    G4VSolid* sol_tet = new G4Tet("sol_tet", 
                  G4ThreeVector(0),
                  G4ThreeVector(1,0,0),
                  G4ThreeVector(0,1,0),
                  G4ThreeVector(0,0,1));
    lv_dTet = new G4LogicalVolume(sol_tet, G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic"), "lv_cTet");
    new G4PVParameterised("pv_doctor", lv_dTet, lv_doctorBox, kUndefined, tetData->GetNumTetrahedron(), new TETParameterisation(tetData));
  }
  else {
    materialMap = tetData->GetMaterialMap();
    numTetMap   = tetData->GetNumTetMap();
    colourMap   = tetData->GetColourMap();
    VisualizePhantom();
  }
  // //
  // // parallel world placement box
  // //
  // G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
  // G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin())*0.5 + G4ThreeVector(5., 5., 5.)*cm; //5-cm-margin
  // G4VSolid* paraBox = new G4Box("phantomBox",halfSize.x(),halfSize.y(),halfSize.z());
  // G4LogicalVolume* lv_phantomBox = new G4LogicalVolume(paraBox,0,"phantomBox");
  // lv_phantomBox->SetVisAttributes(G4VisAttributes::GetInvisible());
  // pv_doctor = new G4PVPlacement(0,center,lv_phantomBox,"phantomBox",lv_world,false,0);

  // //
  // // mother of parallel world parameterized volumes
  // //

  // // G4Material* tissue = G4Material::GetMaterial("G4_TISSUE_SOFT_ICRP");
  // // G4cout<<tetData->GetNumTetrahedron()<<G4endl;
  // // for(size_t i=0;i<tetData->GetNumTetrahedron();i++)
  // // {
  // //   new G4PVPlacement(0,G4ThreeVector(),new G4LogicalVolume(tetData->GetTetrahedron(i),tissue, "tet"),"tet",worldLogical,false,0);
  // // }
  // lv_dTet = new G4LogicalVolume(new G4Tet("tet", G4ThreeVector(), G4ThreeVector(0, 0, 1*cm),
  //                                        G4ThreeVector(0, 1*cm, 0), G4ThreeVector(1*cm, 0, 0)), G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP"),"LVtet");
  // // lv_dTet->SetVisAttributes(G4VisAttributes::GetInvisible());
  // TETParameterisation* param = new TETParameterisation(tetData);
  // new G4PVParameterised("paraPara",lv_dTet, lv_phantomBox, kUndefined, tetData->GetNumTetrahedron(), param);
  // // new G4PVParameterised("param",lv_dTet, lv_phantomBox, kUndefined, 1, param);
}

void ParallelPhantom::VisualizePhantom()
{
  if (isInstall) {
		for (auto itr: materialMap) {
			G4int idx = itr.first;
			delete lvTessMap[idx];
			delete pvTessMap[idx];
		}
		delete pv_doctor;
		lvTessMap.clear();
		pvTessMap.clear();
	}
	
	G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
	G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin()) * 0.5 + G4ThreeVector(1.,1.,1.) * cm;
	G4VSolid* sol_doctorBox = new G4Box("sol_doctorBox", halfSize.x(), halfSize.y(), halfSize.z());
	G4LogicalVolume* lv_doctorBox = new G4LogicalVolume(sol_doctorBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "lv_doctorBox");
	G4VisAttributes* va_doctorBox = new G4VisAttributes();
	va_doctorBox->SetForceWireframe(true);
	lv_doctorBox->SetVisAttributes(va_doctorBox);
	pv_doctor = new G4PVPlacement(0,center,lv_doctorBox,"lv_doctorBox",lv_world,false,0);

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
      if (12200 <= idx && idx <= 12500)
        lvTessMap[idx]->SetVisAttributes(new G4VisAttributes(colourMap[12200]));
			pvTessMap[idx] = new G4PVPlacement(0, G4ThreeVector(), lvTessMap[idx], to_string(idx)+"_"+organName, lv_doctorBox, false, 0);
			cout << "...done" << endl;
		}
	}
	isInstall = true;
}

#include "G4GeometryManager.hh"
void ParallelPhantom::Deform(RotationList vQ, Vector3d root)
{
  tetData->Deform(vQ, root);
  if (!ui) {
    G4ThreeVector center = (tetData->GetPhantomBoxMax() + tetData->GetPhantomBoxMin())*0.5;
    G4ThreeVector halfSize = (tetData->GetPhantomBoxMax() - tetData->GetPhantomBoxMin())*0.5 + G4ThreeVector(5., 5., 5.)*cm; //5-cm-margin
    pv_doctor->SetTranslation(center);
    dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetXHalfLength(halfSize.x());
    dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetYHalfLength(halfSize.y());
    dynamic_cast<G4Box*>(pv_doctor->GetLogicalVolume()->GetSolid())->SetZHalfLength(halfSize.z());
  }
  else {
    VisualizePhantom();
  }
  G4cout << "GeometryHasBeenModified" << G4endl;
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

void ParallelPhantom::ConstructSD()
{
  // Define detector and scorer
	G4SDManager* SDman = G4SDManager::GetSDMpointer();

	// Detector
	G4MultiFunctionalDetector* mfd_doctor = new G4MultiFunctionalDetector("MFD_doctor");
	SDman->AddNewDetector( mfd_doctor );

	// Scorer
	mfd_doctor->RegisterPrimitive(new TETPSEnergyDeposit("eDep", tetData));
	mfd_doctor->RegisterPrimitive(new TETDRFScorer("DRF", tetData));
	
	// Attach the detector to logical volume
	if(!ui) SetSensitiveDetector(lv_dTet, mfd_doctor);
}