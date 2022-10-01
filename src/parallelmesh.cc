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
//

#include "parallelmesh.hh"
#include "meshsd.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4Box.hh"

#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"

//G4ThreadLocal G4bool ParallelMesh::fSDConstructed = false;

ParallelMesh::ParallelMesh(G4String worldName)
:G4VUserParallelWorld(worldName)
{
}

ParallelMesh::~ParallelMesh()
{}

void ParallelMesh::Construct()
{
    auto vis = new G4VisAttributes;
    vis->SetForceWireframe(true);
    //mesh tally box
    G4double meshHalfX = 2.*m, meshHalfY = 2.*m, meshHalfZ = 1.5*m;
  	G4VSolid* meshSolid = new G4Box("meshSolid", meshHalfX, meshHalfY, meshHalfZ);
  	G4LogicalVolume* meshLogical = new G4LogicalVolume(meshSolid,0,"meshLogical");
    meshLogical->SetVisAttributes(vis);
    GetWorld()->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());
	new G4PVPlacement(0,G4ThreeVector(), meshLogical,"meshPhysical", GetWorld()->GetLogicalVolume(), false,0,false);

    G4double halfX=2.5*cm, halfY=2.5*cm, halfZ=2.5*cm;
    ni = meshHalfX / halfX;
    nj = meshHalfY / halfY;
    nk = meshHalfZ / halfZ;

    G4Box* boxX = new G4Box("BoxX", halfX, halfY*nj, halfZ*nk);
    G4LogicalVolume* boxX_log = new G4LogicalVolume(boxX, 0, "boxX_log");
    G4Box* boxY = new G4Box("BoxY", halfX, halfY, halfZ*nk);
    G4LogicalVolume* boxY_log = new G4LogicalVolume(boxY, 0, "boxY_log");
    G4Box* boxZ = new G4Box("BoxZ", halfX, halfY, halfZ);
    boxZ_log = new G4LogicalVolume(boxZ, 0, "boxZ_log");
    boxX_log->SetVisAttributes(vis);
    boxY_log->SetVisAttributes(vis);
    boxZ_log->SetVisAttributes(vis);

    // new G4PVPlacement(0, G4ThreeVector(), boxX_log, "skin", GetWorld()->GetLogicalVolume(), false, 1);
    new G4PVReplica("meshX", boxX_log, meshLogical, kXAxis, ni, halfX*2);
    new G4PVReplica("meshY", boxY_log, boxX_log, kYAxis, nj, halfY*2);
    new G4PVReplica("meshZ", boxZ_log, boxY_log, kZAxis, nk, halfZ*2);
}

void ParallelMesh::ConstructSD()
{
    SetSensitiveDetector(boxZ_log, new MeshSD("meshSD", ni, nj, nk, 125*cm3));
}


