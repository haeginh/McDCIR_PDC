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
#include "ParallelGlass.hh"
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
#include "G4UImanager.hh"

ParallelGlass
::ParallelGlass(G4String parallelWorldName, G4String _glassFile)
:G4VUserParallelWorld(parallelWorldName), glassFile(_glassFile),fConstructed(false), focalLength(810*mm)
{
}

ParallelGlass::~ParallelGlass()
{
}

void ParallelGlass::Construct()
{
  if(fConstructed) return;
  fConstructed = true;
 
  G4String dir = std::getenv("DCIR_PHANTOM_DIR");
  ifstream ifsNode(dir+glassFile + ".node");
  if(!ifsNode.is_open())
    G4Exception("ParallelGlass::Contruct", "", FatalException, ("cannot open"+dir+glassFile + ".node").c_str());
  G4int num, tmp;
  ifsNode>>num>>tmp>>tmp>>tmp;
  vector<G4ThreeVector> nodes;
  G4ThreeVector glassDim;
  for(G4int i=0;i<num;i++)
  {
    G4ThreeVector v;
    ifsNode>>tmp>>v;
    v *= cm;
    nodes.push_back(v);
    if(glassDim.x()<fabs(v.x())) glassDim.setX(fabs(v.x()));
    if(glassDim.y()<fabs(v.y())) glassDim.setY(fabs(v.y()));
    if(glassDim.z()<fabs(v.z())) glassDim.setZ(fabs(v.z()));
  }
  ifsNode.close();

  G4Box* glassSolid = new G4Box("glass", glassDim.x(), glassDim.y(), glassDim.z());
  lv_glass = new G4LogicalVolume(glassSolid, 0, "glassLV");
  glass_rotation_matrix = new G4RotationMatrix;
  pv_glass = new G4PVPlacement(glass_rotation_matrix, G4ThreeVector(), lv_glass, "glassPV", GetWorld()->GetLogicalVolume(), false, 0);
  G4VisAttributes* vis = new G4VisAttributes();
  vis->SetForceWireframe(true);
  pv_glass->GetLogicalVolume()->SetVisAttributes(vis);

  ifstream ifsEle(dir+glassFile + ".ele");
    G4Exception("ParallelGlass::Contruct", "", FatalException, ("cannot open"+dir+glassFile + ".ele").c_str());
  ifsEle>>num>>tmp>>tmp;
  G4int a, b, c, d;
  G4Material* lead = G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb");
  G4VisAttributes* vis2 = new G4VisAttributes();
  vis2->SetColor(G4Color(0., 1., 0., 0.2));
  for(G4int i=0;i<num;i++)
  {
    ifsEle>>tmp>>a>>b>>c>>d;
    G4Tet* tet = new G4Tet("glass",nodes[a], nodes[b], nodes[c], nodes[d]);
    auto logical = new G4LogicalVolume(tet, lead, "glassTet");
    logical->SetVisAttributes(vis2);
    new G4PVPlacement(0, G4ThreeVector(), logical, "glassTetPV", pv_glass->GetLogicalVolume(), false, 0);
  }
  ifsEle.close();

  ConstructCarmDet();
}

void ParallelGlass::ConstructCarmDet()
{
	// This member function didn't considered the set functions of G4VPhyscialVolume.
	carm_rotation_matrix = new G4RotationMatrix;

	// C-arm	
	lv_det = new G4LogicalVolume(new G4Box("pv_det", 42*cm*0.5, 52*cm*0.5, 1*cm), G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_det");
	lv_det->SetVisAttributes(new G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.5)));
	pv_det = new G4PVPlacement(carm_rotation_matrix, G4ThreeVector(), lv_det, "pv_det", GetWorld()->GetLogicalVolume(), false, 0);
}
