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
#ifndef ParallelGlass_h
#define ParallelGlass_h 1

#include "G4VUserParallelWorld.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "G4GeometryManager.hh"
#include "G4SystemOfUnits.hh"
#include <Eigen/Geometry>

class G4LogicalVolume;

class ParallelGlass : public G4VUserParallelWorld
{
public:
  ParallelGlass(G4String parallelWorldName, G4String _galssFile);
  virtual ~ParallelGlass();

  virtual void Construct();

  // Glass
  void SetGlassPose(G4ThreeVector trans, Eigen::Quaterniond q)
  {
    if(pv_glass == NULL)
    {
      glass_rotation_matrix = new G4RotationMatrix();
      pv_glass = new G4PVPlacement(glass_rotation_matrix, G4ThreeVector(), lv_glass, "glassPV",GetWorld()->GetLogicalVolume(), false, 0);
    }

    // G4GeometryManager::GetInstance()->OpenGeometry(pv_glass);
    Eigen::AngleAxisd rot(q.inverse());
    glass_rotation_matrix->setTheta(0);
    glass_rotation_matrix->setAxis(G4ThreeVector(rot.axis()(0), rot.axis()(1), rot.axis()(2)));
    glass_rotation_matrix->setTheta(rot.angle()*rad);
    glass_rotation_matrix->invert();
    pv_glass->SetTranslation(trans);
    // G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_glass);
  }

  void RemoveGlass()
  {
    if(pv_glass == nullptr) return;
    // G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_glass);
    delete glass_rotation_matrix;
    delete pv_glass;
    pv_glass = nullptr;
  }

  // C-arm det
  void SetCarmDetPose(G4double _carm_primary, G4double _carm_secondary, G4double _carm_sid) 
	// 	carm_primary: rao, lao | carm_secondary: cran, caud
	{
    if(pv_det == NULL)
    {
      carm_rotation_matrix = new G4RotationMatrix;
      pv_det = new G4PVPlacement(carm_rotation_matrix, G4ThreeVector(), lv_det, "pv_det", GetWorld()->GetLogicalVolume(), false, 0);
    }
    
		// G4GeometryManager::GetInstance()->OpenGeometry(pv_det);
		carm_rotation_matrix->setTheta(0); //set identity
		carm_rotation_matrix->rotateY(_carm_primary).rotateX(_carm_secondary);
		carm_rotation_matrix->invert();
		pv_det->SetTranslation(carm_rotation_matrix->inverse()*G4ThreeVector(0,0,_carm_sid-focalLength));
		// G4GeometryManager::GetInstance()->CloseGeometry(false, false, pv_det);
  }

  void RemoveCarmDet()
  {
    if(pv_det == NULL) return;
    delete carm_rotation_matrix;
    delete pv_det;
    pv_det;
  }
  void ConstructCarmDet();

private:
  G4String glassFile;
  G4bool fConstructed;

  // Glass
  G4LogicalVolume* lv_glass;
  G4VPhysicalVolume* pv_glass;
	G4RotationMatrix* carm_rotation_matrix; //inverse

	// C-arm
  G4LogicalVolume* lv_det;
	G4VPhysicalVolume* pv_det;
	G4RotationMatrix* glass_rotation_matrix; //inverse
	G4double focalLength;
};

#endif
