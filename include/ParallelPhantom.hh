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
#ifndef ParallelPhantom_h
#define ParallelPhantom_h 1

#include "G4VUserParallelWorld.hh"
#include "G4PVPlacement.hh"
#include "G4RunManager.hh"
#include "globals.hh"
#include "PhantomData.hh"
#include "PhantomAnimator.hh"
#include "G4PVParameterised.hh"

class G4LogicalVolume;
class TETParameterisation;
// class ParallelPhantomMessenger;

class ParallelPhantom : public G4VUserParallelWorld
{
public:
  ParallelPhantom(G4String parallelWorldName, PhantomData* _phantom);
  virtual ~ParallelPhantom();

public:
  virtual void Construct();
  virtual void ConstructSD();

  G4Tet* GetTet(G4int idx){return tetVec[idx];}
  G4int  GetOrganID(G4int idx){return phantomAnimator->GetID(idx);}
  G4VisAttributes* GetVisAtt(G4int idx){return phantomData->GetVisAtt(GetOrganID(idx));}
  G4Material* GetMateiral(G4int idx){
    G4int organID = GetOrganID(idx);
    if(organID>0)
    {
      auto mat = phantomData->GetMaterial(organID);
      if(!mat)
        G4Exception("ParallelPhantom::GetMaterial", "", JustWarning, ("Wrong organID: "+to_string(organID)).c_str());
      return mat;
    }
    else return phantomData->GetMaterial(126);
  }
  G4bool IsForVis(){return phantomData->IsForVis();}
  G4bool GetBoneMassRatio(G4int organID, G4double &rbmR, G4double &bsR)
  {
    return phantomData->GetBoneMassRatio(organID, rbmR, bsR);
  }
 	void GetDRF(G4double energy, G4int organID, G4double &rbmDRF, G4double &bsDRF)
  {
    phantomData->GetDRF(energy, organID, rbmDRF, bsDRF);
  }
  G4double GetVolume(G4int idx){return organVol[idx];}
  G4double GetBoneMassInv(G4int idx){return boneMassInv[idx];}
  G4int GetNumOfSkinFacet() {return numOfSkinFacet;}
  // G4double GetRBMDRF(G4int idx, G4int eIdx){ return phantomData->GetRBMDRF(idx, eIdx);}
	// G4double GetBSDRF (G4int idx, G4int eIdx){ return phantomData->GetBSDRF(idx, eIdx);}

  void LoadPhantom(G4String name, G4String scale);
  bool Deform(RotationList vQ, RowVector3d root); // return true if phantom was installed
  bool InstallPhantomBox();
  bool RemovePhantomBox();
  
  PhantomData*       phantomData;
  Eigen::ArrayXd    doseMass;
private:
  G4String phantomDir;
  G4bool fConstructed;

  // Radiologist
  G4int              numOfTet;
  G4int              numOfSkinFacet;
  PhantomAnimator*   phantomAnimator;
  G4ThreeVector      doctor_translation;
  G4VPhysicalVolume* pv_doctor;
  vector<G4Tet*>     tetVec;
  std::map<G4int, G4double> organVol;
  std::map<G4int, G4double> boneMassInv;
  G4ThreeVector      boxHalfSize, boxCenter;
  G4PVPlacement*     phantomBoxPhy;
  TETParameterisation* phantomParamPhy;  
  G4PVParameterised*   g4ParamPhy;
  G4LogicalVolume*   phantomBoxLog;
  G4LogicalVolume*   tetLog;
  G4double rbmMass, bsMass;

  G4bool phantomInstalled;
};

#endif
