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
#include "G4Tet.hh"
#include "PhantomAnimator.hh"

class G4LogicalVolume;
class G4PVParameterised;
class TETParameterisation;
class ParallelPhantomMessenger;

class ParallelPhantom : public G4VUserParallelWorld
{
public:
  ParallelPhantom(G4String parallelWorldName, PhantomData* amData, PhantomData* _afData);
  virtual ~ParallelPhantom();

public:
  virtual void Construct();
  virtual void ConstructSD();
  void InitializePhantoms(G4String fname);
  void Deform(G4int phantomID, RotationList vQ, Vector3d root);
  void SetNewFrame(G4int phantomList); //binary representation: ...(2)(1)(0)
  G4Tet* GetTet(G4int copyNo)
  {
    auto phantom = thisFrame.upper_bound(copyNo);
    return tetVec[phantom->second][copyNo - phantom->first];
  }
  G4Material* GetMat(G4int copyNo)
  {
    auto phantom = thisFrame.upper_bound(copyNo);
    auto idx = phantomAnimators[phantom->second]->GetMaterialIdx(copyNo - phantom->first);
    return dataVec[phantom->second]->GetMaterial(idx);
  }
  pair<G4int, G4int> GetIdx(G4int copyNo)
  {
    auto phantom = thisFrame.upper_bound(copyNo);
    auto idx = phantomAnimators[phantom->second]->GetMaterialIdx(copyNo - phantom->first);
    return make_pair(phantom->second, idx);
  }
  G4bool GetBoneMassRatio(G4int phantomID, G4int organID, G4double &rbmR, G4double &bsR)
  {
    return dataVec[phantomID]->GetBoneMassRatio(organID, rbmR, bsR);
  }
  PhantomData* GetData(G4int idx) {return dataVec[idx];}
  G4int GetNumberOfPhantoms(){return phantomAnimators.size();}
  std::map<G4int, G4double> GetMassMap(G4int id){return massMaps[id];}
  G4double GetVol(G4int phantomID, G4int organID) {return volMaps[phantomID][organID];}

private:
  G4bool fConstructed;
  ParallelPhantomMessenger* messenger;

  // Radiologist
  PhantomData*           amData;
  PhantomData*           afData;
  std::vector<PhantomAnimator*> phantomAnimators;
  std::vector<std::vector<G4Tet*>> tetVec;
  std::vector<PhantomData*> dataVec;
  std::vector<std::map<G4int, G4double>> massMaps, volMaps;
  std::map<int, int> thisFrame;
  // G4ThreeVector      doctor_translation;
  G4LogicalVolume*   lv_tet;
  TETParameterisation* tetParam;
  G4PVParameterised* paramPV;
  G4VPhysicalVolume* pv_doctor;
};

#endif
