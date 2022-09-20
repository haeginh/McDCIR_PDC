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
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "TETPSEnergyDeposit.hh"
#include "BoneScorer.hh"
// #include "ParallelPhantomMessenger.hh"
#include "G4UImanager.hh"
#include "TETParameterisation.hh"


ParallelPhantom
::ParallelPhantom(G4String parallelWorldName, PhantomData* _phantom)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false), phantomData(_phantom), phantomBoxPhy(0), g4ParamPhy(0), phantomInstalled(false)
{
  // messenger = new ParallelPhantomMessenger(this);
  phantomAnimator = new PhantomAnimator;
  if(!getenv("DCIR_PHANTOM_DIR")){
    G4Exception("ParallelPhantom::ParallelPhantom","",FatalErrorInArgument,
				        G4String("Set env DCIR_PHANTOM_DIR").c_str());
  }
  phantomDir = getenv("DCIR_PHANTOM_DIR");
  doseMass = Eigen::VectorXd::Zero(phantomData->organ2dose.cols());
}

ParallelPhantom::~ParallelPhantom()
{
  // delete messenger;
  delete phantomAnimator;
}

void ParallelPhantom::Construct()
{
  if(fConstructed) return;
  fConstructed = true;
 
  //
  // prepare basic geometries
  //
  GetWorld()->GetLogicalVolume()->SetVisAttributes(G4VisAttributes::GetInvisible());
  tetLog = new G4LogicalVolume(new G4Tet("tet", G4ThreeVector(), G4ThreeVector(0, 0, 1*cm),
                                         G4ThreeVector(0, 1*cm, 0), G4ThreeVector(1*cm, 0, 0)), G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP"),"LVtet");
  G4VSolid* phantomBoxSol = new G4Box(GetName() + "_phantomBox",1*m, 1*m, 1*m);
  phantomBoxLog = new G4LogicalVolume(phantomBoxSol,G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP"), GetName() + "_phantomBox");
  phantomBoxLog->SetSmartless(10);
  auto wireFrame = new G4VisAttributes();
  wireFrame->SetForceWireframe(true);
  phantomBoxLog->SetVisAttributes(wireFrame);
  phantomParamPhy = new TETParameterisation(this);

  // load and install phantom
  LoadPhantom(phantomDir + "AM_BMI26", "original");
  // InstallPhantomBox();
}

void ParallelPhantom::ConstructSD()
{
  G4MultiFunctionalDetector* mfd = new G4MultiFunctionalDetector(GetName()+"_MFD");
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  G4VPrimitiveScorer* scorer = new TETPSEnergyDeposit("edep", this);
  G4VPrimitiveScorer* boneScorer = new BoneScorer("bone", this);
  mfd->RegisterPrimitive(scorer);
  mfd->RegisterPrimitive(boneScorer);
  SetSensitiveDetector(tetLog, mfd);
}

 void ParallelPhantom::LoadPhantom(G4String name, G4String scale){
    phantomAnimator->LoadPhantom(name);
    // phantomAnimator->CalibrateTo(scale);
    RowVector3d max = phantomAnimator->U.colwise().maxCoeff();
    RowVector3d min = phantomAnimator->U.colwise().minCoeff();
    G4ThreeVector g4max = G4ThreeVector(max(0), max(1), max(2))*cm;
    G4ThreeVector g4min = G4ThreeVector(min(0), min(1), min(2))*cm;
    boxCenter = (g4max + g4min)*0.5;
    boxHalfSize = (g4max - g4min)*0.5*1.1;

    numOfTet = phantomAnimator->T.rows();
    G4bool degenChk;
    G4int count(0);
    organVol.clear();
    numOfSkinFacet = 0;
    for(G4int i=0;i<numOfTet;i++)
    {
      G4ThreeVector a, b, c, d;
      phantomAnimator->GetThreeVector(i, a, b, c, d);
      if(i<tetVec.size())
        tetVec[i]->SetVertices(a-boxCenter, b-boxCenter, c-boxCenter, d-boxCenter, &degenChk);
      else
        tetVec.push_back(new G4Tet(GetName()+"_tet",a-boxCenter, b-boxCenter, c-boxCenter, d-boxCenter, &degenChk));
      if(degenChk) count++;
      organVol[phantomAnimator->T(i, 4)]+=tetVec[i]->GetCubicVolume();
      if(phantomAnimator->T(i, 4)<=0) numOfSkinFacet++;
    }

    doseMass.setZero();
    for(G4int i=0;i<phantomData->organ2dose.cols();i++) //organ
    {
      for(G4int j=0;j<phantomData->organ2dose.rows();j++) //dose
      {
        if(phantomData->organ2dose(j,i)==0) continue;
        doseMass(j)+=organVol[i] * phantomData->GetMaterial(i)->GetDensity();
      }
    }
    doseMass(0) = 1;
    doseMass(1) = 1;
    boneMassInv.clear();
    auto vec = phantomData->GetBoneIdxVec();
    for(G4int i:vec)
      boneMassInv[i] = 1./(organVol[i] * phantomData->GetMaterial(i)->GetDensity());
  }
  
bool ParallelPhantom::Deform(RotationList vQ, RowVector3d rootC)
  {
    phantomAnimator->Animate(vQ, rootC);

    RowVector3d max = phantomAnimator->U.colwise().maxCoeff();
    RowVector3d min = phantomAnimator->U.colwise().minCoeff();
    G4ThreeVector g4max = G4ThreeVector(max(0), max(1), max(2))*cm;
    G4ThreeVector g4min = G4ThreeVector(min(0), min(1), min(2))*cm;
    boxCenter = (g4max + g4min)*0.5;
    boxHalfSize = (g4max - g4min)*0.5*1.1;

    numOfTet = phantomAnimator->T.rows();
    G4bool degenChk;
    G4int count(0);
    organVol.clear();
    for(G4int i=0;i<numOfTet;i++)
    {
      G4ThreeVector a, b, c, d;
      phantomAnimator->GetThreeVector(i, a, b, c, d);
      tetVec[i]->SetVertices(a-boxCenter, b-boxCenter, c-boxCenter, d-boxCenter, &degenChk);
      if(degenChk) count++;
      organVol[phantomAnimator->T(i, 4)]+=tetVec[i]->GetCubicVolume();
    }
    if(!phantomBoxPhy) return false;
    else
    {
      auto box = (G4Box*) (phantomBoxLog->GetSolid());
      box->SetXHalfLength(boxHalfSize(0));
      box->SetYHalfLength(boxHalfSize(1));
      box->SetZHalfLength(boxHalfSize(2));
      phantomBoxPhy->SetTranslation(boxCenter);
    }

    doseMass.setZero();
    for(G4int i=0;i<phantomData->organ2dose.cols();i++)
    {
      for(G4int j=0;j<phantomData->organ2dose.rows();j++)
      {
        if(phantomData->organ2dose(j,i)==0) continue;
        doseMass(j)+=organVol[i] * phantomData->GetMaterial(i)->GetDensity();
      }
    }
    doseMass(0) = 1;
    doseMass(1) = 1;
    boneMassInv.clear();
    auto vec = phantomData->GetBoneIdxVec();
    for(G4int i:vec)
      boneMassInv[i] = 1./(organVol[i] * phantomData->GetMaterial(i)->GetDensity());
  }

bool ParallelPhantom::InstallPhantomBox()
  {
    auto box = (G4Box*) (phantomBoxLog->GetSolid());
    box->SetXHalfLength(boxHalfSize(0));
    box->SetYHalfLength(boxHalfSize(1));
    box->SetZHalfLength(boxHalfSize(2));
    if(phantomInstalled) return false;
    g4ParamPhy = new G4PVParameterised(GetName()+"_paramPhy", tetLog, phantomBoxLog, kUndefined, phantomAnimator->T.rows(), phantomParamPhy);
    // g4ParamPhy = new G4PVParameterised(GetName()+"_paramPhy", tetLog, phantomBoxLog, kUndefined, 1000, phantomParamPhy);
    phantomBoxPhy = new G4PVPlacement(0, boxCenter, phantomBoxLog, GetName()+"_phantomBox",GetWorld()->GetLogicalVolume(), false, 0);
    phantomInstalled = true;
    return true;
  }

 bool ParallelPhantom::RemovePhantomBox(){
    if(!phantomInstalled) return false;
    delete g4ParamPhy;
    delete phantomBoxPhy;
    phantomInstalled = false;
    return true;
  }

