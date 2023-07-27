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

#include "G4Box.hh"
#include "G4Tet.hh"
#include "G4LogicalVolume.hh"
#include "G4PVParameterised.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"    
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "TETPSEnergyDeposit.hh"
#include "ParallelPhantomMessenger.hh"
#include "G4UImanager.hh"
#include "TETParameterisation.hh"
#include "BoneScorer.hh"


ParallelPhantom
::ParallelPhantom(G4String parallelWorldName, PhantomData* _amData, PhantomData* _afData)
:G4VUserParallelWorld(parallelWorldName),fConstructed(false), amData(_amData), afData(_afData)
{
  messenger = new ParallelPhantomMessenger(this);
}

ParallelPhantom::~ParallelPhantom()
{
  delete messenger;
}

void ParallelPhantom::Construct()
{
  if(fConstructed) return;
  fConstructed = true;
 
  //
  // World
  //
  G4VPhysicalVolume* ghostWorld = GetWorld();
  G4LogicalVolume* worldLogical = ghostWorld->GetLogicalVolume();

  //
  // parallel world placement box
  //
  // G4ThreeVector center = (PhantomData->GetPhantomBoxMax() + PhantomData->GetPhantomBoxMin())*0.5;
  // G4ThreeVector halfSize = (PhantomData->GetPhantomBoxMax() - PhantomData->GetPhantomBoxMin())*0.5 + G4ThreeVector(5., 5., 5.)*cm; //5-cm-margin
  MatrixXd box;
  for(auto phantom:phantomAnimators)
  {
    box.conservativeResize(box.rows()+2, 3);
    box.bottomRows(2) = phantom->GetBox();
  }
  RowVector3d max = box.colwise().maxCoeff();
  RowVector3d min = box.colwise().minCoeff();
  RowVector3d halfSize = (max-min)*0.5;
  RowVector3d center = (max+min)*0.5;
  G4VSolid* paraBox = new G4Box("phantomBox",halfSize(0)+10*cm,halfSize(1)+10*cm,halfSize(2)+10*cm);
  G4LogicalVolume* lv_phantomBox = new G4LogicalVolume(paraBox,0,"phantomBox");
  lv_phantomBox->SetVisAttributes(G4VisAttributes::GetInvisible());
  pv_doctor = new G4PVPlacement(0,G4ThreeVector(center(0),center(1),center(2)),lv_phantomBox,"phantomBox",worldLogical,false,0);

  //
  // mother of parallel world parameterized volumes
  //

  // G4Material* tissue = G4Material::GetMaterial("G4_TISSUE_SOFT_ICRP");
  // G4cout<<PhantomData->GetNumTetrahedron()<<G4endl;
  // for(size_t i=0;i<PhantomData->GetNumTetrahedron();i++)
  // {
  //   new G4PVPlacement(0,G4ThreeVector(),new G4LogicalVolume(PhantomData->GetTetrahedron(i),tissue, "tet"),"tet",worldLogical,false,0);
  // }
  lv_tet = new G4LogicalVolume(new G4Tet("tet", G4ThreeVector(), G4ThreeVector(0, 0, 1*cm),
                                         G4ThreeVector(0, 1*cm, 0), G4ThreeVector(1*cm, 0, 0)), G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP"),"LVtet");
  // lv_tet->SetVisAttributes(G4VisAttributes::GetInvisible());
  tetParam = new TETParameterisation(this);
  // paramPV = new G4PVParameterised("docParam",lv_tet, lv_phantomBox, kUndefined, PhantomData->GetNumTetrahedron(), param);
  // new G4PVParameterised("param",lv_tet, lv_phantomBox, kUndefined, 1, param);
}

void ParallelPhantom::ConstructSD()
{
  G4MultiFunctionalDetector* mfd = new G4MultiFunctionalDetector("phantom");
  G4SDManager::GetSDMpointer()->AddNewDetector(mfd);
  // G4VPrimitiveScorer* scorer = new TETPSEnergyDeposit("edep", this);
  // mfd->RegisterPrimitive(scorer);
  mfd->RegisterPrimitive(new G4PSEnergyDeposit("edep"));
  mfd->RegisterPrimitive(new BoneScorer("drf", this));
  SetSensitiveDetector(lv_tet, mfd);
}

void ParallelPhantom::InitializePhantoms(G4String name)
{
  if(phantomAnimators.size()) return;
  
  ifstream ifs(name);
  if(!ifs.is_open())
  {
   	G4Exception("ParallelPhantom::InitializePhantoms","",FatalErrorInArgument,
		G4String("      There is no " + name ).c_str());
  }
  vector<pair<G4String, G4String>> phantoms;
  while(1)
  {
    G4String phantom, profile, sex;
    ifs>>phantom>>profile>>sex;
    if(phantom.empty()) continue;
    phantoms.push_back(make_pair(phantom, profile));
    if(sex=="M"||sex=="m") dataVec.push_back(amData);
    else if(sex=="F"||sex=="f") dataVec.push_back(afData);
    else{
     	G4Exception("ParallelPhantom::InitializePhantoms","",FatalErrorInArgument,
	  	G4String("      wrong sex id: " + sex ).c_str());
    }
  }
  ifs.close();
  
  massMaps.resize(phantoms.size());
  volMaps.resize(phantoms.size());
  for(G4int p=0;p<phantoms.size();p++)
  {
    PhantomAnimator* animator = new PhantomAnimator(phantoms[p].first, phantoms[p].second);
    std::vector<G4Tet*> tet;
    for(G4int i=0;i<animator->GetNumOfTet();i++)
    {
      G4ThreeVector a, b, c, d;
      G4bool chk;
      animator->GetTetVertices(i, a, b, c, d);    
      auto tet0 = new G4Tet(to_string(phantomAnimators.size())+"_tet#"+to_string(i),
                  a, b, c, d, &chk); 
      tet.push_back(tet0);
      G4int id = animator->GetMaterialIdx(i);
      massMaps[p][id] 
        += tet0->GetCubicVolume()*dataVec[p]->GetMaterial(id)->GetDensity();
      volMaps[p][id] 
        += tet0->GetCubicVolume();
    }
    phantomAnimators.push_back(animator);
    tetVec.push_back(tet);
  }
}

void ParallelPhantom::Deform(G4int phantomID, RotationList vQ, Vector3d root)
{
  if(phantomID>=phantomAnimators.size())
  {
    G4cerr<<"wront phantom id for ParallelPhantom::Deform"<<G4endl;
    return;
  }
  phantomAnimators[phantomID]->Animate(vQ, root);
}

void ParallelPhantom::SetNewFrame(G4int phantomList)
{
  if(phantomList<0 || phantomList>=1<<phantomAnimators.size())
  {
    G4cerr<<"Wrong phantomList number ("<<phantomList
          <<"), the number of the init. phantom: "<<phantomAnimators.size()<<G4endl;
    exit(10);
  }
  thisFrame.clear();
  G4int num(0);
  MatrixXd box;
  for(G4int i=0;i<phantomAnimators.size();i++)
  {
    if(!(phantomList & 1<<i)) continue;
    num += phantomAnimators[i]->GetNumOfTet();
    thisFrame[num] = i;
    box.conservativeResize(box.rows()+2, 3);
    box.bottomRows(2) = phantomAnimators[i]->GetBox();
  }
  RowVector3d max = box.colwise().maxCoeff();
  RowVector3d min = box.colwise().minCoeff();
  RowVector3d halfSize = (max-min)*0.5;
  RowVector3d center = (max+min)*0.5;
  G4ThreeVector g4center(center(0), center(1), center(2));
  pv_doctor->SetTranslation(g4center);
  ((G4Box*)pv_doctor->GetLogicalVolume()->GetSolid())->SetXHalfLength(halfSize(0)+10*cm);
  ((G4Box*)pv_doctor->GetLogicalVolume()->GetSolid())->SetYHalfLength(halfSize(1)+10*cm);
  ((G4Box*)pv_doctor->GetLogicalVolume()->GetSolid())->SetZHalfLength(halfSize(2)+10*cm);

  G4bool chk;
  for(G4int i=0;i<phantomAnimators.size();i++)
  {
    if(!(phantomList & 1<<i)) continue;
    massMaps[i].clear();
    volMaps[i].clear();
    for(G4int t=0;t<phantomAnimators[t]->GetNumOfTet();t++)
    {
      G4ThreeVector a, b, c, d;
      phantomAnimators[i]->GetTetVertices(t, a, b, c, d); 
      tetVec[i][t]->SetVertices(a-g4center, b-g4center, c-g4center, d-g4center, &chk);
      G4int idx = phantomAnimators[i]->GetMaterialIdx(t);
      massMaps[i][idx] 
        += tetVec[i][t]->GetCubicVolume()*dataVec[i]->GetMaterial(idx)->GetDensity();
      volMaps[i][idx] 
        += tetVec[i][t]->GetCubicVolume();
    }
  }

  if(paramPV) delete paramPV;
  paramPV = new G4PVParameterised("phantomParam",lv_tet, pv_doctor->GetLogicalVolume(), kUndefined, num, tetParam);
}