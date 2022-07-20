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

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4AccumulableManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4NistManager.hh"
#include "DetectorConstruction.hh"
#include "parallelmesh.hh"

RunAction::RunAction()
 : G4UserRunAction()
{
    const G4double milligray = 1.e-3*gray;
    const G4double microgray = 1.e-6*gray;
    const G4double nanogray  = 1.e-9*gray;
    const G4double picogray  = 1.e-12*gray;

    new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
    new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
    new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
    new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);
}

RunAction::~RunAction()
{}

G4Run* RunAction::GenerateRun()
{
    static_cast<ParallelMesh*>(G4RunManager::GetRunManager()->GetUserDetectorConstruction()->GetParallelWorld(0))->GetIJK(ni,nj,nk);
    // generate run
    // fRun = new Run(ni*nj*nk);
    fRun = new Run(ni*nj);
    return fRun;
}

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;
  nps=run->GetNumberOfEventToBeProcessed();
  // nps = 100000000;
  G4RunManager::GetRunManager()->SetPrintProgress(int(nps*0.1));
}

void RunAction::EndOfRunAction(const G4Run* )
{
    if(!IsMaster()) return;
    auto doseMapS=fRun->GetDoseMapS();
    // auto doseMapV=fRun->GetDoseMapV();
    auto doseMapE=fRun->GetDoseMapE();
    // auto doseMapL=fRun->GetDoseMapL();
    // std::ofstream hist(to_string(fRun->GetRunID())+".hist");
    // for(size_t i=0;i<doseMapL->size();i++) hist <<i<<" "<<(*doseMapL)[i]<<endl;
    // hist.close();
    return;

    std::vector<G4float> resultMap(doseMapS->size());
    std::ofstream ofs(std::to_string(fRun->GetRunID())+"AP.map", std::ios::binary);
    std::transform(doseMapS->begin(), doseMapS->end(), resultMap.begin(), [&](G4float d)->G4float{return d/(G4float)nps/gray;});
    ofs.write((char*) (&resultMap[0]), ni*nj*nk*sizeof(G4float));
    // std::transform(doseMapE->begin(), doseMapE->end(), resultMap.begin(), [&](G4float d)->G4float{return d/(G4float)nps/gray;});
    // ofs.write((char*) (&resultMap[0]), ni*nj*nk*sizeof(G4float));
    ofs.close();

    // std::ofstream ofsV("avg_dir.txt");
    // for(G4ThreeVector v:*doseMapV)
    // {
    //     G4ThreeVector u = v.unit();
    //     ofsV<<u.x()<<" "<<u.y()<<" "<<u.z()<<endl;
    // }
    // ofsV.close();

    // std::ofstream ofs("energy_analysis.txt");
    // for(G4int i=0;i<doseMapS->size();i++) ofs<<(*doseMapS)[i]/(*doseMapL)[i]/keV<<G4endl;
    // // for(G4float f:(*doseMapS)) ofs<<f/(doseMapL)/keV<<endl;
    // ofs.close();
    // std::ofstream ofs1("cosine_analysis.txt");
    // for(G4int i=0;i<doseMapE->size();i++) ofs1<<(*doseMapE)[i]/(*doseMapL)[i]<<G4endl;
    // ofs1.close();
    //dose map unit: (/cm2)
//    G4cout<<"dose of DAP meter: "<<G4BestUnit(fRun->GetDap()/nps,"Dose")<<G4endl;
}
