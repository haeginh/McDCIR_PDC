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

#include "Run.hh"

Run::Run(G4int n)
:G4Run()//, dap(0)
{
    fCollID_skin
    = G4SDManager::GetSDMpointer()->GetCollectionID("meshSD/doseS");
    // fCollID_vector
    // = G4SDManager::GetSDMpointer()->GetCollectionID("meshSD/vector");
    fCollID_lens
    = G4SDManager::GetSDMpointer()->GetCollectionID("meshSD/doseE");
    // fCollID_length
    // = G4SDManager::GetSDMpointer()->GetCollectionID("meshSD/doseL");
//    fCollID_dap
//    = G4SDManager::GetSDMpointer()->GetCollectionID("dap/dose");

    doseMapS.resize(n,0);
    // doseMapV.resize(n,G4ThreeVector());
    doseMapE.resize(n,0);
    // doseMapL.resize(120,0);
}

Run::~Run()
{}

void Run::RecordEvent(const G4Event* event)
{
    // Hits collections
    //
    G4HCofThisEvent* HCE = event->GetHCofThisEvent();
    if(!HCE) return;

    auto _doseMapS =
            *static_cast<G4THitsMap<G4float>*>(HCE->GetHC(fCollID_skin))->GetMap();
    // auto _doseMapV =
    //         *static_cast<G4THitsMap<G4ThreeVector>*>(HCE->GetHC(fCollID_vector))->GetMap();
    auto _doseMapE =
            *static_cast<G4THitsMap<G4float>*>(HCE->GetHC(fCollID_lens))->GetMap();
    // auto _doseMapL =
    //         *static_cast<G4THitsMap<G4float>*>(HCE->GetHC(fCollID_length))->GetMap();

    for(auto itr:_doseMapS){
        doseMapS[itr.first]+= *itr.second;
        // doseMapV[itr.first]+= *_doseMapV[itr.first];
        doseMapE[itr.first]+= *_doseMapE[itr.first];
        // doseMapL[itr.first]+= *_doseMapL[itr.first];
    }
//    auto dapMap =
//            *static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_dap))->GetMap();
//    if(dapMap.find(1000)!=dapMap.end()) dap += *dapMap[1000];

}

void Run::Merge(const G4Run* run)
{
    const Run* localRun = static_cast<const Run*>(run);
    // merge the data from each thread
    auto localMapS = localRun->doseMapS;
    // auto localMapV = localRun->doseMapV;
    auto localMapE = localRun->doseMapE;
    // auto localMapL = localRun->doseMapL;

    for(size_t i=0;i<doseMapS.size(); i++){
        doseMapS[i]  += localMapS[i];
        // doseMapV[i]  += localMapV[i];
        doseMapE[i]  += localMapE[i];
        // doseMapL[i]  += localMapL[i];
    }
//    dap += localRun->dap;
    factor = localRun->factor;
    G4Run::Merge(run);
}
