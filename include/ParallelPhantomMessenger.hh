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
// ParallelMessenger.cc
// \author Haegin Han
//

#ifndef SRC_ParallelMessenger_HH_
#define SRC_ParallelMessenger_HH_ 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "functions.h"

class G4UIdirectory;
class G4UIcmdWithAnInteger;
class G4UIcmdWithAString;
class ParallelPhantom;

typedef
  std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >
  VectorList;
class ParallelPhantomMessenger: public G4UImessenger
{
public:
	ParallelPhantomMessenger(ParallelPhantom* phantom);
	virtual ~ParallelPhantomMessenger();

	virtual void SetNewValue(G4UIcommand*, G4String);
	// void ReadPostureData(G4String fileName);

private:
	ParallelPhantom* fPhantom;
	G4UIdirectory*   fPhantomDir;
	G4UIcmdWithAString* fInitCmd; 
	G4UIcmdWithAString* fDeformCmd; 
	G4UIcmdWithAnInteger* fSetFrameCmd; 
	
};

#endif
