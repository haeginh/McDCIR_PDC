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
// DetectorMessenger.cc
// \file   MRCP_GEANT4/External/src/TETModelImport.cc
// \author Haegin Han
//

#ifndef SRC_DetectorMessenger_HH_
#define SRC_DetectorMessenger_HH_ 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include <Eigen/Geometry>
using namespace Eigen;

typedef
  std::vector<Eigen::Quaterniond,Eigen::aligned_allocator<Eigen::Quaterniond> >
  RotationList; 
struct Body
{
    // clock_t time;
    RotationList posture = RotationList(18);
    MatrixXd jointC = MatrixXd::Zero(24, 3);
};

struct DataSet
{
    bool glassChk;
    Affine3d glass_aff = Affine3d::Identity();
    std::map<G4int, Body> bodyMap;
    RowVector3f cArm; // rot, ang, sid
    RowVector3f bed; // long, lat, height
    int kVp;
    float mA;
    int FD;
    float dap;
    bool beamOn;
    clock_t time;
    // RotationList posture = RotationList(18);
    // bool bodyIn;
    // MatrixXd jointC = MatrixXd::Zero(24, 3);
};

class G4UIdirectory;
class G4UIcmdWith3Vector;
class G4UIcmdWithAString;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWithABool;
class G4UIcmdWithAnInteger;
class DetectorConstruction;

class DetectorMessenger: public G4UImessenger
{
public:
	DetectorMessenger(DetectorConstruction* det);
	virtual ~DetectorMessenger();

	virtual void SetNewValue(G4UIcommand*, G4String);
	DataSet ReadAFrameData(std::vector<G4float>);

private:
	DetectorConstruction* fDet;
	G4UIdirectory*        fMachineDir;

	G4UIcmdWith3VectorAndUnit* fTableRefCmd; //table ref. rel. to iso.
	G4UIcmdWith3VectorAndUnit* fTableTransCmd; //trans
	G4UIcmdWithADoubleAndUnit* fTablePivotCmd; //pivot
	G4UIcmdWith3Vector*        fDetCmd; //primary, secondary, SID
	G4UIcmdWithoutParameter*   fRemoveDetCmd;
	G4UIcmdWith3VectorAndUnit* fGlassTransCmd;
	G4UIcmdWith3Vector*        fGlassRotCmd; // axis * angle(in deg)
	G4UIcmdWithoutParameter*   fRemoveGlassCmd;
	G4UIcmdWithABool*          fCurtainCmd;

	G4UIdirectory*         fPhantomDir;
	G4UIcmdWithAnInteger*  fInstallPhantom;
	G4UIcmdWithAnInteger*  fUninstallPhantom;

	G4UIdirectory*        fRecordDir;
	G4UIcmdWithAString*        fReadRecordData;
	G4UIcmdWithAnInteger*      fSetFrameFromRecordData;

	

	G4ThreeVector tableTrans, glassTrans, glassAxis;
	G4double tablePivot, glassTheta;

    std::vector<std::vector<G4float>> recordData;

	G4int maxPeople;
};

#endif
