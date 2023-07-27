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

#ifndef PhantomData_h
#define PhantomData_h 1

#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>

#include "G4UIExecutive.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4String.hh"
#include "G4Tet.hh"
#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "PhantomAnimator.hh"

class PhantomData
{
public:
	PhantomData(G4String filePrifix);
    virtual ~PhantomData();

	// get methods
	G4Material*   GetMaterial(G4int idx)     { return materialMap[idx];}
	// G4int         GetMaterialIndex(G4int idx){ return materialVector[idx]; }
	G4VisAttributes* GetVisAtt(G4int idx) { return visAttMap[idx]; }
	G4bool IsForVis(){return (visAttMap.size()>0);}
	G4bool GetBoneMassRatio(G4int organID, G4double &rbmR, G4double &bsR)  {
		auto iter = rbmRatio.find(organID);
		if(iter == rbmRatio.end()) return false;
		rbmR = iter->second;
		bsR = bsRatio[organID];

		return true;
	}
	std::vector<G4int> GetBoneIdxVec(){
		std::vector<G4int> boneIdxVec;
		for(auto iter:rbmRatio) boneIdxVec.push_back(iter.first);
		return boneIdxVec;
	}
	// G4double GetBSratio(G4int organID)   { return bsRatio[organID];}
	void GetDRF(G4double energy, G4int organID, G4double &rbmDRF, G4double &bsDRF)
	{
		G4double rbmDRF0, rbmDRF1, bsDRF0, bsDRF1;
		auto iter = DRF.upper_bound(energy);
		if(iter==DRF.begin()) iter = next(iter);
		std::tie(rbmDRF0, bsDRF0) = std::prev(iter)->second[organID];
		std::tie(rbmDRF1, bsDRF1) = iter->second[organID];

		G4double factor = (log10(energy/std::prev(iter)->first))/log10(iter->first/std::prev(iter)->first);

		rbmDRF = exp10(log10( rbmDRF0 ) + factor * log10(rbmDRF1/rbmDRF0));
		bsDRF = exp10(log10( bsDRF0 ) + factor * log10(bsDRF1/bsDRF0));
	}
	// G4double GetRBMDRF(G4int idx, G4int eIdx){ return rbmDRF[idx][eIdx];}
	// G4double GetBSDRF (G4int idx, G4int eIdx){ return bsDRF[idx][eIdx];}


private:
	// private methods
    void MaterialRead(G4String);
	void RBMBSRead(G4String);
	void DRFRead(G4String);
	void ColourRead(G4String);

	// G4String filePrefix;

	// G4ThreeVector boundingBox_Min;
	// G4ThreeVector boundingBox_Max;

	// std::vector<G4int>         materialVector;
	// std::map<G4int, G4Colour>  colourMap;
	std::map<G4int, G4VisAttributes*>  visAttMap;
	std::map<G4int, G4double>  rbmRatio;
	std::map<G4int, G4double>  bsRatio;

	std::map<G4double, std::map<G4int, std::pair<G4double, G4double>>> DRF; //rbm, bs
	// std::map<G4int, std::vector<G4double>> rbmDRF;
	// std::map<G4int, std::vector<G4double>> bsDRF;

	std::map<G4int, G4Material*>                             materialMap;

	// std::map<G4int, G4double>           effectiveDoseInfo;

};

#endif