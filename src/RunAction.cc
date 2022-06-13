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
// TETRunAction.cc
// \file   MRCP_GEANT4/External/src/TETRunAction.cc
// \author Haegin Han
//

#include "G4Timer.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include <iostream>
#include "RunAction.hh"

RunAction::RunAction(TETModelImport* _tetData_doctor,
G4String _output, G4Timer* _init, WEBServerConnect* _serverConnect)
:tetData_doctor(_tetData_doctor),
 outputFile(_output), initTimer(_init), serverConnect(_serverConnect)
{
	if(!isMaster) return;
	runTimer = new G4Timer;
	
	if (tetData_doctor != nullptr) {
		ofs_doctor.open(outputFile);
		if (tetData_doctor->DoseWasOrganized()) {
			doctor_massMap = tetData_doctor->GetDoseMassMap();
			for(auto itr: doctor_massMap) doctor_nameMap[itr.first] = tetData_doctor->GetDoseName(itr.first);
		}
		else {
			doctor_massMap = tetData_doctor->GetMassMap();
			for (auto itr: doctor_massMap) doctor_nameMap[itr.first] = tetData_doctor->GetMaterial(itr.first)->GetName();
		}
		doctor_nameMap[-4] = "RBM(DRF)"; doctor_nameMap[-3] = "BS(DRF)";
		doctor_nameMap[-2] = "RBM"     ; doctor_nameMap[-1] = "BS"     ;
	}
	
}

RunAction::~RunAction()
{
	ofs_doctor.close();
}

G4Run* RunAction::GenerateRun()
{
	// generate run
	fRun = new Run(tetData_doctor);
	return fRun;
}


void RunAction::BeginOfRunAction(const G4Run* aRun)
{
	// print the progress at the interval of 10%
	numOfEvent=aRun->GetNumberOfEventToBeProcessed();
	G4RunManager::GetRunManager()->SetPrintProgress(G4int(numOfEvent*0.05));

	if (isMaster) {
		initTimer->Stop();
		runTimer->Start();
	}
}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
	if(!isMaster) return;
	runTimer->Stop();

	runID = aRun->GetRunID();
		
	if (tetData_doctor != nullptr) {
		cout << ">> Doctor Dose" << endl;
		HITSMAP doctor_edepMap = *fRun->GetdoctorEdepMap();
		SetDoses(doctor_edepMap, doctor_massMap, doctor_doses);
		if(tetData_doctor->DoseWasOrganized())
			SetEffectiveDose(doctor_doses, doctor_eff, doctor_eff_DRF);
		PrintResultExternal(ofs_doctor, tetData_doctor, doctor_massMap, doctor_nameMap, doctor_doses, doctor_eff_DRF, doctor_eff);	

		if (serverConnect != nullptr) {
			cout << "Send dose results to web server..." << endl;
			if ( (doctor_nameMap.size()+(3-2))*2 != serverConnect->GetColDoseResults().size() ) {
				cerr << "Check the dose file again. Organ names and Web-server table comlumn names are not matched" << endl; 
				exit(1);
			}
			else 
				serverConnect->SendDoseResultsToWebServer(doctor_nameMap, doctor_doses, doctor_eff_DRF);		
		}
	}

	initTimer->Start();
}

void RunAction::SetDoses(HITSMAP edepMap, map<G4int, G4double> massMap, map<G4int, pair<G4double, G4double>> &doses)
{
	doses.clear();

	// RBM and BS
	for (G4int i=-4; i<0; i++) {
		G4double meanDose    = edepMap[i].first  / numOfEvent;
		G4double squareDose  = edepMap[i].second;
		G4double variance    = ((squareDose/numOfEvent) - (meanDose*meanDose)) / numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
		doses[i] = make_pair(meanDose, relativeE);
	}
	doses[-4].first *= joule/kg;
	doses[-3].first *= joule/kg;

	// Other organs
	for(auto itr : massMap) {
		if (itr.first<0) continue; // continue for RBM and BS
		G4double meanDose    = edepMap[itr.first].first / numOfEvent;
		G4double squareDoese = edepMap[itr.first].second;
		G4double variance    = ((squareDoese/numOfEvent) - (meanDose*meanDose))/numOfEvent;
		G4double relativeE   = sqrt(variance)/meanDose;
        doses[itr.first] = std::make_pair(meanDose/itr.second, relativeE);
	}
}

void RunAction::SetEffectiveDose(map<G4int, pair<G4double, G4double>> &doses, 
								 pair<G4double, G4double> &effective,
								 pair<G4double, G4double> &effective_DRF)
{
	// ICRP-103 Table B.2
	//
	// Group 1: w_T = 0.12
	G4int group1_RBM = -2; // RBM(-2,-4)
	G4int group1_RBM_DRF = -4;
	std::vector<G4int> group1 = {2, 4, 5, 7}; // Colon(2), Lungs(4), Stomach(5), Breasts(7)
	G4int ET1(21), ET2(22); // ET(21,22)
	std::vector<G4int> remainder = {20, 25, 26, 27, 28, 29, 30, 31, 32, 33, 35, 36};
	// adrenals(20), gall bladder(25), heart(26), kidneys(27), lymph nodes(28), muscle(29), oral mucosa(30)
	// pancrease(31), prostate/uterus(32), SI(33), spleen(35), thymus(36)

	// Group 2: w_T = 0.08
	std::vector<G4int> group2 = {8}; // Testes/Ovaries(8)

	// Group 3: w_T = 0.04
	std::vector<G4int> group3 = {9, 11, 13, 14}; // UB(9), Oesophagus(11), Liver(13), Thyroid(14)

	// Group 4: w_T = 0.01
	std::vector<G4int> group4 = {16, 17, 19}; // Brain(16), SalivaryGlands(17), Skin(19, Target)
	G4int group4_BS = -1; // Bone Surface(-1,-3)
	G4int group4_BS_DRF = -3;

	std::vector<std::pair<G4double, G4double>> effDoseComp;
	for(G4int idx:group1)	effDoseComp.push_back(doses[idx]);
	effDoseComp.push_back(doses[group1_RBM]);
	for(G4int idx:group2)	effDoseComp.push_back(doses[idx]);
	for(G4int idx:group3)	effDoseComp.push_back(doses[idx]);
	for(G4int idx:group4)	effDoseComp.push_back(doses[idx]);
	effDoseComp.push_back(doses[group4_BS]);
	for(G4int idx:remainder)	effDoseComp.push_back(doses[idx]);
	effDoseComp.push_back(doses[ET1]); effDoseComp.push_back(doses[ET2]);

	std::vector<G4double> ratios;
    for(size_t i=0;i<group1.size()+1;i++) ratios.push_back(0.12);
    for(size_t i=0;i<group2.size();i++) ratios.push_back(0.08);
    for(size_t i=0;i<group3.size();i++) ratios.push_back(0.04);
    for(size_t i=0;i<group4.size()+1;i++) ratios.push_back(0.01);
    for(size_t i=0;i<remainder.size();i++) ratios.push_back(0.12/(G4double)(remainder.size()+1));
    ratios.push_back(0.12/(G4double)(remainder.size()+1)*0.001);
    ratios.push_back(0.12/(G4double)(remainder.size()+1)*0.999);

    std::vector<std::pair<G4double, G4double>> effDoseComp_DRF(effDoseComp);
	effDoseComp_DRF[0] = doses[group1_RBM_DRF];
	effDoseComp_DRF[13] = doses[group4_BS_DRF];

	G4double sum(0.);
	for(auto r:ratios) sum += r;
	G4cout<<"Sum of the ratios -->"<<sum<<G4endl;

	effective = PropagateError(effDoseComp, ratios);
	effective_DRF = PropagateError(effDoseComp_DRF, ratios);
}

pair<G4double, G4double> RunAction::PropagateError(vector<pair<G4double, G4double>> &doseVec, vector<G4double> &ratio)
{
	typedef pair<G4double, G4double> VALUE;

	// normalize the ratio
	if (doseVec.size() != ratio.size()) {
		G4Exception("TETRunAction::PropagateError","",JustWarning,
				G4String("      uniform ratio was applied " ).c_str());
		std::vector<G4double> uniform(doseVec.size(),1./(G4int)doseVec.size());
		ratio = uniform;
	}
	else {
		G4double sum(0.);
		for (size_t i=0; i< ratio.size(); i++) sum += ratio[i];
		for (auto &r:ratio) r /= sum;
	}

	// set mean dose
	G4double value(0.);
	for (size_t i=0; i<doseVec.size(); i++) value += doseVec[i].first * ratio[i];
	
	// set error
	G4double error(0.);
	for (size_t i=0; i<doseVec.size(); i++) {
		if(isnan(doseVec[i].second)) continue;
		error += pow(doseVec[i].first * doseVec[i].second * ratio[i], 2);
	}
	error = sqrt(error);
	error /= value;
	
	return VALUE(value, error);
}


void RunAction::PrintResultExternal(ostream &out, 
									TETModelImport* tetData,
									map<G4int, G4double> massMap, 
									map<G4int, G4String> nameMap, 
									map<G4int, pair<G4double,G4double>> doses,
									pair<G4double, G4double> effective_DRF,
									pair<G4double, G4double> effective)
{
	// Print run result
	//

	out << G4endl
	    << "=======================================================================" << G4endl
	    << " Run #" << runID << " / Number of event processed : "<< numOfEvent     << G4endl
	    << "=======================================================================" << G4endl
		<< " Init time: " << initTimer->GetRealElapsed() << " s / Run time: "<< runTimer->GetRealElapsed()<<" s"<< G4endl
	    << "=======================================================================" << G4endl
	    << setw(27) << "organ ID| "
		<< setw(15) << "Organ Mass (g)"
        << setw(15) << "Dose (pGy/nps)"
		<< setw(15) << "Relative Error" << G4endl;

	out.precision(3);

	for(G4int i=-4;i<0;i++){
		out << setw(25) << nameMap[i] << "| ";
		out << setw(30) << scientific << doses[i].first/(joule/kg) << setw(15) << fixed << doses[i].second << G4endl;
	}

	for(auto itr : massMap){
		if(tetData->DoseWasOrganized()||itr.first<0) out << setw(25) << nameMap[itr.first]<< "| ";
		else                            out << setw(25) << tetData->GetMaterial(itr.first)->GetName()<< "| ";
		out	<< setw(15) << fixed      << itr.second/g;
        out	<< setw(15) << scientific << doses[itr.first].first/(joule/kg)*1e12;
		out	<< setw(15) << fixed      << doses[itr.first].second << G4endl;
	}

	//effective dose
	out << setw(25) << "eff. dose (DRF)" << "| ";
	out	<< setw(15) << " "                ;
    out	<< setw(15) << scientific << effective_DRF.first/(joule/kg)*1e12;
	out	<< setw(15) << fixed      << effective_DRF.second << G4endl;

	out << setw(25) << "eff. dose" << "| ";
	out	<< setw(15) << " "                ;
    out	<< setw(15) << scientific << effective.first/(joule/kg)*1e12;
	out	<< setw(15) << fixed      << effective.second << G4endl;

	out << "=======================================================================" << G4endl << G4endl;
}