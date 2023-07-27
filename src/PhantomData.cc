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

#include "PhantomData.hh"

PhantomData::PhantomData(G4String _filePrefix)
{
	// set phantom name
	if(!getenv("VIR_DATA_DIR")){
    	G4Exception("PhantomData::PhantomData","",FatalErrorInArgument,
				        G4String("Set env VIR_DATA_DIR").c_str());
  	}
	G4String filePrefix = string(getenv("VIR_DATA_DIR"))+ _filePrefix;

	G4cout << "================================================================================"<<G4endl;
	G4cout << "\t" << filePrefix << " was implemented in this CODE!!   "<< G4endl;
	G4cout << "================================================================================"<<G4endl;

	MaterialRead(filePrefix + ".material");
	// OrganDoseRead(filePrefix + ".organDose");
	// EffectiveDoseRead(filePrefix + ".effDose");
	RBMBSRead(filePrefix + ".RBMnBS");
	DRFRead(filePrefix + ".DRF");
	ColourRead(filePrefix + ".colour");
}

PhantomData::~PhantomData()
{ 
	for(auto &iter:visAttMap)
		delete iter.second;
	for(auto &iter:materialMap)
		delete iter.second;
}

void PhantomData::MaterialRead(G4String materialFile)
{
	// Read material file (*.material)
	//
	std::ifstream ifpMat;

	ifpMat.open(materialFile);
	if(!ifpMat.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("PhantomData::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + materialFile ).c_str());
	}

	G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String MaterialName;
	G4double density;
	std::map<G4int, G4double> densityMap;
	std::map<G4int, G4String> matNameMap;
	// std::vector<G4int>        materialIndex;
	std::map<G4int, std::vector<std::pair<G4int, G4double>>> materialIndexMap;
	while(!ifpMat.eof())
	{
		ifpMat >> read_data;                   //ex) 'C' RBM
		ifpMat >> MaterialName;                //ex)  C 'RBM'
		ifpMat >> read_data;
		density = std::atof(read_data);        //ex) 1.30
		ifpMat >> read_data;                   //ex) g/cm3
		ifpMat >> read_data;
        if(G4String(read_data).empty()) continue;
		token = std::strtok(read_data,"m");
		G4int matID = std::atoi(token);        //ex) m'10'
        // materialIndex.push_back(matID);
		matNameMap[matID]= MaterialName;
		densityMap[matID] = density*g/cm3;

		for(G4int i=0 ;  ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

			zaid = std::floor(std::atof(read_data)*0.001);
			ifpMat >> read_data;
			fraction = -1.0 * std::atof(read_data);
			materialIndexMap[matID].push_back(std::make_pair(G4int(zaid), fraction));
		}
	}
	ifpMat.close();

	// Construct materials for each organ
	//
    G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
    G4NistManager* nistManager = G4NistManager::Instance();

	// for(G4int i=0;i<(G4int)materialIndex.size();i++){
	for(auto iter:materialIndexMap){
		G4int idx = iter.first;
		if(iter.first==0) continue;
		G4Material* mat = new G4Material(matNameMap[idx], densityMap[idx], G4int(iter.second.size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(auto matIter:iter.second){
			if(matIter.first==1) mat->AddElement(elH, matIter.second);
			else mat->AddElement(nistManager->FindOrBuildElement(matIter.first), matIter.second);
		}
		materialMap[idx]=mat;
		// massMap[idx]=densityMap[idx]*volumeMap[idx];
	}
}

void PhantomData::RBMBSRead(G4String bonefile){
	std::ifstream ifs(bonefile);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("PhantomData::RBMBSRead","",JustWarning,
				G4String("      There is no " + bonefile ).c_str());
		return;
	}
	G4int idx;
	G4double rbm, bs;
	G4double rbmSum(0), bsSum(0);
	while(ifs>>idx>>rbm>>bs){
        if(rbmRatio.find(idx)!=rbmRatio.end()) {
            G4cerr<<idx<<" is duplicated in RBMBS file.."<<G4endl;
            exit(0);
        }
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
		rbmSum += rbm;
		bsSum += bs;
	}
	ifs.close();
	for(auto &iter:rbmRatio)
		iter.second = iter.second / rbmSum;
	for(auto &iter:bsRatio)
		iter.second = iter.second / bsSum;
}

void PhantomData::DRFRead(G4String DRFfile){
	std::vector<G4double> energyBin={0.01 ,0.015 ,0.02 ,0.03 ,0.04 ,0.05 ,0.06 ,0.08
									,0.1 ,0.15 ,0.2 ,0.3 ,0.4 ,0.5 ,0.6 ,0.8
									,1 ,1.5 ,2 ,3 ,4 ,5 ,6 ,8 ,10};
	for(G4double &e:energyBin) e*=MeV;

	std::ifstream ifp;
	ifp.open(DRFfile.c_str());

	if(!ifp.is_open()) {
		G4cerr << DRFfile << " not found!!" << G4endl;
		return;
	}

	G4int ID;
    G4double val;
    while (!ifp.eof()) {
        G4String dump;
        getline(ifp, dump);
        std::stringstream ss(dump); dump.clear();
        ss >> dump;
        if(dump.empty()) continue;
        ID = atoi(dump.c_str());
        // if(rbmDRF.find(ID)!=rbmDRF.end()) {
        //     G4cerr<<ID<<" is duplicated in DRF file.."<<G4endl;
        //     exit(0);
        // }
        // rbmDRF[ID]={};
        // bsDRF[ID]={};
    	for (int j=0; j<energyBin.size(); j++) {
            ss >> val;
			DRF[energyBin[j]*MeV][ID].first = val;
    		// rbmDRF[ID].push_back(DRF);
    	}
    	for (int j=0; j<energyBin.size(); j++) {
            ss >> val;
			DRF[energyBin[j]*MeV][ID].second = val;
    		// bsDRF[ID].push_back(DRF);
    	}
    }
    ifp.close();
}

void PhantomData::ColourRead(G4String colourFile)
{
	// Read colour data file (colour.dat)
	//
	std::ifstream ifpColour(colourFile);

	if(!ifpColour.is_open()) {
		// exception for the case when there is no colour.dat file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("Colour data file was not found ").c_str());
	}

	G4cout << "  Opening colour data file "+colourFile <<G4endl;

	G4int organID;
	G4double red, green, blue, alpha;
	while( ifpColour >> organID >> red >> green >> blue >> alpha )
	{
		visAttMap[organID] = new G4VisAttributes(G4Colour(red, green, blue, alpha));
	}

	ifpColour.close();
}

// void PhantomData::PrintMaterialInfomation()
// {
// 	// Print the overall information for each organ
// 	//
// 	G4cout << G4endl
// 		   << std::setw(9)  << "Organ ID"
// 		   << std::setw(11) << "# of Tet"
// 		   << std::setw(11) << "vol [cm3]"
// 		   << std::setw(11) << "d [g/cm3]"
// 		   << std::setw(11) << "mass [g]"
// 		   << "\t" << "organ/tissue"<< G4endl ;
// 	G4cout << "--------------------------------------------------------------------------------"<<G4endl;

// 	std::map<G4int, G4Material*>::iterator matIter;
// 	G4cout<<std::setiosflags(std::ios::fixed);
// 	G4cout.precision(3);
// 	// for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
// 	// {
// 	// 	G4int idx = matIter->first;

// 	// 	G4cout << std::setw(9)  << idx                         // organ ID
// 	// 		   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
// 	// 		   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
// 	// 		   << std::setw(11) << materialMap[idx]
// 	// 		                       ->GetDensity()/(g/cm3)      // organ density
// 	// 		   << std::setw(11) << massMap[idx]/g              // organ mass
// 	// 		   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
// 	// }

// 	for(auto iter:massMap)
// 	{
// 		G4int idx = iter.first;
// 		G4cout << std::setw(9)  << idx                         // organ ID
// 			   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
// 			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
// 			   << std::setw(11) << 1.089      // organ density
// 			   << std::setw(11) << massMap[idx]/mg              // organ mass
// 			   << "\tskin" << G4endl; // organ name 
// 	}
// }
