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

#include "TETModelImport.hh"

TETModelImport::TETModelImport(G4String _phantomName)
:doseOrganized(false)
{
	// set phantom name
	phantomName = _phantomName;

	G4cout << "================================================================================"<<G4endl;
	G4cout << "\t" << phantomName << " was implemented in this CODE!!   "<< G4endl;
	G4cout << "================================================================================"<<G4endl;

	G4String materialFile =  phantomName + ".material";
	G4String doseFile     =  phantomName + ".dose";
	G4String boneFile     =  phantomName + ".RBMnBS";
	G4String drfFile      =  phantomName + ".DRF";

	animator = new PhantomAnimator(phantomName);
	animator->CalibrateTo("S_Moon");
	UpdateBBox();
}

TETModelImport::~TETModelImport()
{ 
	delete animator;
}

void TETModelImport::Deform(RotationList vQ, Vector3d root)
{
	animator->Animate(vQ, root);
	UpdateBBox();

	MatrixXi T = animator->GetT();
	MatrixXd U = animator->GetU();
	for (G4int i=0; i<T.rows(); i++)
		tetVector[i]->SetVertices(RowToG4Vec(U.row(T(i,0))),RowToG4Vec(U.row(T(i,1))),RowToG4Vec(U.row(T(i,2))),RowToG4Vec(U.row(T(i,3))));
}

// void TETModelImport::ReadFrameRecord()
// {
// 	ifstream ifs("../phantoms/record.txt");
// 	if(!ifs.is_open()) {
// 		cerr << "record file is not opened" << endl;
// 		exit(1);
// 	}

// 	int frameNo(0);

// 	string dump;
// 	while (getline(ifs, dump))
// 	{
// 		stringstream ss(dump);
// 		Vector3d position;
// 		ss >> position(0) >> position(1) >> position(2);
// 		rootPosition.push_back(position);

// 		vector<Quaterniond,aligned_allocator<Quaterniond>> qVec;
// 		qVec.clear();
// 		Quaterniond q;
// 		for (int i=0; i<22; i++) {
// 			ss >> q.w() >> q.x() >> q.y() >> q.z();
// 			qVec.push_back(q);
// 		}
// 		boneQuat.push_back(qVec);

// 		frameNo++;
// 	}

// 	ifs.close();
// }


void TETModelImport::DoseRead(G4String doseFile){
	//read dose file : PLEASE be careful not to include dose ID 0
	std::ifstream ifs(doseFile);
	if(!ifs.is_open()) return;
	doseOrganized = true;

	G4String aLine;
	while(!ifs.eof()){
		getline(ifs, aLine);
		if(aLine.empty()) break;

		std::stringstream ss(aLine);
		G4int doseID; ss>>doseID;
		G4String name; ss>>name; doseName[doseID] = name;
		G4int organID;
		while(ss>>organID){
			if(organ2dose.find(organID)==organ2dose.end()) organ2dose[organID] = {doseID};
			else	                                       organ2dose[organID].push_back(doseID);
		}
	}
	ifs.close();

}

void TETModelImport::ConstructTet()
{
	MatrixXi T = animator->GetT();
	MatrixXd U = animator->GetU();
	G4ThreeVector center = (boundingBox_Max + boundingBox_Min)*0.5;
	for(G4int i=0; i<T.rows(); i++)
	{
		materialVector.push_back(T(i,4));

		// save the element (tetrahedron) data as the form of std::vector<G4Tet*>
		tetVector.push_back(new G4Tet("Tet_Solid",
							   		  RowToG4Vec(U.row(T(i,0)))-center,
									  RowToG4Vec(U.row(T(i,1)))-center,
									  RowToG4Vec(U.row(T(i,2)))-center,
									  RowToG4Vec(U.row(T(i,3)))-center));

		// calculate the total volume and the number of tetrahedrons for each organ
		std::map<G4int, G4double>::iterator FindIter = volumeMap.find(materialVector[i]);

		if(FindIter!=volumeMap.end()){
			FindIter->second += tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]]++;
		}
		else {
			volumeMap[materialVector[i]] = tetVector[i]->GetCubicVolume();
			numTetMap[materialVector[i]] = 1;
		}
	}
}

void TETModelImport::MaterialRead(G4String materialFile)
{
	// Read material file (*.material)
	//
	std::ifstream ifpMat;

	ifpMat.open(materialFile);
	if(!ifpMat.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String("      There is no " + materialFile ).c_str());
	}

	G4cout << "  Opening material file '" << materialFile << "'" <<G4endl;

	char read_data[50];
	char* token;
	G4double zaid;
	G4double fraction;
	G4String MaterialName;
	G4double density;

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
        materialIndex.push_back(matID);
		organNameMap[matID]= MaterialName;
		densityMap[matID] = density*g/cm3;

		for(G4int i=0 ;  ; i++)
		{
			ifpMat >> read_data;
			if(std::strcmp(read_data, "C")==0 || ifpMat.eof()) break;

			zaid = (G4int)(std::atoi(read_data)/1000);
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

	for(G4int i=0;i<(G4int)materialIndex.size();i++){
		G4int idx = materialIndex[i];
		if(idx==0) continue;
		G4Material* mat = new G4Material(organNameMap[idx], densityMap[idx], G4int(materialIndexMap[idx].size()), kStateSolid, NTP_Temperature, STP_Pressure);
		for(G4int j=0;j<G4int(materialIndexMap[idx].size());j++){
			if(materialIndexMap[idx][j].first==1) mat->AddElement(elH, materialIndexMap[idx][j].second);
			else mat->AddElement(nistManager->FindOrBuildElement(materialIndexMap[idx][j].first), materialIndexMap[idx][j].second);
		}
		materialMap[idx]=mat;
		massMap[idx]=densityMap[idx]*volumeMap[idx];
	}

	if(DoseWasOrganized()){
		for(auto dm:doseName){
			doseMassMap[dm.first] = 0;
		}
		for(auto od:organ2dose){
			for(auto doseID:od.second)
				doseMassMap[doseID] += massMap[od.first];
		}
	}
}

void TETModelImport::RBMBSRead(G4String bonefile){
	std::ifstream ifs(bonefile);
	if(!ifs.is_open()) {
		// exception for the case when there is no *.material file
		G4Exception("TETModelImport::RBMBSRead","",JustWarning,
				G4String("      There is no " + bonefile ).c_str());
		return;
	}
	G4int idx;
	G4double rbm, bs;
	while(ifs>>idx>>rbm>>bs){
        if(rbmRatio.find(idx)!=rbmRatio.end()) {
            G4cerr<<idx<<" is duplicated in RBMBS file.."<<G4endl;
            exit(0);
        }
		rbmRatio[idx]=rbm;
		bsRatio[idx]=bs;
	}
}

void TETModelImport::DRFRead(G4String DRFfile){
	std::ifstream ifp;
	ifp.open(DRFfile.c_str());

	if(!ifp.is_open()) {
		G4cerr << DRFfile << " not found!!" << G4endl;
		return;
	}

	G4int ID;
    G4double DRF;
    while (!ifp.eof()) {
        G4String dump;
        getline(ifp, dump);
        std::stringstream ss(dump); dump.clear();
        ss >> dump;
        if(dump.empty()) continue;
        ID = atoi(dump.c_str());
        if(rbmDRF.find(ID)!=rbmDRF.end()) {
            G4cerr<<ID<<" is duplicated in DRF file.."<<G4endl;
            exit(0);
        }
        rbmDRF[ID]={};
        bsDRF[ID]={};
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		rbmDRF[ID].push_back(DRF);
    	}
    	for (int j=0; j<25; j++) {
            ss >> DRF;
    		bsDRF[ID].push_back(DRF);
    	}
    }
    ifp.close();
}

void TETModelImport::ColourRead(G4String mtlFile)
{
	// Read colour data file (colour.dat)
	//
	std::ifstream ifpColour;

	ifpColour.open(mtlFile);
	if(!ifpColour.is_open()) {
		// exception for the case when there is no colour.dat file
		G4Exception("TETModelImport::DataRead","",FatalErrorInArgument,
				G4String(mtlFile + " file was not found ").c_str());
	}

	G4cout << "  Opening colour data file "+mtlFile <<G4endl;

	G4int organID;
	G4ThreeVector rgb;
	G4String dump;
	while( ifpColour >> dump ){
		if(dump=="newmtl"){
			ifpColour>>dump;
			organID=GetID(dump);
		}
		else if(dump=="Kd"){
			ifpColour>>rgb;
			colourMap[organID] = G4Colour(rgb.getX(), rgb.getY(), rgb.getZ(), 0.);
		}
	}
	ifpColour.close();
	colourMap[160].SetAlpha(0.1);
	colourMap[140].SetAlpha(0.1);
	colourMap[141].SetAlpha(0.1);
	colourMap[142].SetAlpha(0.1);
}

void TETModelImport::PrintMaterialInfomation()
{
	// Print the overall information for each organ
	//
	G4cout << G4endl
		   << std::setw(9)  << "Organ ID"
		   << std::setw(11) << "# of Tet"
		   << std::setw(11) << "vol [cm3]"
		   << std::setw(11) << "d [g/cm3]"
		   << std::setw(11) << "mass [g]"
		   << "\t" << "organ/tissue"<< G4endl ;
	G4cout << "--------------------------------------------------------------------------------"<<G4endl;

	std::map<G4int, G4Material*>::iterator matIter;
	G4cout<<std::setiosflags(std::ios::fixed);
	G4cout.precision(3);
	for(matIter=materialMap.begin(); matIter!=materialMap.end();matIter++)
	{
		G4int idx = matIter->first;

		G4cout << std::setw(9)  << idx                         // organ ID
			   << std::setw(11) << numTetMap[idx]              // # of tetrahedrons
			   << std::setw(11) << volumeMap[idx]/cm3          // organ volume
			   << std::setw(11) << materialMap[idx]
			                       ->GetDensity()/(g/cm3)      // organ density
			   << std::setw(11) << massMap[idx]/g              // organ mass
			   << "\t"<<materialMap[idx]->GetName() << G4endl; // organ name
	}
}
