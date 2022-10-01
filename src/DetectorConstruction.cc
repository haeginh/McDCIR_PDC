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

#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "TETParameterisation.hh"
#include "Eigen/Core"

using namespace std;

DetectorConstruction::DetectorConstruction()
:worldLogical(0) ,worldPhysical(0), isoZ(113.5*cm), table_ref_pos(-229.35*mm, 447.5*mm, -184.285*mm), table_size(500*mm, 3115*mm, 1.43*mm),//table_default_trans(1.12*m, -0.3*m, 0.785*m),
curtain_size(0.5*mm, 500*mm, 400*mm),table_rotation_center(0, 0, 0), head_margin(10.*cm), curtain_margin(10.*cm)

{
	SetTableSize(G4ThreeVector(55*cm, 290*cm, 1.43*mm));
	SetTableRefPos(G4ThreeVector(-327, 398, -175)*mm);
	// table_rotation_center = table_default_trans;
	// Operating Table (w/ patient, curtain)
	// table_ocr = G4ThreeVector(-280,-960,-10); // Initial position
	// table_trans = G4ThreeVector(-80*mm, 20*mm -10*mm);     // Operating position
	// table_pivot_angle = 0 * deg;              // Z axis rotation
	// table_rotation_center = G4ThreeVector(1200*mm, 0*mm, 820*mm);

	// C-arm det
	// carm_primary   = 20 * deg;   // +LAO, -RAO
	// carm_secondary = 20 * deg;   // +CAU, -CRA

	messenger = new DetectorMessenger(this);
}

DetectorConstruction::~DetectorConstruction()
{
	delete messenger;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
	// ((ParallelGlass*)GetParallelWorld(0))->SetIsoCenter(carm_isocenter);
	SetupWorldGeometry();
	ConstructOperatingTable();
	// ConstructPbGlass();
	// ConstructCarmDet();
	return worldPhysical;
}

void DetectorConstruction::SetupWorldGeometry()
{
	// Define the world box (size: 10*10*5 m3)
    G4double worldHalfX = 5. * m;
    G4double worldHalfY = 5. * m;
    G4double worldHalfZ = 2.5 * m;

	G4VSolid* worldSolid = new G4Box("worldSolid", worldHalfX, worldHalfY, worldHalfZ);
	worldLogical = new G4LogicalVolume(worldSolid,G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"),"worldLogical");
	worldPhysical = new G4PVPlacement(0,G4ThreeVector(), worldLogical,"worldPhysical", 0, false,0,false);
	worldLogical->SetVisAttributes(G4VisAttributes::GetInvisible());
	// G4VisAttributes* va_World = new G4VisAttributes(G4Colour(1.0, 1.0, 1.0));
	// va_World->SetForceWireframe(true);
	// worldLogical->SetVisAttributes(va_World);

	// Define floor
	G4Material* concrete = G4NistManager::Instance()->FindOrBuildMaterial("G4_CONCRETE");
	G4Box* floor = new G4Box("floor", worldHalfX, worldHalfY, (worldHalfZ-isoZ)*0.5);
	G4LogicalVolume* lv_floor = new G4LogicalVolume(floor, concrete, "lv_floor");
	new G4PVPlacement(0, G4ThreeVector(0,0,(-worldHalfZ)+(worldHalfZ-isoZ)*0.5), lv_floor, "pv_floor", worldLogical,0,0);

	G4VisAttributes* vis_floor  = new G4VisAttributes(G4Colour(0.8,0.8,0.8,0.5));
	vis_floor->SetForceAuxEdgeVisible();
	lv_floor->SetVisAttributes(vis_floor);
}


#include "G4GeometryManager.hh"
void DetectorConstruction::ConstructOperatingTable()
{
	auto rot = new G4RotationMatrix;

	// Operating table
	G4Box* table = new G4Box("sol_table", table_size.x()*0.5, table_size.y()*0.5, table_size.z()*0.5);
	G4Material* table_mat = new G4Material("tableMat", 5.5382*g/cm3, 2);
	table_mat->AddElement(G4NistManager::Instance()->FindOrBuildElement(6), 0.5127);
	table_mat->AddElement(G4NistManager::Instance()->FindOrBuildElement(13), 0.4873);
	G4LogicalVolume* lv_table = new G4LogicalVolume(table, table_mat, "lv_table");
	table_center0 = G4ThreeVector(table_ref_pos.x() + table_size.x()*0.5, table_ref_pos.y()-table_size.y()*0.5, table_ref_pos.z()-table_size.z()*0.5);
	pv_table = new G4PVPlacement(rot, table_center0, lv_table, "pv_table", worldLogical, false, 1);

	lv_table->SetVisAttributes( new G4VisAttributes(G4Colour(1.,1.,0.)));

	// phantom box
	G4LogicalVolume* lv_phantomBox = ConstructPatient(patient);
	lv_phantomBox->SetVisAttributes(G4VisAttributes::GetInvisible());

	phantom_center0 = table_center0 + 
	                        G4ThreeVector(0., table_size.y()*0.5-((G4Box*) lv_phantomBox->GetSolid())->GetYHalfLength(),
							                  table_size.z()*0.5 + ((G4Box*) lv_phantomBox->GetSolid())->GetZHalfLength());
	pv_phantom = new G4PVPlacement(rot, phantom_center0,lv_phantomBox, "phantom box", worldLogical, false, 0);

	// Pb Curtain - lead equivalence 0.5 mm Pb
	// G4ThreeVector curtainHalfSize = curtain_size*0.5;
	// G4Box* curtain = new G4Box("sol_curtain", curtainHalfSize.x(), curtainHalfSize.y(), curtainHalfSize.z());
	// curtain_center0.setX(table_center0.x() - table_size.x()*0.5 - curtainHalfSize.x());
	// curtain_center0.setY(table_center0.y() + table_size.y()*0.5 - curtain_margin - curtainHalfSize.y());
	// curtain_center0.setZ(table_center0.z() -table_size.z()*0.5 - curtainHalfSize.z());
	// auto lv_curtain = new G4LogicalVolume(curtain, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_curtain");
	// pv_curtain = new G4PVPlacement(rot, curtain_center0, lv_curtain, "pv_curtain", worldLogical, false, 1);

	// lv_curtain->SetVisAttributes( new G4VisAttributes(G4Colour(0.,0.,1.,0.8)) );
}

G4LogicalVolume* DetectorConstruction::ConstructPatient(G4String _patient)
{
	Eigen::MatrixXd V;
	G4cout<<"Read "+_patient+".node..."<<flush;
	G4String dir = std::getenv("DCIR_PHANTOM_DIR");
	ifstream ifs(dir+_patient+".node");
	if(!ifs.is_open()){
		cout<<"there is no "+dir+_patient+".node!"<<endl;
		exit(100);
	}
	int rowV, tmp;
	ifs>>rowV>>tmp>>tmp>>tmp;
	V.resize(rowV, 3);
	for(int i=0;i<rowV;i++)
	{
		G4double x, y, z;
		ifs>>tmp>>x>>y>>z;
		V.row(i) = Eigen::RowVector3d(x,y,z)*cm;
	}
	ifs.close();
	G4cout<<"done"<<G4endl;

	//move the bbox center to the origin
	Eigen::RowVector3d max= V.colwise().maxCoeff();
	Eigen::RowVector3d min= V.colwise().minCoeff();
	G4cout<<-(max+min)*0.5/cm <<endl;
	V = V.rowwise()-(max+min)*0.5;

	//set the phantom box with no margin
	Eigen::RowVector3d hlafSize = (max-min)*0.5;
	G4VSolid* phantomBox = new G4Box("patientBox", hlafSize.x(), hlafSize.y(), hlafSize.z());
	G4LogicalVolume* lv_phantomBox = new G4LogicalVolume(phantomBox, G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR"), "patientBox");
	
	//read materials
	G4cout<<"Read "+_patient+".material..."<<flush;
	ifstream ifsMat(dir+_patient+".material");
	if(!ifsMat.is_open()){
		cout<<"there is no "+dir+_patient+".material!"<<endl;
		exit(100);
	}
	matMap.clear();
	G4Element *elH = new G4Element("TS_H_of_Water", "H", 1., 1.01*g/mole);
	while(!ifsMat.eof())
	{
		G4String line, first;
		getline(ifsMat, line);
		stringstream ss(line);
		ss>>first;
		if(first=="C")
		{
			G4String name, tmp1, idStr;
			G4double density;
			ss >> name >> density >> tmp1;
			ifsMat >> idStr;
			map<int, double> composition;
			while(1)
			{
				G4String tmp2;
				ifsMat>>tmp2;
				if(tmp2.empty()||tmp2=="C") break;
				G4double fraction;
				ifsMat>>fraction;
				composition[atoi(tmp2.c_str())/1000] = -fraction;
			}
			G4Material* mat = new G4Material(name, density*g/cm3, composition.size());
			for(auto iter:composition)
			{
				if(iter.first==1) mat->AddElement(elH, iter.second);
				else mat->AddElement(G4NistManager::Instance()->FindOrBuildElement(iter.first), iter.second);
			}
			G4int idx = atoi(idStr.substr(1, idStr.length()-1).c_str());
			if(matMap.find(idx)!= matMap.end())
			{
				G4cerr<<"There is duplicated material #"<<idx<<G4endl;
				exit(100);
			}
			matMap[atoi(idStr.substr(1, idStr.length()-1).c_str())] = mat;
		}
	}
	ifsMat.close();
	G4cout<<"done"<<G4endl;

	//construc tets
	G4cout<<"Read "+_patient+".ele..."<<flush;
	ifstream ifsEle(dir+_patient+".ele");
	if(!ifsEle.is_open()){
		cout<<"there is no " + dir +_patient+".ele!"<<endl;
		exit(100);
	}
	tetVec.clear();
	int numT;
	ifsEle>>numT>>tmp>>tmp;
	idArray = Eigen::ArrayXi::Zero(numT);
	// map<int, double> volMap;
	// map<int, int> numMap;
	// G4Material* tissue = G4NistManager::Instance()->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRP");
	for(int i=0;i<numT;i++)
	{
		int a, b, c, d, id;
		ifsEle>>tmp>>a>>b>>c>>d>>id;
		// if(id!=12200) continue;
		tetVec.push_back(new G4Tet("tet", G4ThreeVector(V(a,0),V(a,1),V(a,2)),
		                                 G4ThreeVector(V(b,0),V(b,1),V(b,2)),
										 G4ThreeVector(V(c,0),V(c,1),V(c,2)),
										 G4ThreeVector(V(d,0),V(d,1),V(d,2))));
		// volMap[id] += tet->GetCubicVolume();
		// numMap[id]++;
		idArray(i) = id;
		if(matMap[id]==nullptr) {
			G4cerr<<"Invalid material #"<<id<<G4endl;
			exit(100);
		}
		// G4LogicalVolume* lv_tet = new G4LogicalVolume(tet, matMap[id], "tet");
		// G4LogicalVolume* lv_tet = new G4LogicalVolume(tet, tissue, "tet");
		// new G4PVPlacement(0, G4ThreeVector(), lv_tet, "tet", lv_phantomBox, false, 1);
		// G4VisAttributes* vis = new G4VisAttributes(); vis->SetForceCloud();
		// if(i>2000) lv_tet->SetVisAttributes(G4VisAttributes::GetInvisible());
		// lv_tet->SetVisAttributes(G4VisAttributes(G4Color(1,1,1,0.1)));
	}
	ifsEle.close();
	G4cout<<"done"<<G4endl;

	new G4PVParameterised("phantom_paramPhy", new G4LogicalVolume(tetVec[0], matMap.begin()->second, "tetLog"),
	                     lv_phantomBox, kUndefined, tetVec.size(), new TETParameterisation(this));
	//print summary
	// G4cout << G4endl
	// 	   << std::setw(9)  << "Organ ID"
	// 	   << std::setw(11) << "# of Tet"
	// 	   << std::setw(11) << "vol [cm3]"
	// 	   << std::setw(11) << "d [g/cm3]"
	// 	   << std::setw(11) << "mass [g]"
	// 	   << "\t" << "organ/tissue"<< G4endl ;
	// G4cout << "--------------------------------------------------------------------------------"<<G4endl;

	// std::map<G4int, G4Material*>::iterator matIter;
	// G4cout<<std::setiosflags(std::ios::fixed);
	// G4cout.precision(3);
	// for(matIter=matMap.begin(); matIter!=matMap.end();matIter++)
	// {
	// 	G4int idx = matIter->first;

	// 	G4cout << std::setw(9)  << idx                         // organ ID
	// 		   << std::setw(11) << numMap[idx]              // # of tetrahedrons
	// 		   << std::setw(11) << volMap[idx]/cm3          // organ volume
	// 		   << std::setw(11) << matMap[idx]
	// 		                       ->GetDensity()/(g/cm3)      // organ density
	// 		   << std::setw(11) << volMap[idx] * matMap[idx]->GetDensity()/g              // organ mass
	// 		   << "\t"<<matMap[idx]->GetName() << G4endl; // organ name
	// }
	return lv_phantomBox;
}


// void DetectorConstruction::ConstructPbGlass()
// {
// 	// Pb Glass - 40 x 50 cm tiltable lead acrylic shield, lead equivalence 0.5 mm Pb
// 	G4Box* glass = new G4Box("sol_glass", 200*mm, 250*mm, 0.25*mm);
// 	G4LogicalVolume* lv_glass = new G4LogicalVolume(glass, G4NistManager::Instance()->FindOrBuildMaterial("G4_Pb"), "lv_glass");
// 	lv_glass->SetVisAttributes( new G4VisAttributes(G4Colour(0.,1.,0.,0.8)) );
// 	pv_glass  = new G4PVPlacement(new G4RotationMatrix, G4ThreeVector(), lv_glass, "pv_glass", worldLogical, false, 0);
// }






