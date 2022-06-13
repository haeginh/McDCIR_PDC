// Geant4, RunManager
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"
// Geant4, Geometry
#include "DetectorConstruction.hh"
#include "ParallelGlass.hh"
#include "ParallelPhantom.hh"
// Genat4, Physics
#include "G4PhysListFactory.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
// Geant4, Initialization
#include "ActionInitialization.hh"
// Geant4, etc
#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "Randomize.hh"
// Phantom
#include "TETModelImport.hh"
#include "PhantomAnimator.hh"
// SQL
#include "WEBServerConnect.hh"

// History
//
// Rev 01. 2022.05.18 - Reconstruct the code

void PrintUsage()
{
	G4cout << "==================================================================" << G4endl;
	G4cout << "Monte Carlo Dose Calculation for Interventional Radiology (McDCIR)" << G4endl;
	G4cout << "Dose Calculation Mode: Post-procedure Dose Calculation (PDC)" << G4endl;
	G4cout << "Contributor: Haegin Han, Sungho Moon, Gahee Son" << G4endl << G4endl;

	G4cout << "./calculator [option 1] [name 1] [option 2] [name 2] ..." << G4endl;
	G4cout << "[options]" << G4endl;
	G4cout << "-m [G4 macro]         : Geant4 run macro file" << G4endl;
	G4cout << "-o [Output]           : Dose results file name (.epd: Hp10, .pat: patient dose, .doc: doctor dose)" << G4endl;
	G4cout << "-p [Patient phantom]  : Patient phantom name" << G4endl;
	G4cout << "-d [Doctor phantom]   : Doctor phantom name (If doctor phantom is not set up, doctor will not installed)" << G4endl;
	G4cout << "-v                    : Qt Visualization (Don't use it with '-m' option)" << G4endl;
	G4cout << "-sql                  : Run PDC system (This option will be default mode, Don't use it with '-m' option)" << G4endl;
	G4cout << "-id [Server id]       : Read server # parameter of 'PDC.server' and set them in PDC system. (Use it with '-sql' option)" << G4endl;

	G4cout << "example: ./calculator -m run.mac -o output -p RANDO_M -d M_H175M83800" << G4endl;
	G4cout << "example: ./calculator -p RANDO_M -d M_H175M83800 -v" << G4endl;
	G4cout << "example: ./calculator -p RANDO_M -d M_H175M83800 -id 1 -sql" << G4endl;
}

int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro("run.mac");
	G4String output("output");
	G4String patientPhantom;
	G4String doctorPhantom;
	G4UIExecutive* ui = 0;
	WEBServerConnect* serverConnect = 0;
	G4bool isSQL(false);
	G4bool isReset(false);
	int serverID(0);

	PrintUsage();
	for ( G4int i=1; i<argc; i++ ) 
	{
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[++i];
		}
		else if ( G4String(argv[i]) == "-o" ) {
			output = argv[++i];
		}
		else if ( G4String(argv[i]) == "-p" ) {
			patientPhantom = argv[++i];
		}
		else if ( G4String(argv[i]) == "-d" ) {
			doctorPhantom = argv[++i];
		}
		else if ( G4String(argv[i]) == "-v" )
		{
			ui = new G4UIExecutive(argc, argv, "Qt");
		}
		else if ( G4String(argv[i]) == "-sql" ) {
			isSQL = true;
		}
		else if ( G4String(argv[i]) == "-id") {
			serverID = stoi(argv[++i]);
		}
		else if ( G4String(argv[i]) == "-x") {
			isReset = true; // '-x' option set the system information (Flag, Hospital, Server, Time) to default value.
		}
		else {
			cout << "argument check" << endl;
			return 1;
		}
	}
		
	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(NULL));

	auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);
	runManager->SetNumberOfThreads(1);

	// Initialize Web-server connection
	if (isSQL) {
		serverConnect = new WEBServerConnect(serverID);
		if (isReset) {
			serverConnect->ResetData();
			delete serverConnect;
			return EXIT_SUCCESS;
		}
		G4cout << "START MCDCIR-PDC SYSTEM" << G4endl;
		G4cout << "  1) Server IP   :" << serverConnect->GetMyIPAddress() << G4endl;
		G4cout << "  2) Thread #    :" << serverConnect->GetThreadNo() << G4endl;
		G4cout << "  3) NPS         :" << serverConnect->GetNPS() << G4endl;
		runManager->SetNumberOfThreads(serverConnect->GetThreadNo());
	}

	// Import phantom
	TETModelImport* tetData_patient = nullptr;
	TETModelImport* tetData_doctor  = nullptr;
	if (!patientPhantom.empty())
		tetData_patient = new TETModelImport(patientPhantom, ui);
	else G4cout << "patient phantom is empty" << G4endl;
	if (!doctorPhantom.empty())
		tetData_doctor = new TETModelImport(doctorPhantom, ui);
	else G4cout << "doctor phantom is empty" << G4endl;

	// Set mandatory initialization classes
	auto det = new DetectorConstruction(tetData_patient, ui);
	ParallelPhantom* parallelPhantom;
	if (!doctorPhantom.empty()) {
		parallelPhantom = new ParallelPhantom("parallelPhantom", tetData_doctor, ui);
		det->RegisterParallelWorld(parallelPhantom);
	}
	runManager->SetUserInitialization(det);	

	// Physics List
	G4VModularPhysicsList* physicsList = new G4VModularPhysicsList();
	physicsList->RegisterPhysics(new G4EmLivermorePhysics());
	physicsList->SetCutValue(1*um, "gamma");
	physicsList->SetCutValue(1*um, "e-");
	physicsList->SetCutValue(1*um, "e+");
	if (!doctorPhantom.empty())
		physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallelPhantom", true));
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization(tetData_doctor, output, initTimer, serverConnect));

	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager virtual
	G4UImanager* UImanager = G4UImanager::GetUIpointer();

	if (isSQL) {
		// PDC system
		runManager->Initialize();
		while (true)
		{
			serverConnect->ReadCalculationNumber();
			if ( serverConnect->ProcessData() == SERV_STATUS::ALLOCATED )
				continue;
			
			// Get frame data from web-server and apply to messenger class
			// ===============================================================
			string hospital          = serverConnect->GetHospitalName();
			int frameNo              = serverConnect->GetFrameNo();
			double frameTime         = serverConnect->GetFrameTime();
			int tubeVoltage          = serverConnect->GetTubeVoltage();
			int tubeCurrent          = serverConnect->GetTubeCurrent();
			double DAP               = serverConnect->GetDAP();
			int FD                   = serverConnect->GetFD();
			int RAOLAO               = serverConnect->GetRAOLAO();
			int CRANCAUD             = serverConnect->GetCRANCAUD();
			int SID                  = serverConnect->GetSID();
			Vector3d table_trans     = serverConnect->GetTableTranslation();
			int table_theta          = serverConnect->GetTableTheta();
			Vector3d glass_trans     = serverConnect->GetGlassTranslation();
			Quaterniond glass_quat   = serverConnect->GetGlassQuaternion();
			RotationList doctor_vQ   = serverConnect->GetDoctorRotationList();
			Vector3d     doctor_root = serverConnect->GetDoctorRoot();
			// ===============================================================
			if (!doctorPhantom.empty())
				parallelPhantom->Deform(doctor_vQ, doctor_root);
			UImanager->ApplyCommand("/beam/peak " + to_string(tubeVoltage));
			UImanager->ApplyCommand("/beam/FD FD" + to_string(FD));
			UImanager->ApplyCommand("/beam/c-arm " + to_string(RAOLAO) + " " + to_string(CRANCAUD) + " " + to_string(SID));

			runManager->BeamOn(serverConnect->GetNPS());
			serverConnect->UpdateFlagToDone();

			break;
			sleep(10); // standby
		}
	}
	else 
	{

		if ( ! ui )
		{
			// batch mode
			G4String command = "/control/execute ";
			UImanager->ApplyCommand(command+macro);
		}
		else 
		{
			// interactive mode
			UImanager->ApplyCommand("/control/execute init_vis.mac");
			ui->SessionStart();
			delete ui;
		}
	}

	delete visManager;
	delete runManager;
	delete serverConnect;

    cout << "EXIT_SUCCESS" << endl;
    return EXIT_SUCCESS;
}
