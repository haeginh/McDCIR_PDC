// Geant4
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

#include "DetectorConstruction.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "ActionInitialization.hh"
#include "parallelmesh.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"


int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output;
	G4String phantomName;
	G4UIExecutive* ui = 0;


	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) {
			macro = argv[++i];
		}
		// output file name
		// else if ( G4String(argv[i]) == "-o" ) {
		// 	output = argv[++i];
		// }
		else if ( G4String(argv[i]) == "-v" )
		{
			ui = new G4UIExecutive(argc, argv, "Qt");
		}
		else {
			G4cout << "argument check" << G4endl;
			return 1;
		}
	}

	// default output file name
	// if ( !output.size() ) output = macro + ".out";

	auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// Set mandatory initialization classes
	auto det = new DetectorConstruction();
	det->RegisterParallelWorld(new ParallelMesh("parallel"));
	runManager->SetUserInitialization(det);
	G4VModularPhysicsList* physicsList = new FTFP_BERT;
	physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallel"));
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization());

	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager  virtual
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if ( ! ui ){
		// batch mode
		G4String command = "/control/execute ";
		UImanager->ApplyCommand(command+macro);
	}
	else {
		// interactive mode
		UImanager->ApplyCommand("/control/execute init_vis.mac");
		// runManager->Initialize();
		ui->SessionStart();
		delete ui;
	}

	delete visManager;
	delete runManager;

    G4cout << "EXIT_SUCCESS" << G4endl;
    return EXIT_SUCCESS;
}
