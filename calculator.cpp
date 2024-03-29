// Geant4
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

#include "DetectorConstruction.hh"
// #include "ParallelGlass.hh"
#include "ParallelPhantom.hh"
#include "G4PhysListFactory.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "ActionInitialization.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "PhantomData.hh"
#include "PhantomAnimator.hh"

int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4String macro;
	G4String output;
	G4UIExecutive* ui = 0;

	for ( G4int i=1; i<argc; i++ ) {
		// macro file name
		if ( G4String(argv[i]) == "-m" ) 
			macro = argv[++i];
		// output file name
		else if ( G4String(argv[i]) == "-o" ) 
			output = argv[++i];
		// phantom initialization data
		else if ( G4String(argv[i]) == "-v" )
			ui = new G4UIExecutive(argc, argv);
		else {
			cout << "argument check" << endl;
			return 1;
		}
	}

	// default output file name
	if ( !output.size() ) output = macro + ".out";

G4cout<<1<<G4endl;
	auto* runManager = G4RunManagerFactory::CreateRunManager();
	      G4int nThreads = 4;
       runManager->SetNumberOfThreads(nThreads);
G4cout<<1<<G4endl;

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));
	// Import phantom
	PhantomData* amData = new PhantomData("MRCP_AM");
	PhantomData* afData = new PhantomData("MRCP_AF");
G4cout<<1<<G4endl;

	// Set mandatory initialization classes
	auto det = new DetectorConstruction();
	det->SetPatientName("Patient_M_H175M83800_fixed");
	// det->RegisterParallelWorld(new ParallelGlass("parallelGlass", "./phantoms/glassTet.1"));
	det->RegisterParallelWorld(new ParallelPhantom("parallelPhantom", amData, afData));
	runManager->SetUserInitialization(det);	
	auto factory = new G4PhysListFactory();
    G4VModularPhysicsList* physicsList = factory->GetReferencePhysList("FTFP_BERT_LIV");
	// physicsList->RegisterPhysics(new G4StepLimiterPhysics());
	// physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallelGlass", true));
	physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallelPhantom", true));
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization(output, initTimer));

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

    cout << "EXIT_SUCCESS" << endl;
    return EXIT_SUCCESS;
}
