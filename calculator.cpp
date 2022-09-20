// Geant4
#include "G4UImanager.hh"
#include "G4RunManagerFactory.hh"

#include "DetectorConstruction.hh"
#include "ParallelGlass.hh"
#include "ParallelPhantom.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4ParallelWorldPhysics.hh"
#include "ActionInitialization.hh"

#include "G4Timer.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "Randomize.hh"

#include "PhysicsList.hh"
#include "PhantomData.hh"
#include "PhantomAnimator.hh"

// void SetFrameData(std::vector<G4float>);
int main(int argc, char** argv)
{
	G4Timer* initTimer = new G4Timer;
	initTimer->Start();
	G4UIExecutive* ui = 0;
	G4String macro;

	for ( G4int i=1; i<argc; i++ ) {
		// visualization
		if ( G4String(argv[i]) == "-v" )
		{
			ui = new G4UIExecutive(argc, argv, "Qt");
		}
		else if ( G4String(argv[i]) == "-m" )
		{
			macro = argv[++i];
		}
		else {
			cout << "argument check" << endl;
			return 1;
		}
	}

	auto* runManager = G4RunManagerFactory::CreateRunManager(G4RunManagerType::MT);

	// Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	G4Random::setTheSeed(time(0));

	// phantom-related
	PhantomData* phantomData_AM = new PhantomData("AM");
	G4int maxNum(5);
	if(std::getenv("DCIR_MAX_NUM_OF_PEOPLE")) maxNum = std::atoi(std::getenv("DCIR_MAX_NUM_OF_PEOPLE"));
	std::vector<ParallelPhantom*> parallelworlds; // phantoms  + glass
	for(int i=0;i<maxNum;i++)
		parallelworlds.push_back(new ParallelPhantom("phantom"+std::to_string(i), phantomData_AM));
	auto det = new DetectorConstruction();
	det->SetPatientName("Patient_M_H175M83800_fixed");
	auto physicsList = new PhysicsList();
	for(int i=maxNum-1;i>-1;i--) //register in opposite order
	{
		det->RegisterParallelWorld(parallelworlds[i]);
		physicsList->RegisterPhysics(new G4ParallelWorldPhysics(parallelworlds[i]->GetName(), true));
	}
	auto parallelGlass = new ParallelGlass("parallelGlass", "pbGlass.1");
	det->RegisterParallelWorld(parallelGlass);
	
	// Set mandatory initialization classes
	runManager->SetUserInitialization(det);
	// physicsList->RegisterPhysics(new G4ParallelWorldPhysics("parallelGlass", true));
	runManager->SetUserInitialization(physicsList);
	runManager->SetUserInitialization(new ActionInitialization(parallelworlds, initTimer));

	// runManager->Initialize();
	// Initialize visualization
	G4VisManager* visManager = new G4VisExecutive("Quiet");
	visManager->Initialize();

	// Get the pointer to the User Interface manager  virtual
	G4UImanager* UImanager = G4UImanager::GetUIpointer();
	if ( ! ui ){
		// runManager->Initialize();
		UImanager->ApplyCommand("/control/execute " + macro);
		while(true)
		{
			
			//try to read data from NAS
			//if(quit sign) break;
			//else if(wait sign) {sleep(1); continue;}
			// std::vector<G4float> frameData;
			// SetFrameData(frameData);


			break; //erase later
		}
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
