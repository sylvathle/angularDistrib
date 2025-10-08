#include <iostream>

//#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//#include "QBBC_HP.hh"
#include "QBBC.hh"
#include "QGSP_INCLXX_HP.hh"
#include "QGSP_BERT.hh"
#include "Shielding.hh"
#include "G4HadronicParameters.hh"

#include "geometry.hh"
#include "physics.hh"
#include "action.hh"

int main(int argc, char** argv)
{
	G4RunManager *runManager = new G4RunManager();

	runManager->SetVerboseLevel(0);

	runManager->SetUserInitialization(new MyGeometry());
	runManager->SetUserInitialization(new QBBC(0));
	runManager->SetUserInitialization(new MyActionInitialization());

	G4UIExecutive *ui = nullptr;
	G4UImanager *UImanager = G4UImanager::GetUIpointer();
	if (argc==1) {ui = new G4UIExecutive(argc,argv);}

	if ( ! ui)
	{
		G4String command = "/control/execute ";
		G4String fileName = argv[1];
		UImanager->ApplyCommand(command+fileName);

	}
	else
	{
		G4VisManager *visManager = new G4VisExecutive;
		visManager->SetVerboseLevel(0);
		visManager->Initialize();
		UImanager->ApplyCommand(G4String("/control/execute vis.mac"));
		ui->SessionStart();
		delete visManager;
	}

	return 0;

}
