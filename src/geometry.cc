#include "geometry.hh"

MyGeometry::MyGeometry()
{
}

MyGeometry::~MyGeometry()
{
}

G4VPhysicalVolume *MyGeometry::Construct()
{
	G4NistManager *nist = G4NistManager::Instance();

	// Dimensions of the world
	G4double xWorld = 2.0*m;
	G4double yWorld = 2.0*m;
	G4double zWorld = 2.0*m;

	// Create world of simulation.
	solidWorld = new G4Box("solidWorld", xWorld, yWorld, zWorld);
	logicWorld = new G4LogicalVolume(solidWorld,nist->FindOrBuildMaterial("G4_Galactic"),"logicWorld");
	physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, checkOverLaps);
	
	//Sphere capturing the flux
	rPhantom = 0.200*m;
	solidFlux = new G4Sphere("solidFlux", 0., rPhantom, 0, 360*deg, 0, 180*deg);
	logicFluxSphere = new G4LogicalVolume(solidFlux, nist->FindOrBuildMaterial("G4_Galactic"), "logicFluxSphere");
	physFluxSphere = new G4PVPlacement(0, G4ThreeVector(0.0,0.0,0.0), logicFluxSphere, "physFlux", logicWorld, false, 0, checkOverLaps);

	return physWorld;
}

void MyGeometry::ConstructSDandField()
{
	G4MultiFunctionalDetector* myScorer = new G4MultiFunctionalDetector("myCellScorer");
	G4SDManager::GetSDMpointer()->AddNewDetector(myScorer);
	G4VPrimitiveScorer* totalSurfFlux = new SBG4PSSphereSurfaceFlux("TotalSurfFlux",0);
	myScorer->RegisterPrimitive(totalSurfFlux);
	logicFluxSphere->SetSensitiveDetector(myScorer);
}

