#ifndef GEOMETRY_HH
#define GEOMETRY_HH

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4CSGSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4NistManager.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"

#include "G4SystemOfUnits.hh"
#include "G4MultiFunctionalDetector.hh"

#include "SBG4PSSphereSurfaceFlux.hh"

#include "G4VPrimitiveScorer.hh"
#include "G4RunMessenger.hh"

#include <sstream>


#include "G4SDManager.hh"

#include "G4RunManager.hh"

class MyGeometry : public G4VUserDetectorConstruction
{
	public:		
		MyGeometry();
		~MyGeometry();
		
		virtual G4VPhysicalVolume *Construct();
	
	private:
		G4bool checkOverLaps;
		
	protected:
		//Materials *materials;
	
		G4double rPhantom ;
		//G4LogicalVolume *fScoringVolume;
		//G4GenericMessenger *fMessenger;
		
		G4Sphere *solidFlux;
		G4LogicalVolume *logicFluxSphere;
		G4VPhysicalVolume *physFluxSphere;
		//TETModelImport *fTetData;
		//G4LogicalVolume *logicWorld, *logicAir, *logicDome1, *logicDome2, *logicDome3, *logicDome4, *logicPhantom, *fTetLogic, *fContainer_logic, *logicGround;
		G4Box *solidWorld;
		G4LogicalVolume *logicWorld;
		G4VPhysicalVolume *physWorld;

		// rPhantom stands for the radius of the Human Phantom
		//G4double zWorld;
		//G4ThreeVector      fPhantomSize;
		//G4ThreeVector      fPhantomBoxMin, fPhantomBoxMax;
		//G4int              fNOfTetrahedrons;

		virtual void ConstructSDandField();	

		
};



#endif
