
#include "generator.hh"
#include "Randomize.hh"
//#include "ions.hh"

#include "G4AnalysisManager.hh"

#include <math.h> 


//Constructor, take gcrFlux class by default
MyPrimaryGenerator::MyPrimaryGenerator()//:gcrFlux(new GCRFlux())
{
	fParticleGun = new G4GeneralParticleSource();
}



MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
}

// Class used to generate one particle, it select the ion, the position where it originates, its momentum and energy
void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{
	fParticleGun->GeneratePrimaryVertex(anEvent);
}




