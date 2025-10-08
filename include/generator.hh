#ifndef GENERATOR_HH
#define GENERATOR_HH

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4GeneralParticleSource.hh"

#include "G4ParticleGun.hh"
#include "G4PhysicsFreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleMomentum.hh"
//#include "SBG4SPSEneDistribution.hh"
//#include "G4SPSRandomGenerator.hh"
#include "G4GenericMessenger.hh"

//#include "functions.hh"

#include <vector>
#include <map>
#include <random>
#include <iostream>
#include <cmath>
#include <chrono>

//#include "ions.hh"
//#include "GCRFlux.hh"

#define PI 3.1415926535926535



class G4GeneralParticleSource;

G4double GetAngle(G4ThreeVector,G4ThreeVector);

class MyPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
	public :
		MyPrimaryGenerator();
		~MyPrimaryGenerator();
		
		virtual void GeneratePrimaries(G4Event*);

		// Accessor
		G4double GetBeamSurface() const {return beamSurface;}
		G4double GetRSource() const {return rsource;}
		//G4double GetTotalParticleNumber(G4int Z) const {return gcrFlux->GetParticleFlux(Z);}
		//G4double GetMissionDuration() const {return secondToYear;}
		G4double GetMissionFactor() const {return missionFactor;}
		//G4String GetIonName(G4int i) const { return gcrFlux->GetIonName(i);}
		//G4int GetNIons() const { return gcrFlux->GetNIons();}

		// Return the flux for an ions Z with an arbitrary energy E
		//G4double GetEnergyFlux(G4int Z, G4double E) const {return gcrFlux->GetEnergyFlux(Z,E);};
		//G4int GetNGenerated(G4int i) const {return Npart[i];}
		
		//G4int GetSampleSize() const {return sampleSize;}
		G4int GetNGenerated() const {return Npart;}
		
		
	private :
	
		G4GenericMessenger *fMess;
		
		G4ParticleGun *fParticleGun;
		G4bool halfSphere;
  		G4double rsource, beamSurface, lowPosZ,factorSphere;
		G4double missionFactor;

		//G4double secondToYear;

		//G4SPSRandomGenerator *biasRndm;

		//const GCRFlux *gcrFlux;
		//const G4String fluxFile;

		G4IonTable *ionTable;

		G4int idParticle,sampleSize;
		//G4int nevent,nrejected;

		G4int iNbin;
		G4double ilogemin,ilogemax;
		G4double ilogemin_gen,ilogemax_gen;

		//SBG4SPSEneDistribution *eneGenerator;

		//std::vector<G4double> weightToUse;
		//std::vector<G4double> ratioPart;
		G4int Npart;

		
		std::vector <G4ThreeVector> cmpVect;
		std::vector<G4int> Nangles;

		G4ThreeVector GenMomentum(G4ThreeVector pos);

		// Commands for macros
		void DefineCommands();
		//void SetSampleSize(G4int s) {sampleSize=s;}
		void SetBeamRadius(G4double r);
		//void SetParticleRatio(G4double,G4double);
	
};


#endif
