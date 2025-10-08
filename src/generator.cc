
#include "generator.hh"
#include "Randomize.hh"
//#include "ions.hh"

#include "G4AnalysisManager.hh"

#include <math.h> 


//Constructor, take gcrFlux class by default
MyPrimaryGenerator::MyPrimaryGenerator()//:gcrFlux(new GCRFlux())
{
	fParticleGun = new G4ParticleGun(1);

	/*halfSphere = false;

	//nevent = 0;
	//nrejected = 0;

	rsource = 2.1*m;
	rsource = 3.5*m;

	SetBeamRadius(rsource);*/

	
	//iNbin = 100; //bins.GetNbins();
	//ilogemin = 1; //logMeV //bins.GetMinE(); //1 eV
	//ilogemax = 5; //logMeV //bins.GetMaxE(); //100 GeV

	//ilogemin_gen = 1; //logMeV //bins.GetMinE(); //1 eV
	//ilogemax_gen = 5; //logMeV //bins.GetMaxE(); //100 GeV

	
	// Prepare table of ions
	//ionTable = G4IonTable::GetIonTable();
	//ionTable->CreateAllIon();

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4String particleName = "proton";
	G4ParticleDefinition *particle = particleTable->FindParticle("proton");

	G4ThreeVector pos(0.,0.,1.*m);
	G4ThreeVector mom(0.,0.,1.);

	fParticleGun->SetParticlePosition(pos);
	fParticleGun->SetParticleMomentumDirection(mom);
	fParticleGun->SetParticleEnergy(1*MeV);
	fParticleGun->SetParticleDefinition(particle);

	/*cmpVect.push_back(G4ThreeVector(1,0,0));
	cmpVect.push_back(G4ThreeVector(-1,0,0));
	cmpVect.push_back(G4ThreeVector(0,1,0));
	cmpVect.push_back(G4ThreeVector(0,-1,0));
	cmpVect.push_back(G4ThreeVector(0,0,1));
	cmpVect.push_back(G4ThreeVector(0,0,-1));
	//cmpVect.push_back(G4ThreeVector(1,1,0));
	//cmpVect.push_back(G4ThreeVector(-1,-1,0));
	//cmpVect.push_back(G4ThreeVector(0,1,1));
	//cmpVect.push_back(G4ThreeVector(0,-1,-1));

	for (int i=0;i<cmpVect.size();i++) {Nangles.push_back(0);}*/

	// Preparing random generator
	//biasRndm = new G4SPSRandomGenerator();

	// Prepare energy distribution of all ions based on the modified class SBG4EneDistribution
	//for (G4int Z=1;Z<=gcrFlux->GetNIons();Z++)
	//{
	//G4String ionname(gcrFlux->GetIonName(Z-1));
	//eneGenerator = new SBG4SPSEneDistribution();
	//eneGenerator->SetBiasRndm(biasRndm);

	// Load energy distribution of specific ion Z
	//eneGenerator->SetEnergyDisType("Arb");
    	//    	listEneGenerator[ionname] = eneGenerator;
	//for (G4int ie=0;ie<gcrFlux->GetNbins(Z);ie++)
	//{
	//G4ThreeVector dat(gcrFlux->GetEneVal(Z,ie),gcrFlux->GetFluxVal(Z,ie));
	//listEneGenerator[ionname]->ArbEnergyHisto(dat);
	//	}
	//	listEneGenerator[ionname]->ArbInterpolate("Lin");
	//}
		
	// Number of ions used (usually 28)
	//G4int n_ions_kind = gcrFlux->GetNIons();
	
	// Initialize Nsimulated particle for each kind
	//for (int i=0;i<iNbin;i++) {Npart.push_back(0);}
	//Npart = 0;
	//for (int i=0;i<n_ions_kind;i++) {Npart.push_back(0);} 
	// Preparing ration of simulated particles, equal for all ions by default
	//for (int i=0;i<n_ions_kind;i++) {ratioPart.push_back(G4double(1.0));} 
	// Weigths of simulated particles (normalized ratios)
	//for (int i=0;i<n_ions_kind;i++) {weightToUse.push_back(G4double(1.0/n_ions_kind));} 
	
	// Size of sample, useful if many /run/beamOn are done in the same run.
	//sampleSize = 1e5;

	// Preparing commands concerning the generator for the macros
	//DefineCommands();

	// Generating random seed based on time measured in nanosecond to ensure distinct names 
	//	in case of multiparallel runs
	/*std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
	auto duration = now.time_since_epoch();
	auto nanoseconds = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);

	G4Random::setTheEngine(new CLHEP::RanecuEngine());
	G4long seed = nanoseconds.count()%10000000; //1;//time(NULL);
	G4Random::setTheSeed(seed);*/

}



MyPrimaryGenerator::~MyPrimaryGenerator()
{
	delete fParticleGun;
	//delete fMess;
	//delete Npart
}

//G4int MyPrimaryGenerator::GetNGenerated()
//{
//	return Npart;
//}

// Class used to generate one particle, it select the ion, the position where it originates, its momentum and energy
void MyPrimaryGenerator::GeneratePrimaries(G4Event* anEvent)
{

	// pick an ion in the list based on the probability defined in weightsToUse
	//idParticle = pickIndex(weightToUse);
	//idParticle = 0;
	//G4String particleName = gcrFlux->GetIonName(idParticle);
	//G4cout << "GeneratePrimaries" << G4endl;
	/*G4ParticleDefinition* particleDef = ionTable->GetIon(1,1,0*keV);
	
	// Define its energy
	
	G4double energy = CLHEP::RandFlat::shoot(ilogemin_gen,ilogemax_gen);
	//energy = pow(10,energy)*MeV;
	energy = 1*MeV;

	// Generate polar angle between -pi and pi
	G4double thetapos = CLHEP::RandFlat::shoot(-PI,PI);
	
	// Generate random z between ground and sphere of source
	G4double posz0 = CLHEP::RandFlat::shoot(lowPosZ,rsource);

	// Compute the x and y positions from z and theta
	G4double r = sqrt(rsource*rsource-posz0*posz0);
	G4double posx0 = r * cos(thetapos);
	G4double posy0 = r * sin(thetapos);

	//posx0 = 0;
	//posy0 = 0;
	//posz0 = 1.0*m;
	//posz0 = 282.0948*mm;
	//posz0 = rsource;


	// Define position
	G4ThreeVector pos(posx0,posy0,posz0);
	G4ThreeVector mom = GenMomentum(pos);
	Npart++;*/

	/*for (int i=0;i<3;i++)
	{
		G4cout << mom[i] << " ";
	}*/

	/*for (int i=0;i<cmpVect.size();i++)
	{

		//G4cout << acos(1.0-1.0/(2.0*PI)) << G4endl;
		if (abs(GetAngle(-mom,cmpVect[i]))<acos(1.0-1.0/(2.0*PI))) {Nangles[i]++;}
		//if (abs(GetAngle(-mom,cmpVect[i]))<PI/4.0) {Nangles[i]++;}
		//G4cout << cmpVect[i] << " ";
		//G4cout << GetAngle(-mom,cmpVect[i]) << G4endl;
	}


	for (int i=0;i<cmpVect.size();i++)
	{
		G4cout << (double)Nangles[i]/Npart << " " ;
	}
	G4cout << G4endl;
	*/
	


	// Define momentum within cosine distribution
	
      /*G4double angle = 0;
      G4double magPos=0;
      G4double magMom=0;


	G4double sign = abs(mom[2])/mom[2];
	//pos = G4ThreeVector(0,0,1);
      	for (G4int i=0;i<3;i++)
      	{
		angle += pos[i]*mom[i];
        	magPos += pos[i]*pos[i];
        	magMom += mom[i]*mom[i];
	}

	magPos = sqrt(magPos);
	magMom = sqrt(magMom);

	angle = sign*acos(abs(angle)/(magPos*magMom));*/
	//G4cout << angle << ",";
	//++nevent;

	//G4cout << mom.x() << G4endl;

	// If momentum is going upward repeat the process of spacial selection of particle parameters

	//mom = G4ThreeVector(-posx0,-posy0,-posz0);
	/*G4cout << (G4double) (nrejected)/(G4double)(nevent) << G4endl;

	G4double theta = acos(pos[2]/sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]));
	G4double phi = atan(pos[1]/pos[0]);

	std::ofstream outFile("../genDistrib.txt",std::ios::app);
	outFile << theta << "," << phi << "," << pos[0] << "," << pos[1] << "," << pos[2] << std::endl;
	outFile.close();*/

	

	// By now the particle will be generated and is accounted for simulations
	//G4double logprimkE = log10(energy);	
	//G4int iebin = iNbin*(logprimkE-ilogemin)/(ilogemax-ilogemin);	
	//Npart++;
	//G4cout << "N gen " << iebin << "  " << Npart[iebin] << G4endl;


	// Assign parameters to the gun
	//fParticleGun->SetParticlePosition(pos);
	//fParticleGun->SetParticleMomentumDirection(mom);
	//fParticleGun->SetParticleEnergy(energy);
	//fParticleGun->SetParticleDefinition(particleDef);
	fParticleGun->GeneratePrimaryVertex(anEvent);

	//G4cout << "Gun generated " << G4endl;
}



// Function generating a cosine distribution from the source point
G4ThreeVector MyPrimaryGenerator::GenMomentum(G4ThreeVector pos)
{
	// Normalize pos
	G4double normpos = sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2]);
	G4double normposx = pos[0]/normpos;
	G4double normposy = pos[1]/normpos;
	G4double normposz = pos[2]/normpos;
	//G4cout << normposx << " " << normposy << " " << normposz << G4endl;

	G4double u = CLHEP::RandFlat::shoot(0.0,1.0);
	G4double v = CLHEP::RandFlat::shoot(0.0,1.0);

	//G4cout << "u = " << u << "  v=" << v << G4endl;

	// theta and phi in spherical coordinates
	G4double theta = asin(u);
	//theta = PI / 2.0 * u;
	G4double phi = 2 * PI * v;

	//G4cout << theta << " " << phi << G4endl;

	// generated vector along the z axis
	G4double x = sin(theta)*cos(phi);
	G4double y = sin(theta)*sin(phi);
	G4double z = cos(theta);

	//x = CLHEP::RandFlat::shoot(-1.0,1.0);
	//y = CLHEP::RandFlat::shoot(-1.0,1.0);
	//z = CLHEP::RandFlat::shoot(-1.0,1.0);
	//G4double normXYZ = sqrt(x*x + y*y + z*z);
	//x = x/normXYZ;
	//y = y/normXYZ;
	//z = z/normXYZ;

	//return G4ThreeVector(x,y,z);
	

	//G4cout << x << " " << y << " " << z << G4endl;

	// Calculate vparallel
	G4double xpar = normposy;
	G4double ypar = -normposx;
	G4double zpar = 0.0;
	G4double normpar = sqrt(xpar*xpar + ypar*ypar + zpar*zpar);
	if (normpar==0)
	{
		xpar = -normposz;
		ypar = 0;
		zpar = normposx;
	}
	normpar = sqrt(xpar*xpar + ypar*ypar + zpar*zpar);
	xpar = xpar/normpar;
	ypar = ypar/normpar;
	zpar = zpar/normpar;
	//G4cout << xpar << " " << ypar << " " << zpar << " " << normpar << G4endl;

	// Calculate vperp
	G4double xper = normposy*zpar - normposz*ypar;
	G4double yper = normposz*xpar - normposx*zpar;
	G4double zper = normposx*ypar - normposy*xpar;
	


	//G4double xper = normposx*normposy;
	//G4double yper = normposy*normposy;
	//G4double zper = -normposx*normposz-normposx*normposy;
	G4double normper = sqrt(xper*xper + yper*yper + zper*zper);
	xper = xper/normper;
	yper = yper/normper;
	zper = zper/normper;
	//G4cout << xper << " " << yper << " " << zper << G4endl;

	//New vector
	G4double xrot = x * xpar + y * xper + z * normposx;
	G4double yrot = x * ypar + y * yper + z * normposy;
	G4double zrot = x * zpar + y * zper + z * normposz;
	G4double normrot = sqrt(xrot*xrot + yrot*yrot + zrot*zrot);
	xrot = -xrot/normrot;
	yrot = -yrot/normrot;
	zrot = -zrot/normrot;
	//G4cout << xrot << " " << yrot << " " << zrot << G4endl;
	
	return G4ThreeVector(xrot,yrot,zrot);

	
	//Calculate vparallel
	//G4double randpar = CLHEP::RandFlat::shoot(0.0,1.0);
	/*G4double randpar = CLHEP::RandFlat::shoot(0.0,1.0);
	G4double sizevpar = -sqrt(1-randpar);
	G4double xpar = normposx * sizevpar;
	G4double ypar = normposy * sizevpar;
	G4double zpar = normposz * sizevpar;

	
	// Calculate vperpendicular
	// First we need a base in the perpendicular plane
	// First base
	G4double normb1 = sqrt(zpar*zpar+xpar*xpar);
	G4double b1x = zpar/normb1;
	G4double b1y = 0.0;
	G4double b1z = -xpar/normb1;

	//Second base
	G4double b2x = normposy*b1z - normposz*b1y;
	G4double b2y = normposz*b1x - normposx*b1z;
	G4double b2z = normposx*b1y - normposy*b1x;

	
	G4double sizevperp = sqrt(1.0 - sizevpar*sizevpar);	

	G4double theta = CLHEP::RandFlat::shoot(-PI,PI);
	G4double mx0rot = sizevperp * (b1x * cos(theta) + b2x * sin(theta)) + xpar;
	G4double my0rot = sizevperp * (b1y * cos(theta) + b2y * sin(theta)) + ypar;
	G4double mz0rot = sizevperp * (b1z * cos(theta) + b2z * sin(theta)) + zpar;

	return G4ThreeVector(mx0rot,my0rot,mz0rot);*/

}

/*void MyPrimaryGenerator::SetParticleRatio(G4double Z, G4double rat)
{
	ratioPart[G4int(Z)-1] = rat;
	
	// Normalize ratio part to get weightToUse
	G4double tot = 0.;
	for (G4int i=0;i<weightToUse.size();i++){tot+=ratioPart[i];}
	if (tot==0) return;
	for (G4int i=0;i<weightToUse.size();i++){weightToUse[i]=ratioPart[i]/tot;}
}*/

// Prepare everything that depends on the radius of the beam source
void MyPrimaryGenerator::SetBeamRadius(G4double r)
{
	rsource = r;
	
	//G4cout << "SetBeamRadius  " << rsource << G4endl;
	beamSurface = 4.0 * PI * rsource*rsource/m/m;
	lowPosZ = -rsource;
	factorSphere = 1.0; //Case GCR is taken from all directions

	// If scene is of the Lunar surface (it is), the starting surface is half sphere, the mimal z position is 0, 
	//   and we are only considering half of the possible origins (the other being blocked by the Moon).
	missionFactor = factorSphere*beamSurface*2.0*PI;
}

// Macros commands:
// 	-- sampleSize (= how many particles are generated before values are averaged and recorded)
//	-- radbeam: define the radius of the beam
//	-- ratio of particle generated
void MyPrimaryGenerator::DefineCommands()
{
	//fMess = new G4GenericMessenger(this, "/SIM/scoring/","Set sample size");
	//auto& samplesizeCmd = fMess->DeclareMethod("sampleSize",&MyPrimaryGenerator::SetSampleSize,"Set size of sample");
	//samplesizeCmd.SetParameterName("sampleSize", true);
	//samplesizeCmd.SetDefaultValue("1e5");
	
	fMess = new G4GenericMessenger(this, "/SIM/scoring/","Set radius beam");
	auto& radiusbeamCmd = fMess->DeclareMethodWithUnit("radbeam","mm",&MyPrimaryGenerator::SetBeamRadius,"Set radius beam");
	radiusbeamCmd.SetParameterName("radbeam", true);
	radiusbeamCmd.SetDefaultValue("3004*mm");

	//fMess = new G4GenericMessenger(this, "/SIM/generate/","Set particle ratio");
	//auto& partRatioCmd = fMess->DeclareMethod("rat",&MyPrimaryGenerator::SetParticleRatio,"Set particle ratio");
	//partRatioCmd.SetParameterName("rat", true);
	//partRatioCmd.SetDefaultValue("1 1");
}



G4double GetAngle(G4ThreeVector v1,G4ThreeVector v2)
{
	G4double angle(0.0);
	G4double mag1(0);
	G4double mag2(0);

	for (int i=0;i<3;i++)
	{
		mag1 += v1[i]*v1[i];
		mag2 += v2[i]*v2[i];
		angle += v1[i]*v2[i];
	}

	mag1 = sqrt(mag1);
	mag2 = sqrt(mag2);

	

	//G4cout << v1 << " " << v2 << " " << acos(angle/(mag1*mag2)) << " " << (acos(angle/(mag1*mag2))<PI/4.0) << G4endl;

	return acos(angle/(mag1*mag2));

};







