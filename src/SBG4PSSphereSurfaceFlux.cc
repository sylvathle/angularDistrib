//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
// SBG4PSSphereSurfaceFlux
#include "SBG4PSSphereSurfaceFlux.hh"

#include "G4SystemOfUnits.hh"
#include "G4StepStatus.hh"
#include "G4Track.hh"
#include "G4VSolid.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VPVParameterisation.hh"
#include "G4UnitsTable.hh"
#include "G4GeometryTolerance.hh"
////////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for scoring only Surface Flux.
//  Flux version assumes only for G4Sphere shape.
//
// Surface is defined  at the inside of sphere.
// Direction                  -Rmin   +Rmax
//   0  IN || OUT            ->|<-     |
//   1  IN                   ->|       |
//   2  OUT                    |<-     |
//
// Created: 2005-11-14  Tsukasa ASO, Akinori Kimura.
// 29-Mar-2007  T.Aso,  Bug fix for momentum direction at outgoing flux.
// 2010-07-22   Introduce Unit specification.
// 2010-07-22   Add weighted and divideByAre options
// 2011-02-21   Get correct momentum direction in Flux_Out.
// 2011-09-09   Modify comment in PrintAll().
// 2014-03-03  T.Aso,  To use always positive value for anglefactor.
///////////////////////////////////////////////////////////////////////////////

SBG4PSSphereSurfaceFlux::SBG4PSSphereSurfaceFlux(G4String name, G4int direction,
                                             G4int depth)
  : SBG4PSSphereSurfaceFlux(name, direction, "percm2", depth)
{
}

SBG4PSSphereSurfaceFlux::SBG4PSSphereSurfaceFlux(G4String name, G4int direction,
                                             const G4String& unit, G4int depth)
  : G4VPrimitiveScorer(name, depth)
  , HCID(-1)
  , fDirection(direction)
  , EvtMap(nullptr)
  , weighted(true)
  , divideByArea(true)
{
  G4cout << "Create SBG4PSSphereSurfaceFlux" << G4endl;

  for (G4int i=0;i<Nabins;i++)
  {
	N_in_angle.push_back(0);
  }
  

  DefineUnitAndCategory();
  SetUnit(unit);
}


SBG4PSSphereSurfaceFlux::~SBG4PSSphereSurfaceFlux()
{
	G4cout << "Destroy out" << G4endl;
}

G4bool SBG4PSSphereSurfaceFlux::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  //if (aStep->GetTrack()->GetParticleDefinition()->GetParticleName()!="neutron") return false;
  //G4cout << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << G4endl;
  
  G4StepPoint* preStep = aStep->GetPreStepPoint();

  G4VPhysicalVolume* physVol       = preStep->GetPhysicalVolume();
  G4VPVParameterisation* physParam = physVol->GetParameterisation();
  G4VSolid* solid                  = nullptr;



  if(physParam != nullptr)
  {  // for parameterized volume
    G4int idx =
      ((G4TouchableHistory*) (aStep->GetPreStepPoint()->GetTouchable()))
        ->GetReplicaNumber(indexDepth);
    solid = physParam->ComputeSolid(idx, physVol);
    solid->ComputeDimensions(physParam, idx, physVol);
  }
  else
  {  // for ordinary volume
    solid = physVol->GetLogicalVolume()->GetSolid();
  }

  auto  sphereSolid = (G4Sphere*) (solid);

  G4int dirFlag = IsSelectedSurface(aStep, sphereSolid);
  if(dirFlag > 0)
  {
    if(fDirection == fFlux_InOut || fDirection == dirFlag)
    {
      G4StepPoint* thisStep = nullptr;
      if(dirFlag == fFlux_In)
      {
        thisStep = preStep;
      }
      else if(dirFlag == fFlux_Out)
      {
        thisStep = aStep->GetPostStepPoint();
      }
      else
      {
        return false;
      }

	G4double kE_up = 0.0; // Albedo
	G4double kE_down = 0.0;
	
	//G4cout << thisStep->GetPosition()[2] << " " << thisStep->GetMomentum()[2] << " " << preStep->GetKineticEnergy() << G4endl;

      G4ThreeVector position = thisStep->GetPosition();
      G4ThreeVector momentum = thisStep->GetMomentum();

      /*G4cout << " position: " << position[0] << " ";
      G4cout << position[1] << " ";
      G4cout << position[2] << G4endl;

      G4cout << " momentum: " << momentum[0] << " ";
      G4cout << momentum[1] << " ";
      G4cout << momentum[2] << G4endl;*/

      G4double angle = 0;
      G4double magPos=0;
      G4double magMom=0;


	G4double sign = abs(momentum[2])/momentum[2];
	position = G4ThreeVector(0,0,1);
      	for (G4int i=0;i<3;i++)
      	{
		angle += position[i]*momentum[i];
        	magPos += position[i]*position[i];
        	magMom += momentum[i]*momentum[i];
	}

	magPos = sqrt(magPos);
	magMom = sqrt(magMom);

	angle = sign*acos(abs(angle)/(magPos*magMom));
        //G4cout << magPos << ",";
        G4cout << momentum[2] << ",";

	G4int ibin = G4int (angle/(3.1415926535/2.0)*Nabins);
	N_in_angle[ibin]++;
	/*G4cout << "[";
        for (G4int ii=0; ii<Nabins-1; ii++)
	{
		G4cout << N_in_angle[ii] << ",";
	}
	G4cout << N_in_angle[Nabins-1] << "]";
	G4cout << G4endl;*/
	
	

	if (angle > 3.1415926535/2.0) {G4cout << "BIG ANGLE " << angle << G4endl;}
	if (angle > 0.5719537610547146) {return false;}

      //G4cout << "angle: " << angle << G4endl;

	
 
      if (thisStep->GetMomentum()[2]<0) {kE_down = preStep->GetKineticEnergy();}
      else {kE_up = preStep->GetKineticEnergy();}
     
      /*if (aStep->GetTrack()->GetParticleDefinition() == G4Electron::Electron())
      {
        G4double kE = aStep->GetTrack()->GetKineticEnergy();
        //if ((kE<80*keV) && (kE>50*keV)) 
        {
          aStep->GetTrack()->SetTrackStatus(fStopAndKill);
          return true;
        }
      }*/

      fluxInfo flux = fluxInfo(kE_down,kE_up);
      
      index++;

      EvtMap->add(index, flux);
    }
  }

  return true;
}

G4int SBG4PSSphereSurfaceFlux::IsSelectedSurface(G4Step* aStep,
                                               G4Sphere* sphereSolid)
{
  G4TouchableHandle theTouchable =
    aStep->GetPreStepPoint()->GetTouchableHandle();
  G4double kCarTolerance =
    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();


  if(aStep->GetPreStepPoint()->GetStepStatus() == fGeomBoundary)
  {
    // Entering Geometry
    G4ThreeVector stppos1 = aStep->GetPreStepPoint()->GetPosition();
    G4ThreeVector localpos1 =
      theTouchable->GetHistory()->GetTopTransform().TransformPoint(stppos1);
    G4double localR2 = localpos1.x() * localpos1.x() +
                       localpos1.y() * localpos1.y() +
                       localpos1.z() * localpos1.z();
    G4double InsideRadius = sphereSolid->GetOuterRadius();
    if(localR2 >
         (InsideRadius - kCarTolerance) * (InsideRadius - kCarTolerance) &&
       localR2 <
         (InsideRadius + kCarTolerance) * (InsideRadius + kCarTolerance))
    {
      return fFlux_In;
    }
  }

  return -1;
}

void SBG4PSSphereSurfaceFlux::Initialize(G4HCofThisEvent* HCE)
{
  EvtMap = new G4THitsMap<fluxInfo>(detector->GetName(), GetName());
  index = 0;
  if(HCID < 0)
    HCID = GetCollectionID(0);
  HCE->AddHitsCollection(HCID, (G4VHitsCollection*) EvtMap);
}

void SBG4PSSphereSurfaceFlux::clear() { EvtMap->clear(); }

/*void SBG4PSSphereSurfaceFlux::PrintAll()
{
  G4cout << " MultiFunctionalDet  " << detector->GetName() << G4endl;
  G4cout << " PrimitiveScorer " << GetName() << G4endl;
  G4cout << " Number of entries " << EvtMap->entries() << G4endl;
  for(const auto& [copy, flux] : *(EvtMap->GetMap()))
  {
    G4cout << "  copy no.: " << copy
           << "  Flux  : " << (flux)->kE << " " << (flux)->particleName << " ["
           << GetUnit() << "]" << G4endl;
  }
}*/

void SBG4PSSphereSurfaceFlux::SetUnit(const G4String& unit)
{
  if(divideByArea)
  {
    CheckAndSetUnit(unit, "Per Unit Surface");
  }
  else
  {
    if(unit.empty())
    {
      unitName  = unit;
      unitValue = 1.0;
    }
    else
    {
      G4String msg = "Invalid unit [" + unit + "] (Current  unit is [" +
                     GetUnit() + "] ) for " + GetName();
      G4Exception("SBG4PSSphereSurfaceFlux::SetUnit", "DetPS0016", JustWarning,
                  msg);
    }
  }
}

void SBG4PSSphereSurfaceFlux::DefineUnitAndCategory()
{
  // Per Unit Surface
  new G4UnitDefinition("percentimeter2", "percm2", "Per Unit Surface",
                       (1. / cm2));
  new G4UnitDefinition("permillimeter2", "permm2", "Per Unit Surface",
                       (1. / mm2));
  new G4UnitDefinition("permeter2", "perm2", "Per Unit Surface", (1. / m2));
}
