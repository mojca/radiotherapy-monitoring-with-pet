//
// This code is based on example B3 from Geant4.
// The original licence statement from Geant4 is as follows:
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
/// \file RTPETSteppingAction.cc
/// \brief Implementation of the RTPETSteppingAction class

#include "RTPETSteppingAction.hh"
#include "RTPETDetectorConstruction.hh"
#include "RTPETFileManager.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETSteppingAction::RTPETSteppingAction(
  const RTPETDetectorConstruction* detectorConstruction,
  RTPETFileManager* fileManager)
: G4UserSteppingAction(),
  fDetConstruction(detectorConstruction),
  fFileManager(fileManager)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETSteppingAction::~RTPETSteppingAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// TODO: http://hypernews.slac.stanford.edu/HyperNews/geant4/get/eventtrackmanage/996/1.html
// http://geant4.web.cern.ch/geant4/support/faq.shtml#TRACK-1
void RTPETSteppingAction::UserSteppingAction(const G4Step* step)
{
  // Collect energy and track length step by step
  // see track/src/G4ParticleChange.cc

  // get the track of the current stp
  G4Track* track = step->GetTrack();
  // get the process name
  G4String processName  = step->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
  // get the particle name
  G4String particleName = track->GetDefinition()->GetParticleName();

  // create the energy spectrum on the top and the bottom of the phantom
  if(
    (processName == G4String("Transportation")) &&
    (particleName == G4String("gamma")) &&
    ((track->GetPosition().y()/CLHEP::cm == 10) || (track->GetPosition().y()/CLHEP::cm == -10))
  ) {
    G4double x = track->GetPosition().x()/CLHEP::cm;
    G4double y = track->GetPosition().y()/CLHEP::cm;
    G4double z = track->GetPosition().z()/CLHEP::cm;

    if(y*y == 100 && x*x < 4 && z*z < 4) {
      G4VPhysicalVolume* volume1 = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
      G4VPhysicalVolume* volume2 = step->GetPostStepPoint()->GetTouchable()->GetVolume();
      G4String name1(volume1->GetName());
      G4String name2(".");
      if(volume2) {
        name2 = volume2->GetName();
      }

      G4double energy = step->GetPreStepPoint()->GetKineticEnergy();
      int variant = 0;

      if(name1 == G4String("World") && name2 == G4String("Patient")) {
        // at the entrance to the phantom
        variant = 1;
      } else if(name2 == G4String("World") && name1 == G4String("Patient")) {
        // at the exit of the phantom
        variant = 2;
      }
      fFileManager->AddBeamEnergy(variant, energy);
    }
  }
  // return;

  // processName(s):
  // - Transportation
  // - compt
  // - phot
  // - conv (pair production)
  if((processName != G4String("Transportation"))) {
    // prestep point
    G4StepPoint* preStepPoint = step->GetPreStepPoint();

    G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
    // get volume of the current step
    G4VPhysicalVolume* volume = theTouchable->GetVolume();
    G4String volumeName   = volume->GetName();
    G4String particleName = track->GetDefinition()->GetParticleName();
    G4double density      = volume->GetLogicalVolume()->GetMaterial()->GetDensity();
  //G4double density      = preStepPoint->GetMaterial()->GetDensity();

    G4int icryst, iring, idet;

    // store annihilations
    if(processName == G4String("annihil") && volumeName == G4String("Patient")) {
      fFileManager->AddAnnihilation(track->GetPosition().x()/CLHEP::cm, track->GetPosition().y()/CLHEP::cm, track->GetPosition().z()/CLHEP::cm);
    }
    // store energy deposits
    if(volumeName == G4String("Patient") && step->GetTotalEnergyDeposit() > 0) {
      fFileManager->AddEnergyDeposition(track->GetPosition().x()/CLHEP::cm, track->GetPosition().y()/CLHEP::cm, track->GetPosition().z()/CLHEP::cm, step->GetTotalEnergyDeposit(), density);
    }

    // some energy deposition in the crystal
    if((volumeName == G4String("crystal") && step->GetTotalEnergyDeposit() > 0) || (volumeName == G4String("Patient") && processName == G4String("conv"))) {
      G4int eID = 0;
      const G4Event* evt = G4RunManager::GetRunManager()->GetCurrentEvent();
      if(evt) eID = evt->GetEventID();
      fFileManager->AddEnergyDepositionDetector(eID, track->GetPosition().x()/CLHEP::cm, track->GetPosition().y()/CLHEP::cm, track->GetPosition().z()/CLHEP::cm, step->GetTotalEnergyDeposit(), track->GetTrackID()+100*track->GetParentID(), processName);
    }

//#define RTPET_VERBOSE_PRINT 1
#ifdef RTPET_VERBOSE_PRINT
    G4cout << "\t";
    G4cout << "(" << track->GetTrackID() << "," << track->GetParentID() << "," << particleName << ")";
    G4cout << "(" << processName << ")";
    if(volumeName == G4String("crystal")) {
      icryst = theTouchable->GetCopyNumber();
      iring  = theTouchable->GetCopyNumber(1);
      idet   = theTouchable->GetCopyNumber(2);;
      G4cout << "(" << volumeName << ":" << idet << "," << iring << "," << icryst << ")";
    } else {
      G4cout << "(" << volumeName << ")";
    }
    G4cout << "(E:"
           << step->GetPreStepPoint()->GetKineticEnergy() << "->"
           << step->GetPostStepPoint()->GetKineticEnergy() << ")";
    G4cout << "(x=" << track->GetPosition().x()/CLHEP::cm
           << ",y=" << track->GetPosition().y()/CLHEP::cm
           << ",z=" << track->GetPosition().z()/CLHEP::cm
           << ",dE=" << step->GetPreStepPoint()->GetKineticEnergy() - step->GetPostStepPoint()->GetKineticEnergy()
           << ")";
    G4cout << "\n";
#endif
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
