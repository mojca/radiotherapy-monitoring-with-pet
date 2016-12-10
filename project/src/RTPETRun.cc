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
/// \file RTPETRun.cc
/// \brief Implementation of the RTPETRun class

#include "RTPETRun.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETRun::RTPETRun()
 : G4Run(),
   fCollID_cryst(-1),
   fCollID_patient(-1),
   fPrintModulo(10000),
   fGoodEvents(0),
   fSumDose(0.)
{ }


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETRun::~RTPETRun()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETRun::RecordEvent(const G4Event* event)
{
  if ( fCollID_cryst < 0 ) {
   fCollID_cryst
     = G4SDManager::GetSDMpointer()->GetCollectionID("crystal/edep");
   //G4cout << " fCollID_cryst: " << fCollID_cryst << G4endl;
  }

  if ( fCollID_patient < 0 ) {
   fCollID_patient
     = G4SDManager::GetSDMpointer()->GetCollectionID("patient/dose");
   //G4cout << " fCollID_patient: " << fCollID_patient << G4endl;
  }

  G4int evtNb = event->GetEventID();

  if ((evtNb+1)%fPrintModulo == 0) {
    G4cout << "\n---> end of event: " << (evtNb+1) << G4endl;
  }

  //Hits collections
  //
  G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  if(!HCE) return;

  //Energy in crystals : identify 'good events'
  //
  const G4double eThreshold = 500*keV;
  G4int nbOfFired = 0;

  G4THitsMap<G4double>* evtMap =
    static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_cryst));

  std::map<G4int,G4double*>::iterator itr;
  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    G4double edep = *(itr->second);
    if (edep > eThreshold) nbOfFired++;
    ///G4int copyNb  = (itr->first);
    ///G4cout << "\n  cryst" << copyNb << ": " << edep/keV << " keV ";
  }
  if (nbOfFired == 2) fGoodEvents++;

  //Dose deposit in patient
  //
  G4double dose = 0.;

  evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_patient));

  for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    ///G4int copyNb  = (itr->first);
    dose = *(itr->second);
    // G4cout << "\tdose: " << G4BestUnit(dose,"Dose") << G4endl;
  }
  fSumDose += dose;

  G4Run::RecordEvent(event);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETRun::Merge(const G4Run* aRun)
{
  const RTPETRun* localRun = static_cast<const RTPETRun*>(aRun);
  fGoodEvents += localRun->fGoodEvents;
  fSumDose    += localRun->fSumDose;

  G4Run::Merge(aRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
