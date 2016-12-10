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
/// \file RTPETActionInitialization.cc
/// \brief Implementation of the RTPETActionInitialization class

#include "RTPETActionInitialization.hh"
#include "RTPETPrimaryGeneratorAction.hh"
#include "RTPETRunAction.hh"
#include "RTPETStackingAction.hh"
#include "RTPETSteppingAction.hh"
#include "RTPETFileManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETActionInitialization::RTPETActionInitialization(
  RTPETDetectorConstruction* detConstruction,
  RTPETFileManager* fileManager)
: G4VUserActionInitialization(),
  fDetConstruction(detConstruction),
  fFileManager(fileManager)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETActionInitialization::~RTPETActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETActionInitialization::BuildForMaster() const
{
  SetUserAction(new RTPETRunAction);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETActionInitialization::Build() const
{
  SetUserAction(new RTPETPrimaryGeneratorAction);
  SetUserAction(new RTPETRunAction);
  SetUserAction(new RTPETStackingAction);
  SetUserAction(new RTPETSteppingAction(fDetConstruction,fFileManager));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......