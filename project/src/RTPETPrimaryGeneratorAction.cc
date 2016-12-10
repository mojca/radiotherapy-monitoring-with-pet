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
/// \file RTPETPrimaryGeneratorAction.cc
/// \brief Implementation of the RTPETPrimaryGeneratorAction class

#include "RTPETPrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ChargedGeantino.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

#include <valarray>
// #include <iostream>
#include <cassert>
#include <fstream>
#include <vector>
#include <iterator>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETPrimaryGeneratorAction::RTPETPrimaryGeneratorAction()
 : G4VUserPrimaryGeneratorAction(),
   fParticleGun(0)
{
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);

  // TODO: make this code more elegant
  // TODO: don't forget to change settings in 'RTPETFileManager.cc'
  // readEnergyListFromFile("rtpet-energy-spectrum-06MV.bin");
  readEnergyListFromFile("rtpet-energy-spectrum-15MV.bin");

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle
                    = particleTable->FindParticle("chargedgeantino");
  fParticleGun->SetParticleDefinition(particle);
  fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,0.));
  fParticleGun->SetParticleEnergy(1*eV);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETPrimaryGeneratorAction::~RTPETPrimaryGeneratorAction()
{
  delete fParticleGun;
  // TODO: empty energy_list
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETPrimaryGeneratorAction::readEnergyListFromFile(const char *filename)
{
  std::ifstream file_energy_list(filename, std::ios::binary);
  if(!file_energy_list)
  {
      std::cerr << "Cannot open the file '" << filename << "'" << std::endl;
      exit(1);
  }

  G4float f;
  while(file_energy_list.read(reinterpret_cast<char *>(&f), sizeof(f))) {
    energy_list.push_back(f);
  }

  energy_iterator = energy_list.begin();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

template<typename T>
unsigned my_binary_search(std::valarray<T> values, T value) {
  unsigned left=0, right=values.size()-1, middle;
  // special case for the last element which should't happen anyway
  // or rather: "if(value==1.0)"
  if(value >= values[right]) {
    return right;
  }
  while (true) {
    middle = (left + right)/2;
    if(value < values[left]) {
      return left;
    } else if(value >= values[right]) {
      return right+1;
    } else if(right - left == 1) {
      return right;
    } else {
      if(value < values[middle]) {
        right = middle;
      } else {
        left = middle;
      }
    }
  }
}

double get_x_from_linear_distribution(double y0, double y1, double p)
{
  assert((p >= 0) && (p<=1));

  if(p==0.0) {
    return 0.0;
  } else if(p==1.0) {
    return 1.0;
  } else if (y0 == y1) {
    return p;
  } else {
    double yy0 = 2*y0/(y0+y1);
    double yy1 = 2*y1/(y0+y1);

    double D_a = sqrt(yy0*yy0 + 2*(yy1-yy0)*p) / (yy1-yy0);
    double b_a = -yy0 / (yy1-yy0);

    if (b_a + D_a > 0) {
      return b_a + D_a;
    } else {
      return b_a - D_a;
    }
  }
}

G4double getEnergyFromPiecewiseLinearDistribution(double value)
{
  unsigned l = 10;
  // energy in MeV
  double x_init[] = {0.1, 0.2, 0.25, 0.5, 0.6, 0.95, 1.7, 4.9, 5.6, 5.8};
  // probability (in random units)
  double y_init[] = {100, 200, 400,  400, 300, 200,  100, 10,  3,   1};
  std::valarray<double> x(x_init, l), y(y_init, l);

  std::valarray<double> integral(l);
  double sum;
  assert((value >= 0.0) && (value <= 1.0));

  integral[0] = 0;
  for(unsigned i=1; i<l; i++) {
    // this can probably done with valarray operations
    integral[i] = 0.5*(y[i]+y[i-1])*(x[i]-x[i-1]) + integral[i-1];
  }
  sum = integral[l-1];
  for(unsigned i=0; i<l; i++) {
    integral[i] /= sum;
  }

  // the bin where the value falls
  unsigned bin = my_binary_search(integral, value)-1;
  // scaled probability (again between 0 and 1)
  double   p   = (value - integral[bin]) / (integral[bin+1] - integral[bin]);

  return (x[bin] + (x[bin+1] - x[bin])*get_x_from_linear_distribution(y[bin], y[bin+1], p)) * MeV;
}

G4double RTPETPrimaryGeneratorAction::getEnergyFromList()
{
  // TODO: if we reach the end, start from beginning
  if(energy_iterator == energy_list.end()) energy_iterator = energy_list.begin();

  return MeV * (*energy_iterator++);
  // std::cout << "Energy: " << energy_list.at(2) << "\n";
  // return energy_list.at(2) * MeV;
  // return 6.0*MeV;
}

void RTPETPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  G4ParticleDefinition* particle = fParticleGun->GetParticleDefinition();
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    G4String particleName;
    G4ParticleDefinition* gamma =
      G4ParticleTable::GetParticleTable()->FindParticle(particleName="gamma");
    G4ParticleDefinition* electron =
      G4ParticleTable::GetParticleTable()->FindParticle(particleName="e-");
    fParticleGun->SetParticleDefinition(gamma);
    // fParticleGun->SetParticleEnergy(2.0*MeV);
    // G4double energy = getEnergyFromPiecewiseLinearDistribution(G4UniformRand());
    // fParticleGun->SetParticleEnergy(energy);
    // fParticleGun->SetParticleDefinition(electron);
    // fParticleGun->SetParticleEnergy(6.0*MeV);
  }

  // G4double energy = getEnergyFromPiecewiseLinearDistribution(G4UniformRand());
  G4double energy = getEnergyFromList();
  fParticleGun->SetParticleEnergy(energy);

//#define RTPET_VERBOSE_PRINT 1
#ifdef RTPET_VERBOSE_PRINT
  G4cout << "** New particle with energy " << energy << " **\n";
#endif

  // randomized position
  //
  G4double x0  = 0*cm, y0  = 85*cm, z0  = 0*cm;
  G4double dx0 = 4*cm, dy0 =  0*cm, dz0 = 4*cm;
  x0 += dx0*(G4UniformRand()-0.5);
  z0 += dz0*(G4UniformRand()-0.5);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0,y0,z0));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0,-1,0));

  //create vertex
  //
  fParticleGun->GeneratePrimaryVertex(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
