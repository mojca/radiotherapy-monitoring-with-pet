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
/// \file RTPETDetectorConstruction.cc
/// \brief Implementation of the RTPETDetectorConstruction class

#include "RTPETDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4PSDoseDeposit.hh"
#include "G4VisAttributes.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETDetectorConstruction::RTPETDetectorConstruction()
: G4VUserDetectorConstruction(),
  fCheckOverlaps(true)
{
  DefineMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETDetectorConstruction::~RTPETDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETDetectorConstruction::DefineMaterials()
{
  G4NistManager* man = G4NistManager::Instance();

  G4bool isotopes = false;

  G4Element*  O = man->FindOrBuildElement("O" , isotopes);
  G4Element* Si = man->FindOrBuildElement("Si", isotopes);
  G4Element* Lu = man->FindOrBuildElement("Lu", isotopes);

  G4Material* LSO = new G4Material("Lu2SiO5", 7.4*g/cm3, 3);
  LSO->AddElement(Lu, 2);
  LSO->AddElement(Si, 1);
  LSO->AddElement(O , 5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* RTPETDetectorConstruction::Construct()
{
#define DETECTOR_MONOLITHIC 1
#define DETECTOR_NARROW_HOLE 1

  // Gamma detector Parameters
  //
  G4double cryst_dPhi;
  G4double cryst_dX, cryst_dY;
  G4double cryst_dZ = 2*cm;

  G4double ring_R1  = 50*cm;
  G4double ring_R2  = ring_R1 + cryst_dZ;
  G4double ring_dR2 = 0;

  G4double ring_h   = 40*cm;

  G4double ring_Phi1_deg, ring_Phi2_deg;
  G4double ring_Phi1,     ring_Phi2;
  G4double ring_dPhi;

// both need to be positive, sum has to be less than 180
#ifdef DETECTOR_NARROW_HOLE
  ring_Phi1_deg = 5.0;
  ring_Phi2_deg = 15.0;
  // ring_Phi2_deg = 5.0;
#else
  ring_Phi1_deg = 20.0;
  ring_Phi2_deg = 30.0;
#endif

  ring_Phi1 = (  ring_Phi1_deg/180.0)*pi;
  ring_Phi2 = (1-ring_Phi2_deg/180.0)*pi;
  ring_dPhi = ring_Phi2-ring_Phi1;

  G4int nb_rings, nb_cryst;

#ifdef DETECTOR_MONOLITHIC
  nb_rings = 1;
  nb_cryst = 1;
#else
  nb_rings = 10;
  nb_cryst = 20;

  cryst_dY = 2 * ring_R1 * std::tan(0.5 * ring_dPhi / nb_cryst);
  ring_dR2 = ring_R2 * (std::sqrt(1 + std::pow(0.5 * cryst_dY / ring_R2, 2)) - 1);
#endif
  cryst_dX   = ring_h    / nb_rings;
  cryst_dPhi = ring_dPhi / nb_cryst;

#ifndef DETECTOR_MONOLITHIC
  std::cout << "crystal dimensions: (x=" << cryst_dX/cm << "cm,y=" << cryst_dY/cm << "cm) "
            << "ring: " << (ring_R2 + ring_dR2)/cm << "cm\n";
#endif

  //
  G4NistManager* nist     = G4NistManager::Instance();
  G4Material* default_mat = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* cryst_mat   = nist->FindOrBuildMaterial("Lu2SiO5");

  //
  // World
  //
  G4double world_sizeXY = 2.4 * (ring_R2 + ring_dR2);
  G4double world_sizeY  = 2.0*m;
  G4double world_sizeZ  = 1.2 * std::max(ring_h, 10*cm); // 10cm: patient

  G4Box* solidWorld =
    new G4Box("World",                       // its name
       0.5*world_sizeXY, 0.5*world_sizeY, 0.5*world_sizeZ); // its size

  G4LogicalVolume* logicWorld =
    new G4LogicalVolume(solidWorld,          // its solid
                        default_mat,         // its material
                        "World");            // its name

  G4VPhysicalVolume* physWorld =
    new G4PVPlacement(0,                     // no rotation
                      G4ThreeVector(),       // at (0,0,0)
                      logicWorld,            // its logical volume
                      "World",               // its name
                      0,                     // its mother  volume
                      false,                 // no boolean operation
                      0,                     // copy number
                      fCheckOverlaps);       // checking overlaps

  //
  // single ring
  //
  G4Tubs* solidRing =
    new G4Tubs("Ring", ring_R1, ring_R2 + ring_dR2, 0.5*cryst_dX, 0, ring_dPhi);

  G4LogicalVolume* logicRing =
    new G4LogicalVolume(solidRing,           // its solid
                        default_mat,         // its material
                        "Ring");             // its name

  //
  // one half of the detector
  //
  G4Tubs* solidDetector =
    new G4Tubs("Detector", ring_R1, ring_R2 + ring_dR2, 0.5*ring_h, 0, ring_dPhi);

  G4LogicalVolume* logicDetector =
    new G4LogicalVolume(solidDetector,       // its solid
                        default_mat,         // its material
                        "Detector");         // its name

#ifdef DETECTOR_MONOLITHIC
  //
  // one crystal in form of the ring
  //
  G4Tubs* solidCryst =
    new G4Tubs("crystal", ring_R1, ring_R2, 0.5*ring_h, 0, ring_dPhi);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          // its solid
                        cryst_mat,           // its material
                        "CrystalLV");        // its name

  // place crystal in the ring
  new G4PVPlacement(0,                       // no rotation
                    G4ThreeVector(),         // at (0,0,0)
                    logicCryst,              // its logical volume
                    "crystal",               // its name
                    logicRing,               // its mother volume
                    false,                   // no boolean operation
                    0,                       // copy number
                    fCheckOverlaps);         // checking overlaps

#else
  //
  // define crystal (TODO)
  //
  //G4double gap = 0.5*mm;        //a gap for wrapping
  //G4double dX = cryst_dX - gap, dY = cryst_dY - gap;
  G4Box* solidCryst =
    new G4Box("crystal", cryst_dX/2, cryst_dY/2, cryst_dZ/2);

  G4LogicalVolume* logicCryst =
    new G4LogicalVolume(solidCryst,          // its solid
                        cryst_mat,           // its material
                        "CrystalLV");        // its name

  // place crystals in the ring
  //
  for (G4int icrys = 0; icrys < nb_cryst ; icrys++) {
    G4double phi = (icrys+0.5)*cryst_dPhi;
    G4RotationMatrix rotm  = G4RotationMatrix();
    rotm.rotateY(90*deg);
    rotm.rotateZ(phi);
    G4ThreeVector uz = G4ThreeVector(std::cos(phi), std::sin(phi),0.);
    G4ThreeVector position = (ring_R1+0.5*cryst_dZ)*uz;
    G4Transform3D transform = G4Transform3D(rotm,position);

    new G4PVPlacement(transform,             // rotation, position
                      logicCryst,            // its logical volume
                      "crystal",             // its name
                      logicRing,             // its mother volume
                      false,                 // no boolean operation
                      icrys,                 // copy number
                      fCheckOverlaps);       // checking overlaps
  }
#endif

  //
  // place rings within detector
  //
  G4double OG = -0.5*(ring_h + cryst_dX);
  for (G4int iring = 0; iring < nb_rings ; iring++) {
    OG += cryst_dX;
    // new G4PVPlacement(new G4RotationMatrix(ring_Phi1 - (iring%2)*(ring_Phi1+ring_Phi2),0,0), //no rotation
    new G4PVPlacement(0,                     // no rotation
                      G4ThreeVector(0,0,OG), // position
                      logicRing,             // its logical volume
                      "ring",                // its name
                      logicDetector,         // its mother volume
                      false,                 // no boolean operation
                      iring,                 // copy number
                      fCheckOverlaps);       // checking overlaps
  }

  //
  // place both parts of detector in world
  //
  for(G4int idet = 0; idet < 2; idet++) {
    new G4PVPlacement(new G4RotationMatrix(0.5*pi + ring_Phi1 - (idet%2)*(ring_Phi1+ring_Phi2),0,0), // rotation
                      G4ThreeVector(),       // at (0,0,0)
                      logicDetector,         // its logical volume
                      "Detector",            // its name
                      logicWorld,            // its mother volume
                      false,                 // no boolean operation
                      idet,                  // copy number
                      fCheckOverlaps);       // checking overlaps
  }

  //
  // patient
  //
  // G4double patient_radius = 8*cm;
  // G4double patient_dZ = 10*cm;
  G4Material* patient_mat      = nist->FindOrBuildMaterial("G4_BRAIN_ICRP");

// #define PATIENT_COMPLEX 1
// split into multiple parts for the sake of observing the beam energy spectrum
// #define PATIENT_SPLIT 1
#define PATIENT_BONE_LUNG 1 // first bone, then lung

#ifdef PATIENT_COMPLEX
  G4Material* patient_mat_bone = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
//G4Material* patient_mat_bone = nist->BuildMaterialWithNewDensity("Brain_Bone", "G4_BRAIN_ICRP", 1.850*g/cm3);
  // TODO: mix with air!!!
//G4Material* patient_mat_lung = nist->FindOrBuildMaterial("G4_LUNG_ICRP");
  G4Material* patient_mat_lung = nist->BuildMaterialWithNewDensity("Lung_0.3", "G4_LUNG_ICRP", 0.3*g/cm3);
//G4Material* patient_mat_lung = nist->BuildMaterialWithNewDensity("Brain_Lung", "G4_BRAIN_ICRP", 0.3*g/cm3);
#endif

  // G4Tubs* solidPatient =
  //   new G4Tubs("Patient", 0., patient_radius, 0.5*patient_dZ, 0., twopi);
  G4double patient_dZ = 20*cm;
  G4double patient_dX = 10*cm;
  G4double patient_dY = 10*cm;
  G4Box* solidPatient = new G4Box("Patient", patient_dY/2, patient_dZ/2, patient_dX/2);

  G4LogicalVolume* logicPatient =
    new G4LogicalVolume(solidPatient,        // its solid
                        patient_mat,         // its material
                        "PatientLV");        // its name

#ifdef PATIENT_COMPLEX
  G4double patient_dZ_bone = 4*cm;
  G4double patient_dZ_lung = 4*cm;

  G4Box* solidPatientBone = new G4Box("PatientBone", patient_dY/2, patient_dZ_bone/2, patient_dX/2);
  G4Box* solidPatientLung = new G4Box("PatientLung", patient_dY/2, patient_dZ_lung/2, patient_dX/2);

  G4LogicalVolume* logicPatientBone =
    new G4LogicalVolume(solidPatientBone,    // solid
                        patient_mat_bone,    // material
                        "PatientBoneLV");    // name
  G4LogicalVolume* logicPatientLung =
    new G4LogicalVolume(solidPatientLung,    // solid
                        patient_mat_lung,    // material
                        "PatientLungLV");    // name

#ifdef PATIENT_BONE_LUNG
  G4ThreeVector boneCenter(0,  5*cm, 0);
  G4ThreeVector lungCenter(0, -5*cm, 0);
#else
  G4ThreeVector boneCenter(0, -5*cm, 0);
  G4ThreeVector lungCenter(0,  5*cm, 0);
#endif
  new G4PVPlacement(0,                       // no rotation
                    boneCenter,              // at y=5cm
                    logicPatientBone,        // its logical volume
                 // "PatientBone",           // its name
                    "Patient",               // its name
                    logicPatient,            // its mother volume
                    false,
                    0,                       // copy number
                    fCheckOverlaps);
  new G4PVPlacement(0,                       // no rotation
                    lungCenter,              // at y=-5cm
                    logicPatientLung,        // its logical volume
                 // "PatientLung",           // its name
                    "Patient",               // its name
                    logicPatient,            // its mother volume
                    false,
                    0,                       // copy number
                    fCheckOverlaps);

#endif

  //
  // place patient in world
  //
  new G4PVPlacement(0,                       // no rotation
                    G4ThreeVector(),         // at (0,0,0)
                    logicPatient,            // its logical volume
                    "Patient",               // its name
                    logicWorld,              // its mother  volume
                    false,                   // no boolean operation
                    0,                       // copy number
                    fCheckOverlaps);         // checking overlaps

  //
  // Metal plate
  //
  G4double plate_thickness1 =  1   *mm;
  G4double plate_thickness2 = 14.25*mm;
  G4double plate_thickness3 =  3   *mm;
  G4double plate_size_XZ    =  5   *cm;
  G4double plate_thickness_sum = plate_thickness1 + plate_thickness2 + plate_thickness3;

  G4double plate_pos1 = 80*cm;
  G4double plate_pos2 = plate_pos1 - 0.5*(plate_thickness1+plate_thickness2);
  G4double plate_pos3 = plate_pos2 - 0.5*(plate_thickness2+plate_thickness3);

  G4Material* plate_mat1 = nist->FindOrBuildMaterial("G4_W");
  G4Material* plate_mat2 = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* plate_mat3 = nist->FindOrBuildMaterial("G4_Fe");

  G4Box* solidPlate1 = new G4Box("Plate1", plate_size_XZ/2, plate_thickness1/2, plate_size_XZ/2);
  G4Box* solidPlate2 = new G4Box("Plate2", plate_size_XZ/2, plate_thickness2/2, plate_size_XZ/2);
  G4Box* solidPlate3 = new G4Box("Plate3", plate_size_XZ/2, plate_thickness3/2, plate_size_XZ/2);

  G4LogicalVolume* logicPlate1 = new G4LogicalVolume(solidPlate1, plate_mat1, "Plate1LV");
  G4LogicalVolume* logicPlate2 = new G4LogicalVolume(solidPlate2, plate_mat2, "Plate2LV");
  G4LogicalVolume* logicPlate3 = new G4LogicalVolume(solidPlate3, plate_mat3, "Plate3LV");

  // Visualization attributes
  //
  logicRing->SetVisAttributes (G4VisAttributes::Invisible);
  logicDetector->SetVisAttributes (G4VisAttributes::Invisible);

  G4Colour lightYellow (1.0, 1.0, 0.5);
  G4VisAttributes* visCrystAtt = new G4VisAttributes(lightYellow);
  logicCryst->SetVisAttributes(visCrystAtt);

  G4VisAttributes* visPatientAtt = new G4VisAttributes(G4Colour::Blue());
  logicPatient->SetVisAttributes(visPatientAtt);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETDetectorConstruction::ConstructSDandField()
{
  G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  // declare crystal as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* cryst = new G4MultiFunctionalDetector("crystal");
  G4VPrimitiveScorer* primitiv1 = new G4PSEnergyDeposit("edep");
  cryst->RegisterPrimitive(primitiv1);
  SetSensitiveDetector("CrystalLV",cryst);

  // declare patient as a MultiFunctionalDetector scorer
  //
  G4MultiFunctionalDetector* patient = new G4MultiFunctionalDetector("patient");
  G4VPrimitiveScorer* primitiv2 = new G4PSDoseDeposit("dose");
  patient->RegisterPrimitive(primitiv2);
  SetSensitiveDetector("PatientLV",patient);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
