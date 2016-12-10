/// \file RTPETFileManager.hh
/// \brief Holds parameters and files for storing data during the run

// See also G4VScoreWriter in runAndEvent example

#ifndef RTPETFileManager_h
#define RTPETFileManager_h 1

#include "globals.hh"

#include <iostream>
#include <fstream>
#include <string.h>

class RTPETFileManager
{
public:
  RTPETFileManager();
  virtual ~RTPETFileManager();

public:
  G4int GetEnergy() const { return fEnergy; }

  // std::ofstream& GetFileAnnihilationCoordinates() { return fFileAnnihilationCoordinates; }
  // const std::ofstream& GetFileAnnihilationCoordinates() const { return fFileAnnihilationCoordinates; }
  void AddAnnihilation(G4double x, G4double y, G4double z);
  void AddEnergyDeposition(G4double x, G4double y, G4double z, G4double dE, G4double density);
  void AddEnergyDepositionDetector(G4int eventNr, G4double x, G4double y, G4double z, G4double dE, G4int trackId, G4String processName);
  void AddBeamEnergy(G4int variant, G4double energy);

private:
  G4int fEnergy;
  G4int fEventsExp;
  std::string fNameSimulation;

  // G4bool fStoreAnnihilationCoordinates;
  // std::ofstream fFileAnnihilationCoordinates;

  bool fStoreAnnihilationYHisto;
  bool fStoreDoseYHisto;
  bool fStoreAnnihilationXZSlice;
  bool fStoreEnergyDepositionDetector;
  bool fStoreBeamEnergy;

  double fBinAnnihilationY_ymin;
  double fBinAnnihilationY_ymax;

  double fBinDoseY_ymin;
  double fBinDoseY_ymax;

  double fBinAnnihilationXZHistoThick_ymin;
  double fBinAnnihilationXZHistoThick_ymax;
  double fBinAnnihilationXZHistoThin_ymin;
  double fBinAnnihilationXZHistoThin_ymax;

  unsigned int fBinAnnihilationY_nbins;
  unsigned int fBinDoseY_nbins;
  unsigned int fBinAnnihilationXZHistoThick_nbins;
  unsigned int fBinAnnihilationXZHistoThin_nbins;

  unsigned int *fBinAnnihilationY_bins;
  double       *fBinEnergyY_bins;
  double       *fBinDoseY_bins;

  static const int fBinBeamEnergy_nbins = 500;
  uint32_t         fBinBeamEnergy_bins[2][fBinBeamEnergy_nbins];

  unsigned int fBinAnnihilationXZHistoThick_bins[10][10]; // TODO
  unsigned int fBinAnnihilationXZHistoThin_bins[50][50]; // TODO

  FILE *fileEnergyDepositionDetector = NULL;
};

#endif
