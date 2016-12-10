///
/// \file RTPETFileManager.cc
/// \brief Implementation of the RTPETFileManager class

#include "RTPETFileManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include <cmath>
#include <iostream>
#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETFileManager::RTPETFileManager() {
  // TODO - this info should come from the gun, it shouldn't be hardcoded
  // TODO: don't forget to change settings in 'RTPETPrimaryGeneratorAction.cc'
  fEnergy    = 15;
  // fEnergy = 6;
  // 10^(fEventsExp)
  fEventsExp =  8;
  // fEventsExp =  7;
  // fEventsExp =  6;

  fNameSimulation = "brain";
  // fNameSimulation = "bonelung";
  // fNameSimulation = "lungbone";
  // fNameSimulation = "fakebonelung";

  fStoreAnnihilationYHisto       = true;
  fStoreDoseYHisto               = true;
  fStoreEnergyDepositionDetector = true;
  fStoreBeamEnergy               = true;

  fBinAnnihilationY_ymin  = -10.0;
  fBinAnnihilationY_ymax  =  10.0;
  if(fEventsExp <= 6) {
    fBinAnnihilationY_nbins =  40;
  } else if(fEventsExp >= 8) {
    fBinAnnihilationY_nbins = 200;
  } else {
    fBinAnnihilationY_nbins = 100;
  }

  fBinAnnihilationY_bins = new unsigned int[fBinAnnihilationY_nbins];
  for(unsigned int i=0; i<fBinAnnihilationY_nbins; i++) fBinAnnihilationY_bins[i] = 0;

  fBinDoseY_ymin  = -10.0;
  fBinDoseY_ymax  =  10.0;
  fBinDoseY_nbins =   200;

  fBinDoseY_bins   = new double[fBinDoseY_nbins];
  fBinEnergyY_bins = new double[fBinDoseY_nbins];
  for(unsigned int i=0; i<fBinDoseY_nbins; i++) {
    fBinDoseY_bins[i]   = 0;
    fBinEnergyY_bins[i] = 0;
  }

  fStoreAnnihilationXZSlice = true;
  if(fStoreAnnihilationXZSlice) {
    fBinAnnihilationXZHistoThick_ymin  = 6.0;
    fBinAnnihilationXZHistoThick_ymax  = 7.0;
    fBinAnnihilationXZHistoThick_nbins =  10;

    fBinAnnihilationXZHistoThin_ymin   = 6.4;
    fBinAnnihilationXZHistoThin_ymax   = 6.6;
    fBinAnnihilationXZHistoThin_nbins  =  50;
    // fBinAnnihilationXZHistoThick_bins  = new unsigned int[fBinAnnihilationXZHistoThick_nbins][fBinAnnihilationXZHistoThick_nbins];

    // TODO: memcpy
    for(unsigned int i=0; i<fBinAnnihilationXZHistoThick_nbins; i++) {
      for(unsigned int j=0; j<fBinAnnihilationXZHistoThick_nbins; j++) {
        fBinAnnihilationXZHistoThick_bins[i][j] = 0;
      }
    }
    for(unsigned int i=0; i<fBinAnnihilationXZHistoThin_nbins; i++) {
      for(unsigned int j=0; j<fBinAnnihilationXZHistoThin_nbins; j++) {
        fBinAnnihilationXZHistoThin_bins[i][j] = 0;
      }
    }
  }
  if(fStoreEnergyDepositionDetector) {
    // TODO
    std::string path_out("/tmp");
    char charEnergy[4];
    char charEventsExp[4];
    snprintf(charEnergy, 3, "%02d", fEnergy);
    snprintf(charEventsExp, 2, "%d", fEventsExp);

    std::string filename(path_out + "/data_" + charEnergy + "MV_1e" + charEventsExp + "_" + fNameSimulation + "_detector.bin");
    fileEnergyDepositionDetector = fopen(filename.c_str(), "wb");
  }

  if(fStoreBeamEnergy) {
    // std::string path_out("/tmp");
    // char charEnergy[4];
    // char charNum[4];
    // snprintf(charEnergy, 3, "%02d", fEnergy);

    for(unsigned int i = 0; i < 2; i++) {
      // snprintf(charNum, 3, "%02d", i);
      // std::string filename(path_out + "/data_"  + charEnergy + "MV_beam_energy_" + charNum + ".bin");
      // fileBeamEnergy[i] = fopen(filename.c_str(), "wb");
      for(unsigned int j = 0; j < fBinBeamEnergy_nbins; j++) {
        fBinBeamEnergy_bins[i][j] = 0;
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RTPETFileManager::~RTPETFileManager()
{
  std::string path_out("/tmp");
  char charEnergy[4];
  char charEventsExp[4];
  snprintf(charEnergy, 3, "%02d", fEnergy);
  snprintf(charEventsExp, 2, "%d", fEventsExp);

  /*
    fileHandle.open("...");
    if(fileHandle.is_open()) {
      fileHandle << "hello\n";
    }
  */

  if(fStoreAnnihilationYHisto) {
    // 20
    double yspan = fBinAnnihilationY_ymax-fBinAnnihilationY_ymin;
    // 0.5
    double dy = yspan/fBinAnnihilationY_nbins;

    std::string filename(path_out + "/data_" + charEnergy + "MV_1e" + charEventsExp + "_" + fNameSimulation + "_annihilation_y_hist.txt");
    FILE *f = fopen(filename.c_str(), "wt");
    fprintf(f, "# depth [cm]\tnumber of annihilations in the i-th bin [depth,depth+dy)\n");
    for(unsigned int i=0; i<fBinAnnihilationY_nbins; i++) {
      fprintf(f, "%4.1f\t%3d\n", i*dy, fBinAnnihilationY_bins[i]);
    }
    fprintf(f, "%4.1f\t  0\n", fBinAnnihilationY_nbins*dy);
    fclose(f);
  }
  if(fStoreDoseYHisto) {
    // 20
    double yspan = fBinDoseY_ymax-fBinDoseY_ymin;
    // 0.1
    double dy = yspan/fBinDoseY_nbins;

    dy = yspan/fBinDoseY_nbins;
    double totaldose = 0;
    std::string filename(path_out + "/data_" + charEnergy + "MV_1e" + charEventsExp + "_" + fNameSimulation + "_dose_y_hist.txt");
    FILE *f = fopen(filename.c_str(), "wt");
    fprintf(f, "# depth [cm]\tdeposited energy in the i-th bin [depth,depth+dy) [GeV/dm^3]\tdeposited dose [nGy]\n");
    for(unsigned int i=0; i<fBinDoseY_nbins; i++) {
      fprintf(f, "%4.1f\t%g\t%g\n", i*dy, fBinEnergyY_bins[i], fBinDoseY_bins[i]);
      totaldose += fBinDoseY_bins[i];
    }
    fprintf(f, "%4.1f\t0.0\t0.0\n", fBinDoseY_nbins*dy);
    fclose(f);
    // this doesn't work properly
    // G4cout << "total dose: " << totaldose/fBinDoseY_nbins * 0.16 << "\n";
  }
  if(fStoreAnnihilationXZSlice) {
    std::string filenameThick(path_out + "/data_" + charEnergy + "MV_1e" + charEventsExp + "_" + fNameSimulation + "_annihilation_xz_hist_1cm.txt");
    FILE *f = fopen(filenameThick.c_str(), "wt");
    fprintf(f, "# y=[%.0lf, %.0lf)\n", fBinAnnihilationXZHistoThick_ymin, fBinAnnihilationXZHistoThick_ymax);
    for(unsigned int i=0; i<fBinAnnihilationXZHistoThick_nbins; i++) {
      fprintf(f, "%d", fBinAnnihilationXZHistoThick_bins[i][0]);
      for(unsigned int j=1; j<fBinAnnihilationXZHistoThick_nbins; j++) {
        fprintf(f, "\t%d", fBinAnnihilationXZHistoThick_bins[i][j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
    std::string filenameThin(path_out + "/data_" + charEnergy + "MV_1e" + charEventsExp + "_" + fNameSimulation + "_annihilation_xz_hist_2mm.txt");
    f = fopen(filenameThin.c_str(), "wt");
    fprintf(f, "# y=[%.1lf, %.1lf)\n", fBinAnnihilationXZHistoThin_ymin, fBinAnnihilationXZHistoThin_ymax);
    G4cout << "looping " << fBinAnnihilationXZHistoThin_nbins << "\n";
    for(unsigned int i=0; i<fBinAnnihilationXZHistoThin_nbins; i++) {
      fprintf(f, "%d", fBinAnnihilationXZHistoThin_bins[i][0]);
      for(unsigned int j=1; j<fBinAnnihilationXZHistoThin_nbins; j++) {
        fprintf(f, "\t%d", fBinAnnihilationXZHistoThin_bins[i][j]);
      }
      fprintf(f, "\n");
    }
    fclose(f);
  }
  if(fStoreEnergyDepositionDetector) {
    if(fileEnergyDepositionDetector != NULL) {
      fclose(fileEnergyDepositionDetector);
    }
  }

  if(fStoreBeamEnergy) {
    char filename[200];
    snprintf(filename, 200, "/tmp/data_%02dMV_beam_energy.dat", fEnergy);
    // FILE *f = fopen(filename, "wb");
    // fwrite(fBinBeamEnergy_bins[i], sizeof(fBinBeamEnergy_bins[i][0]), fBinBeamEnergy_nbins, f);
    FILE *f = fopen(filename, "wt");

    for(unsigned int i = 0; i < fBinBeamEnergy_nbins; i++) {
      fprintf(f, "%u\t%u\n", fBinBeamEnergy_bins[0][i], fBinBeamEnergy_bins[1][i]);
    }
    fclose(f);
  }

  delete [] fBinAnnihilationY_bins;
  delete [] fBinDoseY_bins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RTPETFileManager::AddAnnihilation(G4double x, G4double y, G4double z)
{
  // check if it fail within the 4x4 region
  if(fabs(x) <= 2 && fabs(z) <= 2 && fabs(y) <= 10) {
    unsigned int bin_number = (unsigned int)(floor(fBinAnnihilationY_nbins*(10.0-y)/20.0));
    if(bin_number >= fBinAnnihilationY_nbins) {
      bin_number = fBinAnnihilationY_nbins-1;
    }
    fBinAnnihilationY_bins[bin_number]++;

    // int bin_number = (y-fBinAnnihilationY_ymin)/(fBinAnnihilationY_ymax-fBinAnnihilationY_ymin);
    // for 1D graph of annihilations along the y axes
    // if(fFileAnnihilationCoordinates.is_open()) {
    //   fFileAnnihilationCoordinates
    // }
  }
  // add a tiny layer
  if(fStoreAnnihilationXZSlice) {
    double xmin = -5, xmax = 5;
    double zmin = -5, zmax = 5;
    unsigned int xbin, zbin;
    // thick
    unsigned int n = fBinAnnihilationXZHistoThick_nbins;
    if(y >= fBinAnnihilationXZHistoThick_ymin && y < fBinAnnihilationXZHistoThick_ymax) {
      // G4cout << "event y=" << y << "\n";
      xbin = (unsigned int)(floor(n*(x-xmin)/(xmax-xmin)));
      zbin = (unsigned int)(floor(n*(z-zmin)/(zmax-zmin)));
      if(xbin >= n) xbin = n;
      if(zbin >= n) zbin = n;
      fBinAnnihilationXZHistoThick_bins[xbin][zbin]++;
      // G4cout << "x=" << xbin << " z=" << zbin << "\n";
    }
    // thin
    n = fBinAnnihilationXZHistoThin_nbins;
    if(y >= fBinAnnihilationXZHistoThin_ymin && y < fBinAnnihilationXZHistoThin_ymax) {
      // G4cout << "event y=" << y << "\n";
      xbin = (unsigned int)(floor(n*(x-xmin)/(xmax-xmin)));
      zbin = (unsigned int)(floor(n*(z-zmin)/(zmax-zmin)));
      if(xbin >= n) xbin = n;
      if(zbin >= n) zbin = n;
      fBinAnnihilationXZHistoThin_bins[xbin][zbin]++;
      G4cout << "x=" << xbin << " z=" << zbin << "\n";
    }
  }
}

void RTPETFileManager::AddEnergyDeposition(G4double x, G4double y, G4double z, G4double dE, G4double density)
{
  // check if it fail within the 4x4 region
  if(fabs(x) <= 2 && fabs(z) <= 2 && fabs(y) <= 10) {
    unsigned int bin_number = (unsigned int)(floor(fBinDoseY_nbins*(10.0-y)/20.0));

    double dy = (fBinDoseY_ymax-fBinDoseY_ymin)/fBinDoseY_nbins;
    double volume = 0.001*4*4*dy; // [dm^3]

    // 1 MeV = 0.1602176565 picoJoules
    double MeV_to_nJ = 0.1602176565e-3;
    double factor    = MeV_to_nJ * CLHEP::g / (volume * density * CLHEP::cm3);

    if(bin_number >= fBinDoseY_nbins) {
      bin_number = fBinDoseY_nbins-1;
    }
    fBinEnergyY_bins[bin_number] += 0.001 * dE / volume; // GeV/dm^3
    fBinDoseY_bins[bin_number]   += dE * factor;
  }
}

void RTPETFileManager::AddEnergyDepositionDetector(G4int eventNr, G4double x, G4double y, G4double z, G4double dE, G4int trackId, G4String processName)
{
  G4double angle = (fmod(atan2(y,x)/pi + 2.5, 2.0) - 1.0) * 180;

  float_t vars[3]; //= {(float_t)angle, (float_t)z, (float_t)dE};
  vars[0] = angle;
  vars[1] = z;
  vars[2] = dE;

  if(processName == G4String("conv")) {
    vars[0] = 0;
    vars[1] = 0;
    vars[2] = 0;
    // std::cout << eventNr << " PAIR!!!\n";
  }

  // TODO!!!
  // eventNr angle      z             dE
  // 1-N     [-180,180) [-20,20] [cm] (0,...) [MeV]
  if(fStoreEnergyDepositionDetector) {
    // std::cout << "edepcryst: " << eventNr << "\t" << x << "\t" << y << "\t" << z << "\t" << angle << "\t" << dE << "\t" << trackId << "\t" << processName << "\n";
    // printf("edepcryst: %6d %6.2f %6.2f %6.2f %7.2f %9.5f %d\n", eventNr, x, y, z, angle, dE, sizeof(vars[0]));
    fwrite(&eventNr, sizeof eventNr, 1, fileEnergyDepositionDetector);
  //fwrite(&angle,   sizeof angle,   1, fileEnergyDepositionDetector);
  //fwrite(&z,       sizeof z,       1, fileEnergyDepositionDetector);
  //fwrite(&dE,      sizeof dE,      1, fileEnergyDepositionDetector);
    fwrite(vars,     sizeof vars[0], 3, fileEnergyDepositionDetector);
  }
}

void RTPETFileManager::AddBeamEnergy(G4int variant, G4double energy)
{
  if(fStoreBeamEnergy) {
    if(variant == 1 || variant == 2) {
      G4double dE  = 1.0*fEnergy/fBinBeamEnergy_nbins;
      G4int    bin = int(floor(energy / dE));
      // G4cout << "** " << variant << " " << energy << " " << dE << " " << bin << "\n";
      if (bin >= 0 && bin < fBinBeamEnergy_nbins)
        fBinBeamEnergy_bins[variant-1][bin]++;
    }
  }
}
