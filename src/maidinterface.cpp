#include <iostream>

// Declare the Fotran function:
extern "C" {
void maid_( // Inputs:
    float *beam_energy, float *W, float *Q2, float *costheta, float *phi,
    int *helicity, int *model_opt, int *channel_opt, int *resonance_opt,
    // Ouputs:
    float *sigma0, float *sigma_t, float *sigma_tt, float *sigma_l,
    float *sigma_lt, float *sigma_ltp, float *asym_p);

// Define the C++ primary interface function:
void maid_dsigma_all( // Inputs:
    float beam_energy, float W, float Q2, float costheta, float phi,
    int helicity, int model_opt, int channel_opt, int resonance_opt,
    // Returns:
    float *sigma0, float *sigma_t, float *sigma_tt, float *sigma_l,
    float *sigma_lt, float *sigma_ltp, float *asym_p) {
  // Just a wrapper for the Fortran function

  // Create temporary variables:
  float *tmpbeam_energy = new float;
  *tmpbeam_energy = beam_energy;
  float *tmpW = new float;
  *tmpW = W;
  float *tmpQ2 = new float;
  *tmpQ2 = Q2;
  float *tmpcostheta = new float;
  *tmpcostheta = costheta;
  float *tmpphi = new float;
  *tmpphi = phi;
  int *tmphelicity = new int;
  *tmphelicity = helicity;
  int *tmpmodel_opt = new int;
  *tmpmodel_opt = model_opt;
  int *tmpchannel_opt = new int;
  *tmpchannel_opt = channel_opt;
  int *tmpresonance_opt = new int;
  *tmpresonance_opt = resonance_opt;

  maid_(tmpbeam_energy, tmpW, tmpQ2, tmpcostheta, tmpphi, tmphelicity,
        tmpmodel_opt, tmpchannel_opt, tmpresonance_opt,

        sigma0, sigma_t, sigma_tt, sigma_l, sigma_lt, sigma_ltp, asym_p);
  delete tmpbeam_energy;
  delete tmpW;
  delete tmpQ2;
  delete tmpcostheta;
  delete tmpphi;
  delete tmphelicity;
  delete tmpmodel_opt;
  delete tmpchannel_opt;
  delete tmpresonance_opt;
  return;
}

// Define the simple interface dsigma function:
float maid_dsigma(float beam_energy, float W, float Q2, float costheta,
                  float phi, int helicity, int model_opt, int channel_opt,
                  int resonance_opt) {
  // Setup needed variables
  float result;
  float *sigma0 = new float;
  float *sigma_t = new float;
  float *sigma_tt = new float;
  float *sigma_l = new float;
  float *sigma_lt = new float;
  float *sigma_ltp = new float;
  float *asym_p = new float;
  // Retrieve Fortran results:
  maid_dsigma_all( // Inputs:
      beam_energy, W, Q2, costheta, phi, helicity, model_opt, channel_opt,
      resonance_opt,
      // Returns:
      sigma0, sigma_t, sigma_tt, sigma_l, sigma_lt, sigma_ltp, asym_p);
  result = *sigma0;
  // Cleanup memory:
  delete sigma0;
  delete sigma_t;
  delete sigma_tt;
  delete sigma_l;
  delete sigma_lt;
  delete sigma_ltp;
  delete asym_p;
  // Return the total cross section:
  return result;
}
}
