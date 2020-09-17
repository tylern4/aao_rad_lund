/*
  This header declares the C++ functions which can access the dsigma
  function supplied by Fortran.

  Argument descriptions:

  beam_energy:   Experiment beam energy in GeV
  W:             Center of mass energy in GeV
  Q2:            Minus the squared electron four-momentum transfer
  costheta:      cos(theta*) (* means center of mass)
  phi:           phi* (in degrees, also center of mass)
  helicity:      -1, 0, or 1 for electron beam helicity
  model_opt:     1=A0, 4=MAID98, 5=MAID2000
  channel_opt:   1=pi0, 2=pi-, 3=pi+
  resonance_opt: Selects the resonance, use 0 for full cross section

*/
extern "C" {
// Function which only returns the full cross section (not individual
// terms):
float maid_dsigma(float beam_energy, float W, float Q2, float costheta, float phi, int helicity, int model_opt,
                  int channel_opt, int resonance_opt);

// Function which accesses all terms of the cross section:
void maid_dsigma_all(  // Inputs:
    float beam_energy, float W, float Q2, float costheta, float phi, int helicity, int model_opt, int channel_opt,
    int resonance_opt,
    // Returns:
    float* sigma0, float* sigma_t, float* sigma_tt, float* sigma_l, float* sigma_lt, float* sigma_ltp, float* asym_p);
}
