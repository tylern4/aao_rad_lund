#include <exception>
#include <iostream>
#include "maid.h"

int main() {
  float W_bins[] = {1.1, 1.125, 1.15, 1.175, 1.2, 1.225, 1.25, 1.275, 1.3, 1.325, 1.35, 1.375,
                    1.4, 1.425, 1.45, 1.475, 1.5, 1.525, 1.55, 1.575, 1.6, 1.625, 1.65, 1.675,
                    1.7, 1.725, 1.75, 1.775, 1.8, 1.825, 1.85, 1.875, 1.9, 1.925};

  float Q2_bins[] = {0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  float cos_bins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  float phi_bins[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

  double dsigma = 0;

  int w_size = 32;
  int q2_size = 7;
  int c_size = 10;
  int p_size = 10;
  std::cout << "w,q2,cos_theta,phi,y\n";

  for (int w = 0; w < w_size; w++) {
    for (int q = 0; q < q2_size; q++) {
      for (int c = 0; c < c_size; c++) {
        for (int p = 0; p < p_size; p++) {
          std::cout << W_bins[w] << ",";
          std::cout << Q2_bins[q] << ",";
          std::cout << cos_bins[c] << ",";
          std::cout << phi_bins[p] << ",";
          dsigma = maid_dsigma(  // beam_energy
              4.81726,
              // W
              W_bins[w],
              // Q^2
              Q2_bins[q],
              // cos(theta*)
              cos_bins[c],
              // phi* (degrees)
              phi_bins[p],
              // helicity
              1,
              // model_opt
              5,
              // channel_opt
              3,
              // resonance_opt
              0);
          std::cout << dsigma << "\n";
        }
      }
    }
  }
}
