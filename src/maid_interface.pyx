# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cython

import glob

cdef extern from "maid.h":
    float maid_dsigma(float beam_energy, float W, float Q2, float costheta, float phi, int helicity, int model_opt,
                  int channel_opt, int resonance_opt)



@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.infer_types(False)
def maid_2007_Npi(beam_energy, W, Q2, costheta, phi):
    helicity=0
    model_opt=5
    channel_opt=3
    resonance_opt=0
    return maid_dsigma(beam_energy, W, Q2, costheta, phi, helicity, model_opt,channel_opt, resonance_opt)

    


