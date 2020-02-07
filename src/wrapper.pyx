# distutils: language = c++
from libcpp.string cimport string
from libcpp.vector cimport vector
cimport cython

import glob

cdef extern from "wrapper.hxx":
    void aao_rad_wrapper_(int *th_opt, int *flag_ehel, float *reg1, float *reg2,
                      float *reg3, float *reg4, int *npart, int *epirea,
                      float *mm_cut, float *t_targ, float *r_targ,
                      float *vertex_x, float *vertex_y, float *vz, float *ebeam,
                      float *q2_min, float *q2_max, float *ep_min,
                      float *ep_max, float *delta, int *nmax, float *fmcall,
                      float *sigr_max)



@cython.boundscheck(False)
@cython.cdivision(True)
@cython.wraparound(False)
@cython.infer_types(False)
def aao_generator(int th_opt=5, int flag_ehel=1, float reg1=0.2, float reg2=0.12,
                      float reg3=0.2, float reg4=0.2, int npart=3, int epirea=3,
                      float mm_cut=0.2, float t_targ=5.0, float r_targ=0.5,
                      float vertex_x=0.0, float vertex_y=0.0, float vz=0.0, float ebeam=4.8,
                      float q2_min=0.9, float q2_max=3.5, float ep_min=0.1,
                      float ep_max=4.25, float delta=0.0, int nmax=500, float fmcall=0.0,
                      float sigr_max=0.005):

                        cdef int *_th_opt=&th_opt
                        cdef int *_flag_ehel=&flag_ehel
                        cdef float *_reg1=&reg1
                        cdef float *_reg2=&reg2
                        cdef float *_reg3=&reg3
                        cdef float *_reg4=&reg4
                        cdef int *_npart=&npart
                        cdef int *_epirea=&epirea
                        cdef float *_mm_cut=&mm_cut
                        cdef float *_t_targ=&t_targ
                        cdef float *_r_targ=&r_targ
                        cdef float *_vertex_x=&vertex_x
                        cdef float *_vertex_y=&vertex_y
                        cdef float *_vz=&vz
                        cdef float *_ebeam=&ebeam
                        cdef float *_q2_min=&q2_min
                        cdef float *_q2_max=&q2_max
                        cdef float *_ep_min=&ep_min
                        cdef float *_ep_max=&ep_max
                        cdef float *_delta=&delta
                        cdef int *_nmax=&nmax
                        cdef float *_fmcall=&fmcall
                        cdef float *_sigr_max=&sigr_max

                        aao_rad_wrapper_(_th_opt, _flag_ehel, _reg1, _reg2, _reg3, _reg4, _npart, 
                            _epirea, _mm_cut, _t_targ, _r_targ, _vertex_x, _vertex_y, _vz,
                            _ebeam, _q2_min, _q2_max, _ep_min, _ep_max, _delta, _nmax, 
                            _fmcall, _sigr_max)
    


