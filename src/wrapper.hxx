#ifndef WRAPPER_H
#define WRAPPER_H

extern "C" {
void aao_rad_wrapper_(int *th_opt, int *flag_ehel, float *reg1, float *reg2,
                      float *reg3, float *reg4, int *npart, int *epirea,
                      float *mm_cut, float *t_targ, float *r_targ,
                      float *vertex_x, float *vertex_y, float *vz, float *ebeam,
                      float *q2_min, float *q2_max, float *ep_min,
                      float *ep_max, float *delta, int *nmax, float *fmcall,
                      float *sigr_max);
}

#endif