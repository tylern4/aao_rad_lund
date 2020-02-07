from build.wrapper import aao_generator
import sys
import argparse
import timeit

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='aao_rad')

    parser.add_argument('-theory', action="store", default=5, type=int)
    parser.add_argument('-helicity', action="store", default=0, type=int)
    parser.add_argument('-reg1', action="store", default=0.20, type=float)
    parser.add_argument('-reg2', action="store", default=0.12, type=float)
    parser.add_argument('-reg3', action="store", default=0.20, type=float)
    parser.add_argument('-reg4', action="store", default=0.20, type=float)
    parser.add_argument('-npart', action="store", default=2, type=int)
    parser.add_argument('-epirea', action="store", default=3, type=int)
    parser.add_argument('-mm_cut', action="store", default=0.2, type=float)
    parser.add_argument('-t_targ', action="store", default=5.0, type=float)
    parser.add_argument('-r_targ', action="store", default=0.43, type=float)
    parser.add_argument('-vertex_x', action="store", default=0.0, type=float)
    parser.add_argument('-vertex_y', action="store", default=0.0, type=float)
    parser.add_argument('-vz', action="store", default=-3.0, type=float)
    parser.add_argument('-ebeam', action="store", default=10.646, type=float)
    parser.add_argument('-q2_min', action="store", default=1.0, type=float)
    parser.add_argument('-q2_max', action="store", default=10.0, type=float)
    parser.add_argument('-ep_min', action="store", default=1.4, type=float)
    parser.add_argument('-ep_max', action="store", default=10.6, type=float)
    parser.add_argument('-delta', action="store", default=0.005, type=float)
    parser.add_argument('-nmax', action="store", default=500, type=int)
    parser.add_argument('-fmcall', action="store", default=0.0, type=float)
    parser.add_argument('-sigr_max', action="store", default=0.005, type=float)
    args = parser.parse_args()

    print(
        f"aao_rad: running {args.nmax} events with a beam energy of {args.ebeam}\n")

    th_opt = args.theory
    flag_ehel = args.helicity
    reg1 = args.reg1
    reg2 = args.reg2
    reg3 = args.reg3
    reg4 = args.reg4
    npart = args.npart
    epirea = args.epirea
    mm_cut = args.mm_cut
    t_targ = args.t_targ
    r_targ = args.r_targ
    vertex_x = args.vertex_x
    vertex_y = args.vertex_y
    vz = args.vz
    ebeam = args.ebeam
    q2_min = args.q2_min
    q2_max = args.q2_max
    ep_min = args.ep_min
    ep_max = args.ep_max
    delta = args.delta
    nmax = args.nmax
    fmcall = args.fmcall
    sigr_max = args.sigr_max

    start = timeit.timeit()
    aao_generator(th_opt=th_opt,
                  flag_ehel=flag_ehel,
                  reg1=reg1,
                  reg2=reg2,
                  reg3=reg3,
                  reg4=reg4,
                  npart=npart,
                  epirea=epirea,
                  mm_cut=mm_cut,
                  t_targ=t_targ,
                  r_targ=r_targ,
                  vertex_x=vertex_x,
                  vertex_y=vertex_y,
                  vz=vz,
                  ebeam=ebeam,
                  q2_min=q2_min,
                  q2_max=q2_max,
                  ep_min=ep_min,
                  ep_max=ep_max,
                  delta=delta,
                  nmax=nmax,
                  fmcall=fmcall,
                  sigr_max=sigr_max)
