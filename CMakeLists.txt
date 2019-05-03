cmake_minimum_required(VERSION 3.13)
project(aao_rad Fortran C)

set(CMAKE_C_FLAGS "-Ofast")
set(CMAKE_Fortran_FLAGS "-Ofast")

add_library(libs
        src/lenocc.c
        src/unixtime.c
        src/aao.f90
        src/abend.f90
        src/cgln_amps.f90
        src/daresbury.f90
        src/dsigma.f90
        src/fint.f90
        src/helicity_amps.f90
        src/interp.f90
        src/kerset.f90
        src/legendre.f90
        src/maid_lee.f90
        src/main.f90
        src/multipole_amps.f90
        src/read_sf_file.f90
        src/revinm.f90
        src/splie2.f90
        src/splin2.f90
        src/spline.f90
        src/splint.f90
        src/xsection.f90
        )

add_executable(aao_rad src/aao_rad.f90)
set_property(TARGET aao_rad PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(aao_rad libs)
