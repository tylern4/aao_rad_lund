cmake_minimum_required(VERSION 3.5)
project(aao_rad CXX C Fortran)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(CMAKE_C_FLAGS "-Ofast")
set(CMAKE_Fortran_FLAGS "-Ofast")

set(CMAKE_MODULE_PATH
   ${CMAKE_MODULE_PATH}
   ${PROJECT_SOURCE_DIR}/cmake)

include_directories("src")
add_library(aao_libs
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
        src/multipole_amps.f90
        src/read_sf_file.f90
        src/revinm.f90
        src/splie2.f90
        src/splin2.f90
        src/spline.f90
        src/splint.f90
        src/xsection.f90
        #src/wrapper_aao_rad.f90
        src/dvmpw.f90
        src/maidinterface.cpp
        src/maid.F
        )

add_executable(aao_rad src/aao_rad.f90)
set_property(TARGET aao_rad PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(aao_rad aao_libs)

add_executable(maid_csv src/maid_csv.cxx)
set_property(TARGET maid_csv PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries(maid_csv aao_libs)

find_package(Cython)
IF(Cython_FOUND)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  find_package(PythonExtensions REQUIRED)
  find_package(NumPy REQUIRED)

  include_directories(${PYTHON_INCLUDE_DIRS} ${PYTHON_NUMPY_INCLUDE_DIR} ${NumPy_INCLUDE_DIR})

  file(GLOB CYTHON_TARGETS src/*.pyx )
  foreach(pyx ${CYTHON_TARGETS} )
        get_filename_component(pyxname ${pyx} NAME)
        string( REPLACE ".pyx" "" pyxname ${pyxname} )
        add_cython_target(${pyxname} ${pyx} CXX PY3 OUTPUT_VAR _pyxname)
        add_library(${pyxname} MODULE ${_pyxname})
        python_extension_module(${pyxname})
        target_link_libraries(${pyxname} aao_libs ${PYTHON_LIBRARIES})
  endforeach( pyx ${CYTHON_TARGETS} )
ENDIF()

