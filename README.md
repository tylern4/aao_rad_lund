# aao_rad_lund

To build with makefile:
```
mkdir -p lib
mkdir -p bin
make
```

To build with cmake:
```
mkdir build
cd build
cmake ..
make
```

Make sure you have the maid tables:
```
tar -xvf parms.tar.gz
export CLAS_PARMS=${PWD}/parms
```

To run aao_rad from the build or bin folder:
```
./aao_rad_lund < test.inp
```

To run with command line arguments:
```
export PATH=${PATH}:${PWD}/bin
./aao_rad --trig 1000 --experiment rgb --q2min 0.25
```

