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
./aao_rad < test.inp
```


