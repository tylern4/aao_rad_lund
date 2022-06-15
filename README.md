# aao_rad_lund

To build:
```
mkdir -p lib
mkdir -p bin
make
```

Make sure you have the maid tables:
```
tar -xvf parms.tar.gz
export CLAS_PARMS=${PWD}/parms
```

To run aao_rad from the folder:
```
bin/aao_rad < test.inp
```


