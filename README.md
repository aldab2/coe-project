# coe-project

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal.c -o gauss_seidal
./gauss_seidal 1400 0
```


```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_wave.c -o gauss_seidal_wave
./gauss_seidal_wave 600 0
```

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_rb_double.c -o gauss_seidal_rb_double
./gauss_seidal_rb_double 8
```


```
mpicc gauss_seidal_rb_mpi.c -o gauss_seidal_rb_mpi
time mpirun -np 2 ./gauss_seidal_rb_mpi 400
```

0.000783, 13.153779, 75.560532, 45.865013, 
53.276724, 21.895919, 4.704462, 67.886472, 
67.929641, 93.469290, 38.350208, 51.941637, 
83.096535, 3.457211, 5.346164, 52.970019, 

0.000783, 13.153779, 75.560532, 45.865013, 
53.276724, 39.665380, 55.123023, 67.886472, 
67.929641, 37.107985, 37.379705, 51.941637, 
83.096535, 3.457211, 5.346164, 52.970019,