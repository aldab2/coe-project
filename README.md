# coe-project

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal.c -o gauss_seidal
./gauss_seidal 400 0
```

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_wave.c -o gauss_seidal_wave
./gauss_seidal_wave 400 0
```

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_rb_double.c -o gauss_seidal_rb_double
./gauss_seidal_rb_double 400 0
```


```
mpicc gauss_seidal_rb_mpi.c -o gauss_seidal_rb_mpi
time mpirun -np 2 ./gauss_seidal_rb_mpi 8 0
```