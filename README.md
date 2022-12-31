# Gauss Sidel Using Red Black (OpenMP and MPI)

```
icx -fopenmp  gs_colored.c -o gs_colored
./gs_colored 
```

```
icx -fopenmp  gs_colored_seq.c -o gs_colored
./gs_colored <N> <NUM_THREADS>
```

```
mpicc gs_colored_mpi.c -o gs_colored-mpi./par_seq
mpirun -np <N> ./gs_colored-mpi <NUM_THREADS>
```

// Run through gprof
```
gcc  -Wall -pg -fopenmp par.c -o par
./par
gprof par gmon.out > profile.txt
cat profile.txt
```