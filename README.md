# coe-project

```
/opt/homebrew/opt/llvm/bin/clang guess_seq_arr.c -o guess_seq_arr
gcc guess_seq_arr.c -o guess_seq_arr
./guess_seq_arr
```

```
/opt/homebrew/opt/llvm/bin/clang guess_seq_arr.c -o guess_seq_arr
gcc seq.c -o seq
./seq
```

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp par_seq.c -o par_seq
./par_seq
```


```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp par.c -o par
./par
```

// Run through gprof
```
gcc  -Wall -pg -fopenmp par.c -o par
./par
gprof par gmon.out > profile.txt
cat profile.txt
```