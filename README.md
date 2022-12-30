# coe-project

```
/opt/homebrew/opt/llvm/bin/clang gauss_seidal.c -o gauss_seidal
gcc gauss_seidal.c -o gauss_seidal
./gauss_seidal 5
```


```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_wave.c -o gauss_seidal_wave
./gauss_seidal_wave 5
```

```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_rb_double.c -o gauss_seidal_rb_double
./gauss_seidal_rb_double 5
```
