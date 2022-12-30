# coe-project

```
/opt/homebrew/opt/llvm/bin/clang guess_seq_arr.c -o gauss-seidal
gcc gauss_seidal.c -o gauss_seidal
./gauss_seidal 2000
```


```
/opt/homebrew/opt/llvm/bin/clang -L/opt/homebrew/opt/libomp/lib -fopenmp gauss_seidal_wave.c -o gauss_seidal_wave
./gauss_seidal_wave 2000
```