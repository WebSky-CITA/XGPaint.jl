# Running on NERSC

Request a standard Haswell node (32 physical cores, 64 hyperthreads)
```
salloc -N 1 -C haswell -q interactive -t 03:00:00 --mem=0 --exclusive
```

Run Julia with 64 threads,
```
julia -t 64 cib_planck2013_chunked.jl
```
