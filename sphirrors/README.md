SPHERICAL MIRROR FRACTAL GENERATOR

This program generates fractal images using raytracing applied to systems of spherical mirrors.

Run `run.jl` to use.

The program will automatically detect and, if necessary, download the required Julia packages from the Julia Standard Registry. This is done in `init.jl`.

The raytracing & reflection functions are defined in `sphirror_fcns.jl`. The parameters are defined and raytracing functions are called from `main.jl`.

By default, the program uses a square pyramid of spherical mirrors. The arrangement can be customized by editing the vector of sphere center coordinates, `C`, within `main.jl`. By default, the spheres all have the same radius.