# scorpio

Scorpio MHD code with two fluid ambipolar diffusion solver.

Read the [Scorpio Wiki](https://github.com/SFG-CUHK/scorpio/wiki) for more details and test cases.



## Requirements

1. MPICH installed with enabling fortran compiler 
2. HDF5 installed with mpi & fortran enabled. 
3. FFTW installed with mpi & fortran enabled. 

## Get started

Download the source code using the standard git clone:

```
> git clone https://github.com/SFG-CUHK/scorpio.git
```

Go to the `./src` directory and type:

```
> make
```

to compile the executable file, a file named `Scorpio` will be created in the `./src` if the compilation success.

Run the program by the following line.

```
mpiexec -n <np> ./Scorpio
```

You can specific the number of processes to use using the `-n <np>` flag of `mpiexec`, see [here](https://www.mpich.org/static/docs/v3.1/www1/mpiexec.html) for more details.


After running the program, output file with format `g<test_case_id>_<number_of_output>.h5` in `.h5` format will be generated.


## Upgrade

In the working directory, type:
```
> git pull
```
