# ParallelProgrammingProject
This project is about parallelize and accelerate an application described in the paper ["Tunable Approximation to Control Time-to-Solution in an HPC Molecular Docking Mini APP"](https://arxiv.org/abs/1901.06363). The aim of the application is to find the better pose of a molecule through all the of possible rotations of the fragments of the molecule.

## Technologies used
- C
- Cuda Toolkit V11
- OpenMP.

## Description of the files
All the .cu files in this repo contains the same application parallelized in different versions:
- `kernel_np.cu`: the serial implementation of the application. 
- `kernel_p_is_ligand_feasible.cu`: parallelized function is_ligan_feasible.
- `kernel_p_rotate_molecola.cu`: parallelized function rotate_molecola.
- `kernel_p_rotate_molecola_2.cu`: parallelized function rotate_molecola with different teqniques. 
- `kernel_p_measure_expansion.cu`: parallelized function measure_expansion.
- `kernel_p_place_best_angle.cu`: parallelized function place_best_angle.
- `kernel_p_place_best_angle_2.cu`: parallelized function place_best_angle with different teqniques. 
- `kernel_p_place_best_angle_3.cu`: parallelized function place_best_angle and concurrently processing the molecules using OMP threads. 

In the molecules folder there are the molecules in mol2 format. 

## How to compile
After cloning this repo, compile every file with `nvcc name_of_the_file.cu -o name_of_the_executable`. 
File `kernel_p_place_best_angle_3.cu` must be compiled with: `nvcc -Xcompiler -openmp kernel_p_place_best_angle_3.cu -o name_of_the_executable`
