#!/bin/bash

module load apps/miniforge/24.7.1/bin
conda create --prefix=condaenv
conda activate /users/2850244m/condaenv
conda install conda-forge::gsl