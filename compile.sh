#!/bin/bash

module purge
module load gcc


gcc -O2 gramSc.c -o gramSc -lblas -lm 
