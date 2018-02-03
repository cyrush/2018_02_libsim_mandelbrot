#!/bin/bash
#########################################
# 
# file: go_mandelbrot_batch.sh
#
# script that runs mandelbrot libsim 
# example in batch mode, rendering and 
# exporting data without using 
#VisIt's UI
#########################################
# set up visit and lib paths
#export VISITHOME=path/to/visit/
#on osx, if using standard app bundle
export VISITHOME=path/to/VisIt.app/Contents/Resources/

# set visit version
export VISITVERSION=2.13.0
# set visit architecture
#
export VISITARCH=darwin-x86_64
#export VISITARCH=linux-x86_64

# for OSX you will need to set DYLD_LIB_PATH
export DYLD_LIBRARY_PATH=${VISITHOME}/${VISITVERSION}/${VISITARCH}/lib

# build the examples
make

# clean up old sim files
rm -rf ~/.visit/simulations/*
# clean up any render and export results
rm mandel*amr*.png
rm -rf mandel*amr*export*

# run the sim
./mandelbrot_batch -trace batch_trace.txt -dir ${VISITHOME}
