#!/bin/bash
#########################################
# 
# file: go_mandelbrot.sh
#
# script that runs mandelbrot libsim 
# example which connects to VisIt's UI
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

# copy custom ui for sim into folder where visit can find it
cp mandelbrot.ui ~/.visit/ui

# run the sim
./mandelbrot -trace trace.txt -dir ${VISITHOME}
