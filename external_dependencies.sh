#!/bin/bash


EXTERNAL_PATH="${1:-"/home/geoscore/Desktop/GEO++/public/external_build"}"

echo "-DCMAKE_PREFIX_PATH=${EXTERNAL_PATH}/Main_Install/eigen3;${EXTERNAL_PATH}/Main_Install/triangle;${EXTERNAL_PATH}/Main_Install/tetgen;${EXTERNAL_PATH}/Main_Install/vtk;${EXTERNAL_PATH}/Main_Install/googletest;${EXTERNAL_PATH}/Main_Install/lapack;${EXTERNAL_PATH}/Main_Install/metis;${EXTERNAL_PATH}/Main_Install/voro;/home/geoscore/Desktop/GEO++/public/gedim/release/GeDiM/GeDiM/lib/cmake/GeDiM"
