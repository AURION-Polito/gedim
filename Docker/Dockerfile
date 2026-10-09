## Set Base Stage
FROM ubuntu:24.04 AS build_windows

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
apt-utils \
ssh \
software-properties-common \
apt-transport-https ca-certificates gnupg software-properties-common wget \
bash-completion

RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test \
&& apt update \
&& DEBIAN_FRONTEND=noninteractive apt install -y \
gcc-13 g++-13 gcc-13-base gfortran-13 libopenmpi-dev \
gcc-mingw-w64-x86-64 \
g++-mingw-w64-x86-64 \
gfortran-mingw-w64-x86-64 \
binutils-mingw-w64-x86-64 

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
&& update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 \
&& update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-13 100

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
&& update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 \
&& update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-13 100

RUN apt update && \
apt install -y cmake && \
apt install -y cmake-curses-gui && \
apt install -y git

WORKDIR /main_folder
