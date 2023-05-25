FROM ubuntu:22.04

#WORKDIR /root

# workaround for tzdata asking for user input
# see https://dev.to/grigorkh/fix-tzdata-hangs-during-docker-image-build-4o9m for details
ENV TZ=Europe/Paris
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone


##################################
# INSTALLING OS AND BASIC DEPENDENCIES
##################################
# update ubuntu and install basic dependencies
RUN apt update && \
    apt install -y \
        apt-utils \
        build-essential \
        libboost-iostreams-dev \
        libtbb-dev \
        libblosc-dev \
        zlib1g-dev \
        cmake \
        git

##################################
# Installing OpenVDB
##################################

# clone git repo in home folder
RUN git clone https://github.com/AcademySoftwareFoundation/openvdb.git /openvdb && \
    cd /openvdb && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && make install

##################################
# Building the bridge
##################################

