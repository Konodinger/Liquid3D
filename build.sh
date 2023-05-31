#!/bin/bash

# Define color variables
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
NC='\033[0m'
ORANGE='\033[38;5;172m'
BLUE='\033[34m'       
MAGENTA='\033[35m'       
CYAN='\033[36m'       
WHITE='\033[37m'       
BOLDBLACK='\033[1m\033[30m' 
BOLDRED='\033[1m\033[31m' 
BOLDGREEN='\033[1m\033[32m' 
BOLDYELLOW='\033[1m\033[33m' 
BOLDBLUE='\033[1m\033[34m' 
BOLDMAGENTA='\033[1m\033[35m' 
BOLDCYAN='\033[1m\033[36m' 
BOLDWHITE='\033[1m\033[37m' 

echo -e "${BOLDCYAN}Building 3D IISPH Simulator${NC}"
cd simulation/v2\ -\ IISPH\ 3D/
mkdir build
cmake -B build -S .
cmake --build build
echo -e "${BOLDGREEN}3D IISPH Simulator built successfuly${NC}"

cd ../../
echo -e "\n"

echo -e "${BOLDCYAN}Building OpenVDB Bridge${NC}"
cd openvdbCppBridge
mkdir build
cmake -B build -S .
cmake --build build
echo -e "${BOLDGREEN}OpenVDB Bridge built successfuly${NC}"

cd ../../
