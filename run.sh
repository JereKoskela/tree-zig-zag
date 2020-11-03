#!/bin/bash

mkdir ./Results

time ./Zig-zag/Infinite_sites/simulate ./Configs/zigzag-infinitesites-default.cfg > ./Results/zigzag-infinitesites-default.txt
time ./Zig-zag/Infinite_sites/simulate ./Configs/zigzag-infinitesites-long.cfg > ./Results/zigzag-infinitesites-long.txt
time ./Zig-zag/Infinite_sites/simulate ./Configs/zigzag-infinitesites-large.cfg > ./Results/zigzag-infinitesites-large.txt

time ./Zig-zag/Infinite_sites/simulate ./Configs/hybrid-infinitesites-default.cfg > ./Results/hybrid-infinitesites-default.txt
time ./Zig-zag/Infinite_sites/simulate ./Configs/hybrid-infinitesites-long.cfg > ./Results/hybrid-infinitesites-long.txt
time ./Zig-zag/Infinite_sites/simulate ./Configs/hybrid-infinitesites-large.cfg > ./Results/hybrid-infinitesites-large.txt

time ./Metropolis/Infinite_sites/simulate ./Configs/metropolis-infinitesites-default.cfg > ./Results/metropolis-infinitesites-default.txt
time ./Metropolis/Infinite_sites/simulate ./Configs/metropolis-infinitesites-long.cfg > ./Results/metropolis-infinitesites-long.txt
time ./Metropolis/Infinite_sites/simulate ./Configs/metropolis-infinitesites-large.cfg > ./Results/metropolis-infinitesites-large.txt

time ./Zig-zag/Finite_sites/simulate ./Configs/zigzag-finitesites-default.cfg > ./Results/zigzag-finitesites-default.txt
time ./Zig-zag/Finite_sites/simulate ./Configs/zigzag-finitesites-long.cfg > ./Results/zigzag-finitesites-long.txt
time ./Zig-zag/Finite_sites/simulate ./Configs/zigzag-finitesites-large.cfg > ./Results/zigzag-finitesites-large.txt

time ./Zig-zag/Finite_sites/simulate ./Configs/hybrid-finitesites-default.cfg > ./Results/hybrid-finitesites-default.txt
time ./Zig-zag/Finite_sites/simulate ./Configs/hybrid-finitesites-long.cfg > ./Results/hybrid-finitesites-long.txt
time ./Zig-zag/Finite_sites/simulate ./Configs/hybrid-finitesites-large.cfg > ./Results/hybrid-finitesites-large.txt

time ./Metropolis/Finite_sites/simulate ./Configs/metropolis-finitesites-default.cfg > ./Results/metropolis-finitesites-default.txt
time ./Metropolis/Finite_sites/simulate ./Configs/metropolis-finitesites-long.cfg > ./Results/metropolis-finitesites-long.txt
time ./Metropolis/Finite_sites/simulate ./Configs/metropolis-finitesites-large.cfg > ./Results/metropolis-finitesites-large.txt
