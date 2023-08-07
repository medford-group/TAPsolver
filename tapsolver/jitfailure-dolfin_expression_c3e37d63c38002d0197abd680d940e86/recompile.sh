#!/bin/bash
# Execute this file to recompile locally
c++ -Wall -shared -fPIC -std=c++11 -O3 -fno-math-errno -fno-trapping-math -ffinite-math-only -I/Users/ayonge/anaconda3/envs/tapsolver-2/include -I/Users/ayonge/anaconda3/envs/tapsolver-2/include/eigen3 -I/Users/ayonge/anaconda3/envs/tapsolver-2/.cache/dijitso/include dolfin_expression_c3e37d63c38002d0197abd680d940e86.cpp -L/Users/ayonge/anaconda3/envs/tapsolver-2/lib -L/Users/ayonge/anaconda3/envs/tapsolver-2/Users/ayonge/anaconda3/envs/tapsolver-2/lib -L/Users/ayonge/anaconda3/envs/tapsolver-2/.cache/dijitso/lib -Wl,-rpath,/Users/ayonge/anaconda3/envs/tapsolver-2/.cache/dijitso/lib -lpmpi -lmpi -lmpicxx -lpetsc -lslepc -lm -ldl -lz -lpthread -lhdf5 -lboost_timer -ldolfin -Wl,-install_name,/Users/ayonge/anaconda3/envs/tapsolver-2/.cache/dijitso/lib/libdijitso-dolfin_expression_c3e37d63c38002d0197abd680d940e86.so -olibdijitso-dolfin_expression_c3e37d63c38002d0197abd680d940e86.so