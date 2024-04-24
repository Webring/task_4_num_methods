#!/bin/bash
g++ cpp/gauss_on_eigen.cpp -O3 -o cpp/build/gauss
g++ cpp/QR_givens.cpp -O3 -o cpp/build/qrgivens
g++ cpp/svd.cpp -O3 -o cpp/build/svd


echo "compile done"
