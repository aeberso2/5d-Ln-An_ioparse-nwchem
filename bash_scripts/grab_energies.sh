#!/bin/bash

FUNCS1="svwn bp86 blyp pw91 pbe tpss m06l b3p86 x3lyp b97-1 b3lyp pbe0 mpw1k bhlyp tpssh m06 m06.2x"
FUNCS2="svwn bp86 blyp pw91 pbe tpss m06l b3p86 x3lyp b97-1 b3lyp pbe0 tpssh"
echo "DFT"
    for i in ${FUNCS2}; do
grep "Total DFT"      *.$i.[raDP][peKP].out | awk '{print$5}' | head -n1; done

#grep "Total DFT+MP2" *.b2plyp.*.out | awk '{print$4}'
echo "SODFT"
for i in ${FUNCS2}; do
    grep "Total DFT"      *.$i.[raDP][peKP].out | awk '{print$5}' | tail -n1; done

