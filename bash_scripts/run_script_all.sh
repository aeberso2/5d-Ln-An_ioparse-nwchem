#!/bin/bash

FUNCS1="svwn bp86 blyp pw91 pbe tpss m06l b3p86 x3lyp b97-1 b3lyp pbe0 mpw1k bhlyp tpssh m06 m06.2x m11 b2plyp"

FUNCS2="svwn bp86 blyp pw91 pbe tpss m06l b3p86 x3lyp b97-1 b3lyp pbe0 tpssh"

## WilsonCluster 
#for i in $FUNCS1; do run nwchem *.$i.rp.nw 4c_fast 4 $1; done

## CRUNTCH2
#for i in $FUNCS1; do nwchem *.$i.rp.nw -nc 16 -l; done

## ICER
#for i in $FUNCS1; do nwchem *.$i.rp.nw -n 4 -p 20 -m 60; done