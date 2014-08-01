#!/bin/bash

inbix --numeric-file testdata10.num --pheno testdata10.pheno --1 --dcgain --out testdata10


inbix --numeric-file testdata10.num --pheno testdata10.pheno --1 --regain --out testdata10 --regain-matrix-transform abs

inbix --regain-file testdata10.block.regain --rank-by centrality --out testdata10 --rank-centrality-gamma 0.85

inbix --regain-file testdata10.block.regain --modularity --out testdata10

