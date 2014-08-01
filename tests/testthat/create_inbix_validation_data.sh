#!/bin/bash
#
# create_inbix_validation_data.sh - Bill White 

# dcGAIN
inbix \
--numeric-file testdata10.num \
--pheno testdata10.pheno --1 \
--dcgain \
--dcgain-abs \
--out testdata10

# dmGAIN
inbix \
--numeric-file testdata10.num \
--pheno testdata10.pheno --1 \
--dmgain --dmgain-abs \
--out testdata10

# reGAIN
inbix \
--numeric-file testdata10.num \
--pheno testdata10.pheno --1 \
--regain \
--regain-matrix-transform abs \
--out testdata10 

# reGAIN + SNPrank
inbix \
--regain-file testdata10.block.regain \
--rank-by centrality \
--out testdata10 

# reGAIN + modularity
inbix \
--regain-file testdata10.block.regain \
--modularity \
--out testdata10
