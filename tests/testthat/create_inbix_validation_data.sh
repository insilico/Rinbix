#!/bin/bash
#
# create_inbix_validation_data.sh - Bill White 

inbix_prefix="testdata10"
inbix_num="${inbix_prefix}.num"
inbix_pheno="${inbix_prefix}.pheno"

# dcGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--dcgain \
--dcgain-abs \
--out ${inbix_prefix}

# dmGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--dmgain --dmgain-abs \
--out ${inbix_prefix}

# ----------------------------------------------------------------------------
# four scenarios of binary parameters: abs/no abs, z-value stat or beta coeff
# ----------------------------------------------------------------------------

# abs/z-value
# reGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--regain \
--regain-matrix-transform abs \
--out ${inbix_prefix}-abs-zval

# reGAIN + SNPrank
inbix \
--regain-file ${inbix_prefix}-abs-zval.block.regain \
--rank-by centrality \
--out ${inbix_prefix}-abs-zval

# reGAIN + modularity
inbix \
--regain-file ${inbix_prefix}-abs-zval.block.regain \
--modularity \
--out ${inbix_prefix}-abs-zval

# ----------------------------------------------------------------------------
# no abs/z-value
# reGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--regain \
--out ${inbix_prefix}-noabs-zval

# reGAIN + SNPrank
inbix \
--regain-file ${inbix_prefix}-noabs-zval.block.regain \
--rank-by centrality \
--out ${inbix_prefix}-noabs-zval

# reGAIN + modularity
inbix \
--regain-file ${inbix_prefix}-noabs-zval.block.regain \
--modularity \
--out ${inbix_prefix}-noabs-zval

# ----------------------------------------------------------------------------
# abs/beta coeff
# reGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--regain \
--regain-matrix-transform abs \
--regain-use-beta-values \
--out ${inbix_prefix}-abs-beta

# reGAIN + SNPrank
inbix \
--regain-file ${inbix_prefix}-abs-beta.block.regain \
--rank-by centrality \
--out ${inbix_prefix}-abs-beta

# reGAIN + modularity
inbix \
--regain-file ${inbix_prefix}-abs-beta.block.regain \
--modularity \
--out ${inbix_prefix}-abs-beta

# ----------------------------------------------------------------------------
# no abs/beta coeff
# reGAIN
inbix \
--numeric-file $inbix_num \
--pheno $inbix_pheno --1 \
--regain \
--regain-use-beta-values \
--out ${inbix_prefix}-noabs-beta

# reGAIN + SNPrank
inbix \
--regain-file ${inbix_prefix}-noabs-beta.block.regain \
--rank-by centrality \
--out ${inbix_prefix}-noabs-beta

# reGAIN + modularity
inbix \
--regain-file ${inbix_prefix}-noabs-beta.block.regain \
--modularity \
--out ${inbix_prefix}-noabs-beta
