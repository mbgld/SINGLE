# set up
library(Seurat)
Split.outdir = "./Results/"

# FB: #FF6600
FB.sbj <- subset(ITG.sfj, idents = 'FB')

# EC: #CC9933
EC.sbj <- subset(ITG.sfj, idents = 'EC')

# MY: #009900
MY.sbj <- subset(ITG.sfj, idents = 'MY')

# MA: #999900
MA.sbj <- subset(ITG.sfj, idents = 'MA')

# TC: #00CCFF
TC.sbj <- subset(ITG.sfj, idents = 'TC')

# BC: ##66CC99
BC.sbj <- subset(ITG.sfj, idents = 'BC')

# EP: '#6699FF'
EP.sbj <- subset(ITG.sfj, idents = 'EP')

# CA: '#FF66CC'
CA.sfj <- subset(ITG.sfj, idents ='CA')

# EPICA: EP + CA 
EPICA.sbj <- subset(ITG.sfj, idents = c('CA', 'EP'))

################################################################################
### Save End-Products ##########################################################
################################################################################

save.image(paste0(Split.outdir, "{ }.RData"))
gc(); q()
