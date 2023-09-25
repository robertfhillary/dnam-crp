##################################################
##########  9.0 OSCA REML MODELS #################
##################################################

#### Run on p17 due to version differences - p17 accepts combined multi-orm option #####

# Set working directory 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/
  
# Epigenetics alone (ORM)
osca_Linux --reml \
--orm GS_EWAS/bvals-norm20k-adj \
--pheno GS_EWAS/crp.phen \
--out GS_EWAS/orm_reml 
  
# Genetics alone (GRM)
gcta64 --reml \
--grm gs_grm \
--pheno GS_EWAS/crp.phen \
--thread-num 8 \
--out GS_EWAS/grm_reml

# Combined Genetics and Epigenetics (ORM+GRM)
osca_Linux --reml \
--multi-orm myorm.flist.txt \
--pheno GS_EWAS/crp.phen \
--out GS_EWAS/multireml 
