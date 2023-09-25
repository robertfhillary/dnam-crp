######################################
####### 8.0 OSCA EWAS MODELS #########
######################################

## Basic ## 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/
  osca_Linux \
--linear \
--befile bvals-norm20k-18413-831733 \
--pheno crp.phen \
--fast-linear \
--qcovar crp_quant.qcov \
--covar crp_fact.cov \
--out OSCA_Outputs/crp_standard_05042023 \
--methylation-beta


## FULL ## 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/
  osca_Linux \
--linear \
--befile bvals-norm20k-18413-831733 \
--pheno crp.phen \
--fast-linear \
--qcovar crp_quant_full.qcov \
--covar crp_fact_full.cov \
--out OSCA_Outputs/crp_standard_full \
--methylation-beta


## BASIC + BMI ONLY ## 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/
  osca_Linux \
--linear \
--befile bvals-norm20k-18413-831733 \
--pheno crp.phen \
--fast-linear \
--qcovar crp_quant_bmi.qcov \
--covar crp_fact.cov \
--out OSCA_Outputs/crp_standard_bmi \
--methylation-beta

## BASIC + SMK ONLY ## 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/
  osca_Linux \
--linear \
--befile bvals-norm20k-18413-831733 \
--pheno crp.phen \
--fast-linear \
--qcovar crp_quant_smk.qcov \
--covar crp_fact.cov \
--out OSCA_Outputs/crp_standard_smk \
--methylation-beta


## BASIC + BMI AND SMK ## 
cd /Cluster_Filespace/Marioni_Group/Rob/CRP/GS_EWAS/
  osca_Linux \
--linear \
--befile bvals-norm20k-18413-831733 \
--pheno crp.phen \
--fast-linear \
--qcovar crp_quant_bmismk.qcov \
--covar crp_fact.cov \
--out OSCA_Outputs/crp_standard_bmismk \
--methylation-beta

