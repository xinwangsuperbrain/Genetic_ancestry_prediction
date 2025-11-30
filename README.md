# Introduction
These codes are for the paper:   Explicitly modeling genetic ancestry to improve polygenic prediction accuracy for height in a large, admixed cohort of US Latinos: Findings from HCHS/SOL

# Contents
* TrainingPGS_Using_YengoGWAS_includingUKB_inUKBEUR.m

This script trains a polygenic score (PGS) for height in UK Biobank participants of European ancestry. The model uses SNPs identified in the Yengo et al. (2022) height GWAS, which includes European-ancestry samples with UK Biobank participants.

* SOLINCA_Height_Prediction_using_PGS_training_inUKBEUR_withYengoGWAS_includingUKB.m

This script applied the height PGS trained from TrainingPGS_Using_YengoGWAS_includingUKB_inUKBEUR.m and incorporating Principal Components (PCs) to HSHC/SOL participants 

