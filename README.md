# bsc-thesis
Code for my BSc thesis (A dynamic budget approach to identify a fast-slow life history continuum in microorganisms)

Order of running code (see code-file for specific instructions):
1. Calculating von Bertalanffy growth rates per species (id1-id14) using growth measurements (VONBERTALANFFY.R)
2. Running DEB-IPM (DEPIPMShrink_multiple_species.m)
3. Estimation of good & bad feeding levels using output DEB-IPM
4. Creating Excel-file with good feeding levels & Excel-file with bad feeding levels (Good_Bad_environment.m)
5. Running phylo PCA using parameters DEB-IPM for good feeding levels & bad feeding levels (PhyloPCA_Highfeedinglevel.R)
6. Running perturbation analysis (DEPIPMShrink_sens_multiple_species.m)
