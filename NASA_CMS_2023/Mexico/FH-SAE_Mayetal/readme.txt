Data and code for 'A spatially varying model for biomass density across the contiguous United States"

'UShex_SFH_RSE.R' contains the code, and is well-commented, but those unfamiliar with INLA
might want to refer to the free, online book https://becarioprecario.bitbucket.io/spde-gitbook/
in order to fully understand the code. The code, as is, fits the GEDI+TCC Spatial Fay-Herriot model, 
but predictors and spatial effects can be omitted/introduced according to the users desires. 
The later part of the script executes the 10-fold cross-validation study.

'GEDI_HEX.RData' contains the data necessary to fit the GEDI model.
'TCC_HEX.RData' contains the data necessary to fit the TCC model.
Both are needed to fit the GEDI+TCC model. The variables inside are described in the comments
of the R script.

CONUSbiohex2020 is a folder with the shapebundle from Menlove & Healey, 
and is also available at https://daac.ornl.gov/CMS/guides/FIA_Forest_Biomass_Estimates.html#datadescraccess.
It's not necessary for the analysis in the R script, but is useful for spatial plots,
and for comparing to (or using) the post-stratification estimates.
