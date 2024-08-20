####################### CODE DOCMENTATION FOR GLOBAL FOREST CLASSIFICATION ###########################
CONTACT: 
Dr. Neha Hunka
Assistant Research Professor
University of Maryland
nhunka@umd.edu
#####################################################################################################

The following process is set up for the classification of the world's forests into primary, young secondary and old 
secondary forests, as per the 2019 Refinement to the 2006 IPCC Guidelines for National Greenhouse Gas Inventories Volume 4 Agriculture, Forestry and Other Land Use, Table 4.7 for natural forests.

1. Various EO-derived and spatial datasets are downloaded from source (wget or curl commands)
2. Layers are spatially resampled and aligned to an approx. 30 m grid (GDAL commands)
3. A Boolean set of conditions is applied to layers to classify into forest statuses/conditions (AWS DPS algorithm)

All steps are reproducible for batch processing on the AWS DPS cloud-computing system that supports the NASA MAAP. 
For ease of use, step 1 and step 2 are broken down per 10 x 10 degree tile and described in the file 
NOTES_data_download_and_preprocessing.ipynb such that they are implementable on local machines using R / Python. 

Step 3 is provided as a DPS algorithm, which means that every 10 x 10 degree tile across the globe runs in parallel 
on AWS. For ease of understanding, it is recommended to start with the file FOREST_Classification/IPCC_GEDI_Table4.7.py. 
The Boolean combination used for the global forest classification is contained entirely within this file. It can be 
run from the command like and executed for a single tile if needed. 

#####################################################################################################