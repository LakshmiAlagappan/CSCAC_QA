# CSCAC_QA
CSCAC is part of my industrial PhD thesis in NUS School of Computing. CSCAC (Class-specific Correction and Classification) is a framework to perform both standardization of spectra to mitigate batch effects as well as classification of samples accurately. Here, we show the application of CSCAC on a peanut oil-maize oil mixtures dataset. Here, we use the CSCAC to perform quality check to determine if the test spectra passes the quality check or not when compared against a target. 

1. Download the repository
2. Data folder: Contains NIR spectra belonging to train, transfer and test. (5 test batches)
3. helper.R: Contains all the helper and auxilary functions.
4. preprocessing.R: Contains the steps to preprocess the spectra and split the dataset into training, transfer and testing spectra 
5. Run main_flow.R 
