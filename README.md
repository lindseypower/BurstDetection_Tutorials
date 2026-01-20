# BurstDetection_Tutorials

This repository provides a set of burst detection tutorials to accompany our comprehensive methodological review entitlted: *A Neuroscientist's Guide to Neural Burst Detection*. 

The *Tutorials* folder contains step-by-step worked examples demonstrating how seven common burst detection methods can be applied to single or multi-channel neurophysiological data. Tutorials are available for the following methods: 

 1. Amplitude Thresholding*
 2. Periodic-Aperiodic Parameterization of Transient Oscillations (PAPTO)*
 3. Extended Better Oscillation Detection (eBOSC)
 4. Cycle-by-Cycle Oscillation Detection
 5. Sliding Window Matching*
 6. Time-delay embedding Hidden Markov Modelling (tde-HMM)
 7. Convoutional Dictionary Learning (CDL)*

*In some cases, additional utility scripts are required to run the burst detection method. In these cases, the required utility scripts are indicated at the start of the tutorial and the scripts are available for download in the *Utils* folder. Ensure that the required utility scripts are saved in the same folder as the associated tutorial before attempting to run the tutorial. 

The following python package versions were used to generate this set of tutorials. All packages can be installed using *pip install*. Please refer to the import list at the start of each tutorial for information on the required python packages for a given tutorial.

* numpy==1.26.4
* matplotlib==3.10.1
* pandas==2.2.3
* scipy==1.15.2
* sklearn==1.6.1
* seaborn==0.13.2
* joblib==1.4.2
* mne==1.8.0
* alphacsc==0.4.1
* bycycle==1.1.0
* osl_dynamics==2.1.0
* ebosc==0.10.dev0

In addition, the *Tutorials* folder contains a preprocessing tutorial which walks the user through the necessary steps to download and preprocess a sample of BIDS-formatted MEG data and extract clean time series for use in burst detection. This tutorial uses the *Brainstorm* software available at https://neuroimage.usc.edu/bst/download.php . 

Sample data used in this tutorial is taken from the OMEGA database (https://doi.org/10.1016/j.neuroimage.2015.04.028) and can be downloaded from OpenNeuro (https://openneuro.org/datasets/ds000247/versions/1.0.2). Specific information on how to download and process this data can be found in the Brainstorm Preprocessing Tutorial. For convenience, the extracted source time series used in tutorials 1-6 is included above in the *Data* folder. Tutorial 7 can be completed directly on the raw BIDS-formatted data obtaind from OpenNeuro. 

Note that this resource is a compilation of methods previously published in the scientific literature and therefore is largely *not* original work. Therefore, when using these methods, be sure to cite the appropriate original source. Links to documentation and publications from the original developers are embedded in the individual tutorials. More information on specific methods can also be found in our associated review article (*link here*) 
