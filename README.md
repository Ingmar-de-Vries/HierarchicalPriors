UNDER CONSTRUCTION - UNDER CONSTRUCTION - UNDER CONSTRUCTION

# Hierarchical priors enable neural prediction of perceived biological motion

DOI code: TBD

This repository contains Matlab code accompanying the scientific manuscript (preprint) available at 

De Vries, I.E.J., De Lange, F.P., Wurm, M.F. Hierarchical priors enable neural prediction of perceived biological motion [https://doi.org/10.1038/s41467-023-39355-y](https://doi.org/10.1101/2025.10.09.681210)

If you use the data on OSF or the code here, please cite the above article.

The code in this repository allows for exact replication of the dRSA pipeline as presented in the article, but is also meant as inspiration for people interested in using dRSA to answer their own research questions. We're currently building a dRSA toolbox that will be shared on GitHub later 2026. In principle, dRSA can be applied to many different sensory modalities and contexts (e.g., naturalistic sound scenes, music, language), and on any signal with high enough temporal resolution (M/EEG, ECoG, eyetracking, etc.). Additionally, it should be straightforward to implement different dissimilarity measures for the neural and model RDMs, and a different similarity measure for the dRSA (i.e., here the principal component regression approach). We have found this similarity measure to be most effective for the current experiment, as tested with the simulations, but we have not extensively tested different dissimilarity measures.  

For details regarding this experiment, stimuli and analysis code, please see methods section of the article. Please contact me for any further questions at i.e.j.de.vries@gmail.com or ingmar.devries@unitn.it

Note that the larger data files belonging to this repository are stored on a public OSF repository (DOI: TBD; or look for Ingmar de Vries - HierarchicalPriors). The OSF repository includes: 
  - MEG data of 1 example participant to test this code. See more information in the "dataset_readme.txt" file included in this repository. 
  - The 18 model RDMs used in the reported study (9 for the normal and inverted conditions, 9 for the scrambled condition).

Note that this custom-written code uses several functions from the Brainstorm (tested version: 3) and Fieldtrip (tested version: 20191113) toolboxes, and was written and tested in Matlab 2023b. Additionally, the first 2 MEG preprocessing steps made use of MNE-python (for automatic bad channel detection using MNE python’s find_bad_channels_maxwell algorithm and for detecting the middle run (in terms of head location) to which to spatially realign the other runs to during MaxFiltering), after which Maxfiltering was applied using Neuromag’s MaxFilter implementation (version 2.2) of Signal Source Separation (SSS). 

The code in this GitHub repository is structured as follows:

  -	Experiment
    - In the “experiment” subdirectory, you will find the experiment script “HierarchicalPriors_MEGexperiment.m” 
    - You need Psychophysics Toolbox Version 3 (PTB-3) to run this experiment. 
    - This experiment can in principle be run as a behaviour-only experiment, include eyetracking, or include eyetracking and MEG. However, I have only tested the latest version of this experiment in the MEG lab at CIMeC, using Matlab 2012b. You might need to make minor adjustments for your setup. 
    - In the subdirectory “experiment/stimuli”, you will find the 14 unique 5-second-long ballet dancing videos used in the experiment, plus the corresponding and temporally aligned 3D kinematic marker locations at 100 Hz, stored in Matlab matrices. 
    - The experiment script makes use of the following helper scripts also present in the experiment directory:
      - “angle2pix.m” – transform degrees of visual angle to pixels on screen
      - “HierarchicalPriors_CreateCatchTrials.m” – create pool of catch trials that the experiment script randomly picks from on each run.
      - "HierarchicalPriors_CatchTrialPool.mat" - a file containing the catch trial info that results from running the "HierarchicalPriors_CreateCatchTrials.m" script. 
   
  - Analysis and plotting of behavioural results
    - In the "behaviour" subdirectory, you will find the following script:
      - "HierarchicalPriors_BehaviouralAnalysis.m" - behavioural analysis and plotting of Figure 1d.
      - plotSpread helper functions for nice scatter plots.

  -	Pre-processing of MEG and eyetracking data
    - Pre-processing of eyetracking data was done using custom written script "HierarchicalPriors_PP0_eyetracking_asc2ft.m", which takes raw Eyelink data in asc format as input, and gives pre-processed eyetracking data in Fieldtrip format as output. This is subsequently used to create an eyetracker RDM per individual subject (see RDM section below). The folder also contains a script called "HierarchicalPriors_eyetracking_missingsamples.m" that counts the total amount of missing samples in the eyetracker data to report this information in the manuscript.
    - The first 2 MEG pre-processing steps were done using custom MNE-python scripts located in the subdirectory "preprocessing":
      - "HierarchicalPriors_PP1_headPositionHistory.py" - finds the run with the middle head position to spatially realign the other runs to during Maxfiltering.
      - "HierarchicalPriors_PP2_badchannels.py" - automatic bad channel detection of up to a maximum of 12 bad channels, which were interpolated during Maxfiltering.
      - "ct_sparse.fif" - crosstalk file needed for the automatic bad channel detection and Maxfiltering.
      - "sss_cal_3045_180914_upright.dat" - calibration file needed for the automatic bad channel detection and Maxfiltering.
    - After this Maxfiltering was applied using Neuromag's MaxFilter implementation (version 2.2) of Signal Source Separation (SSS) provided by Elekta Neuromag (MEGIN).  
    - The rest of the MEG pre-processing was done using the Brainstorm toolbox version 3 using GUI operations, which were transformed into Matlab scripts where possible:
      - "HierarchicalPriors_PP3_CAT12segmentation.m" - segments MRI scans into cortical surface with 15000 sources using CAT12 in Brainstorm.
      - "HierarchicalPriors_PP4_importRAW.m" - import (or better: link to) raw MEG data, and refine MEG-MRI co-registration with extra head points. 
      - "HierarchicalPriors_PP5_addEvents.m" - read events from trigger channel and give appropriate names.
      - "HierarchicalPriors_PP6_checkVidOnset.m" - only sometimes necessary, i.e., sometimes triggers were erroneously stored double. If that's the case, this script helps finding those duplicates so they can be removed manually in Brainstorm GUI. But only happened in very rare cases. 
      - "HierarchicalPriors_PP7_filters.m" - notch filter, downsample, and create powerspectra for sanity check.
      - "HierarchicalPriors_PP8_ICA.m" - run ICA for ocular and cardiac artifacts, separately for magneto- and gradiometers.
      - "HierarchicalPriors_PP9_detectArtifacts.m" - automatic bad segment detection based on large amplitude at low frequency (1-7 Hz, i.e., blinks, movements, etc.), or high frequency (40-240 Hz; i.e., muscle activity) 
      - "HierarchicalPriors_PP10_epoch_singletrialDCcorrection.m" - epoch and single-trial baseline correction. Here there's also a visual inspection to check whether the automatically detected trials should indeed all be removed, e.g., sometimes a bad segment is in our padding windows, or might not be bad in our eyes. Those can be marked here.
      - "HierarchicalPriors_PP11_export2FT.m" - export from Brainstorm to Fieldtrip format.
      - "HierarchicalPriors_PP12_realign2photodiode.m" - realign single trials to photodiode. This script is called from PP11. 
      - "HierarchicalPriors_PP13_computeInversionKernel.m" - apply minimum norm estimation (MNE) and store resulting inversion kernel to transform sensor level data to source level data outside of Brainstorm (which I do in the main dynamic RSA analysis script). 

  -	Create model RDMs
    - The 18 model RDMs themselves can be found in the OSF repository. 
    - RDMs based on video data, kinematic marker data, and eyetracker data are created outside of the main "HierarchicalPriors_dRSA.m" script (see section "Run dRSA analysis" below), and if necessary up- or downsampled to 100 Hz. The pre-created RDMs therefore have size 14x14x500x500 (i.e., stim1 x stim2 x timestim1 x timestim2). This is done to save computation time in the dRSA pipeline, i.e., the RDMs only need to be loaded in, not computed each time the dRSA script is run. 
    - In the "modelRDMs" subdirectory, you'll find the following scripts:
      - "HierarchicalPriors_DynamicModelRDMs_eyeTracker.m" - create dynamic RDM of individual subject eyetracking data.
      - "HierarchicalPriors_DynamicModelRDMs_pixelwise.m" - create dynamic RDM of smoothed grayscale pixelwise luminance values. 
      - "HierarchicalPriors_video2vector.m" - create smoothed grayscale vector representation of videos. Called from "DynamicPredictions_DynamicModelRDMs_pixelwise.m".
      - "HierarchicalPriors_DynamicModelRDMs_opticalflow.m" - create dynamic RDM of optical flow vectors.
      - "HierarchicalPriors_video2opticalflow.m" - create optical flow vector representation of videos. Called from "DynamicPredictions_DynamicModelRDMs_opticalflow.m".
      - "HierarchicalPriors_DynamicModelRDMs_kinematic.m" - create 6 dynamic RDMs of kinematic marker data, i.e., view-dependent and view-invariant posture, motion and acceleration.
      - "procrustes_constrain_rotationZaxis_IdV.m" - modified version of Matlab's procrustes.m, which now constrains rotation to vertical (Z) axis, because that is how we define viewpoint invariant body posture, motion and acceleration. Note that my modified version is correct, but currently very time inefficient, effectively more than doubling the total computation time. This was a later modification and I'm sure I can find a much faster implementation. Feel free to have a look in the script and suggest a faster implementation! Called from "HierarchicalPriors_DynamicModelRDMs_kinematic.m".
      - "HierarchicalPriors_exampleFigureModels.m" - plots illustrations of the different models for a single frame of 2 videos. It was used for creating Figure S4 in the article. This script also contains information about where each of the 13 kinematic markers were located on the ballet dancer's body.

  - Run dRSA analysis, statistics, and plotting
    - In the "dynamicRSA" subdirectory, you'll find the following scripts:
      - "HierarchicalPriors_master_dRSA.m" - the main analysis pipeline from which all other functions are called.
      - "cluster_shell.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
      - "HierarchicalPriors_defineSourceROIs.m" - create ROIs based on (combinations of) parcels of HCP atlas.
      - "HierarchicalPriors_checkAtlases.m" - just sanity check that correct atlas and inversion kernel will be selected in main analysis
      - "HierarchicalPriors_dRSA.m" - main analysis script, which is called from "HierarchicalPriors_master_dRSA.m"
      - "regressionBorderPerModel_smRDM30msec.mat" - file containing regression borders used to regress out model itself to attenuate effects of model autocorrelation. These borders are determined by the simulations (see methods section in article and explanation in "HierarchicalPriors_dRSA.m" for details). 
      - "HierarchicalPriors_STATS_ERFdynamicRSA_ROIsource.m" - run statistics on ROI-based dRSA results, and compute peak latency and representational spread (RS) index. This function is called from main script "DynamicPredictions_pipeline.m". 
      - "modelautocorr_slopes.mat" - file containing dRSA curves resulting from PCR on simulated data. This is used to compute the representational spread (RS) index (see methods section in article and explanation in "DynamicPredictions_STATS_ERFdynamicRSA_ROIsource.m" for details).
      - "HierarchicalPriors_runFTstats.m" - shell around Fieldtrip functions for running cluster-based permutation tests on 2D dRSA matrix or on averaged dRSA lag-plot. This function is called from "HierarchicalPriors_STATS_ERFdynamicRSA_ROIsource.m". See scripts for details. 
      - "HierarchicalPriors_PLOTS_dRSA.m" - plot ROI-based results, in article: figure 2a, 3, and S1.
      - "brewermap.m" - creates nice colormaps that are colorblind friendly. Not my code, for all colormaps and source code see: https://colorbrewer2.org/
      - "boundedline.m" - creates nice shading around lines, e.g., with a measure of distribution across subjects (here standard error). Not my code, for source code see https://github.com/kakearney/boundedline-pkg

- Run dRSA simulations, and plotting
  - Note that the simulations are also called from the main "HierarchicalPriors_master_dRSA.m" script. In the "simulations" subdirectory, you'll find the following scripts:
    - "cluster_shell_simulations.m" - used for sending analysis as parallel jobs to a computing cluster (e.g., with different subjects and ROIs in parallel).
    - "DHierarchicalPriors_dRSA_simulations.m" - run simulations.
    - "HierarchicalPriors_PLOTS_dRSA_simulations.m" - plot simulations. 
