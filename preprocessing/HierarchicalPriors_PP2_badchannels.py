# -*- coding: utf-8 -*-
"""
A first look at the data
"""

# load modules
import os.path as op
import os
import sys
import numpy as np
import glob

import mne
import matplotlib.pyplot as plt
import scipy

# set paths
rootdir = r'\\XXX\HierarchicalPriors'
code_path = os.path.join(rootdir,'code\preprocessing_MNEpython')
    
# load crosstalk and calibration files
crosstalk_file = os.path.join(code_path,'ct_sparse.fif')
cal_file = os.path.join(code_path,'sss_cal_3045_180914_upright.dat')
    
# loop over subjects
allsubs = ['MROS','MRSN','LANR','LCBR','MRDG','SMHS','CRSA','CIBN','FARN','RGBR','LRSU','NNGG','POML','LRVL','SLBG','AAMN','RTKM','SNMN','GIMR','ATFC','JNFN','RBTM','SLML','ANDA','LLMS','TTMD','ANDL','RTDC','PTCS','GGSI','AGCI','ANVC']# ['HLKR','MRBU','CAOM','DNCC','AKAl','LCTY','MRVN','DBBU','PTCR','NCGR','CICE']
for isub in range(len(allsubs)):
    datadir = os.path.join(rootdir,'data\MEG\ELEKTAoutput',allsubs[isub])  
    
    # loop over runs    
    allruns = glob.glob(os.path.join(datadir,'*.fif'))
    allchans2remove = np.zeros((len(allruns),2))
    for irun in range(len(allruns)):
       
        # load data
        file2load = allruns[irun]
        
        data = mne.io.read_raw_fif(file2load, allow_maxshield=True,preload=True,verbose=True)
        
        # # first check of power spectra
        # %matplotlib qt
        # data.plot_psd(fmax=60);
        
        # # check raw data
        # %matplotlib qt
        # data.plot(duration=10,title='tsss20');
                
        # if more than the max of 12 to-be-rejected sensors, rerun loop with lower limit
        chans2remove = []
        
        alllims = [3.5,3.6,3.7,3.8,3.9,4,4.5,5,5.5,6]
        ilim = 0
        while len(chans2remove) == 0 or len(chans2remove) > 12:
            
            # current limit, always start at 4 and only reduce if necessary
            current_limit = alllims[ilim]
            
            # identify noisy sensors to remove before Maxfiltering
            auto_noisy_chs, auto_flat_chs, auto_scores = mne.preprocessing.find_bad_channels_maxwell(
                data,
                limit=current_limit,
                cross_talk=crosstalk_file, 
                calibration=cal_file,
                return_scores=True, 
                verbose=True)
            
            chans2remove = auto_noisy_chs+auto_flat_chs
            ilim += 1
        
        del data # save storage space  
        
        # store relevant variables, chans2remove and current_limit, for all runs in single variable per subject
        allchans2remove[irun,0] = len(chans2remove)
        allchans2remove[irun,1] = current_limit

        # write chans2remove to .txt file so it can be read in directly by neuromag software for Maxfiltering
        with open(file2load[:len(file2load)-4]+'.badchannels.txt', 'w') as f:
            for chan in chans2remove:
                f.write(chan)
                f.write('\n')
        
    # save relevant variables
    from scipy.io import savemat
    savemat(os.path.join(datadir,'numchan2remove.mat'),mdict={'allchans2remove': allchans2remove})

# # add bad channel info to data 
# data1.info['bads'].extend(auto_noisy_chs + auto_flat_chs)

# # apply maxfilter
# data1_tsss = mne.preprocessing.maxwell_filter(
#     data1, 
#     st_duration=10,
#     cross_talk=crosstalk_file,
#     calibration=cal_file,
#     verbose=True)

# # check data after tsss
# %matplotlib qt
# data1_tsss.plot(duration=10,title='tsss_10badchan');
