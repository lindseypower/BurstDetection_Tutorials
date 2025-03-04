"""
Modified from: https://github.com/tbardouille/papto_camcan/blob/main/01a_PAPTO_find_events.py 
Original Authors: Brendan Brady, Timothy Bardouille
Modified by: Lindsey Power, March 2025
"""
# Import libraries
import os
import sys
import time
import mne
import math
import numpy as np
import scipy.signal as ss
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import logging
from spectralevents_functions import *
import multiprocessing as mp
import warnings
from scipy import stats
from fooof import FOOOF
from fooof import FOOOFGroup
import random
from scipy import signal

def energyvec(f,s,Fs,width):
    
	# Modified from Shin et al., eLife, 2017

    # Return a vector containing the energy as a
    # function of time for frequency f. The energy
    # is calculated using Morlet's wavelets. 
    # s : signal
    # Fs: sampling frequency
    # width : width of Morlet wavelet (>= 5 suggested).


    dt = 1/Fs
    sf = f/width
    st = 1/(2 * np.pi * sf)

    t= np.arange(-3.5*st, 3.5*st, dt)
    m = morlet(f, t, width)

    y = np.convolve(s, m)
    y = (dt * np.abs(y))**2
    lowerLimit = int(np.ceil(len(m)/2))
    upperLimit = int(len(y)-np.floor(len(m)/2)+1)
    y = y[lowerLimit:upperLimit]

    return y

def morlet(f,t,width):

	# Modified from Shin et al., eLife, 2017

    # Morlet's wavelet for frequency f and time t. 
    # The wavelet will be normalized so the total energy is 1.
    # width defines the ``width'' of the wavelet. 
    # A value >= 5 is suggested.
    #
    # Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

    sf = f/width
    st = 1/(2 * np.pi * sf)
    A = 1/np.sqrt((st/2 * np.sqrt(np.pi))) 
    y = A * np.exp(-t**2 / (2 * st**2)) * np.exp(1j * 2 * np.pi * f * t)

    return y


def interpolate_60Hz_notch(PSD, PSD_fVec):


	#linear interpolate over 60 Hz notch filter
	m = (PSD[0][61]-PSD[0][57])/(62-58)
	b = PSD[0][61] - m*62
	for i in range(57,62):
		PSD[0][i] = m*PSD_fVec[i]+b

	return PSD

# This one only operates on one sample
def TFR_via_morlet_wavelet(s, fVec, Fs):

    # Adapted from spectralevents_ts2tfr function from Shin et al (2017)

    width = 10
    s = s.T

    # obtain time vector (tVec) from timecourse (tVec starting with t=0s)
    numSamples = s.shape[1]
    tVec = np.arange(numSamples)/Fs

    # find number of frequencies for convolution
    numFrequencies = len(fVec)

    # generate TFR row by row
    TFR = []
    B = np.zeros((numFrequencies, numSamples))
    # Frequency loop
    for j in np.arange(numFrequencies):
        B[j,:] = energyvec(fVec[j], signal.detrend(s[0,:]), Fs, width)
    TFR.append(B)

    return TFR, tVec

def spectralevents_ts2tfr (S,fVec,Fs,width):
    # spectralevents_ts2tfr(S,fVec,Fs,width);
    #
    # Calculates the TFR (in spectral power) of a time-series waveform by 
    # convolving in the time-domain with a Morlet wavelet.                            
    #
    # Input
    # -----
    # S    : signals = time x Trials      
    # fVec    : frequencies over which to calculate TF energy        
    # Fs   : sampling frequency
    # width: number of cycles in wavelet (> 5 advisable)  
    #
    # Output
    # ------
    # t    : time
    # f    : frequency
    # B    : phase-locking factor = frequency x time
    #     
    # Adapted from Ole Jensen's traces2TFR in the 4Dtools toolbox.
    #
    # See also SPECTRALEVENTS, SPECTRALEVENTS_FIND, SPECTRALEVENTS_VIS.

    S = S.T
    numTrials = S.shape[0]
    numSamples = S.shape[1]
    numFrequencies = len(fVec)

    tVec = np.arange(numSamples)/Fs

    TFR = []
    # Trial Loop
    for i in np.arange(numTrials):
        B = np.zeros((numFrequencies, numSamples))
        # Frequency loop
        for j in np.arange(numFrequencies):
            B[j,:] = energyvec(fVec[j], signal.detrend(S[i,:]), Fs, width)
        TFR.append(B)

    TFR = np.asarray(TFR)

    return TFR, tVec, fVec

# Function to compare threshold choice for spectral event code 
def find_bestThreshold_papto(data, thresholds, Fs, width, fVec,fmin,fmax):
    columns = ['threshold', 'coef']
    df_plot = pd.DataFrame(columns=columns)

    # Calculate PSD
    PSD, PSD_fVec = mne.time_frequency.psd_array_welch(data, Fs, fmin=1, fmax=80, n_fft=1000, n_overlap=900)

    # FOOOF modeling
    fm = FOOOF(peak_width_limits=(2,10), max_n_peaks=4, aperiodic_mode='fixed', min_peak_height=0.05, peak_threshold=1.5)
    fm.fit(PSD_fVec, np.squeeze(PSD), [fmin,fmax])
    exponent = fm.get_params('aperiodic_params', 'exponent')
    offset = fm.get_params('aperiodic_params', 'offset')

    # Reshapes the data into 2 second epochs 
    chan_data = data[1000:-1000]
    remainder = len(chan_data)%2000 
    if remainder > 0:
        chan_data = chan_data[:-remainder]
    chan_data = chan_data.reshape(int(len(chan_data)/2000),2000)
    TFR, tVec, fVec = spectralevents_ts2tfr(chan_data.T, fVec, Fs, width) #Shape: [283,16,2000]

    # make a list of 2d arrays (one for each trial). i.e. this is a list of spectrograms for this subject
    trial_tfrs = []
    for trial in range(TFR.shape[0]-1): # Drop the last epoch because its skewing the correlations
        trial_tfrs.append(TFR[trial])

    # Calculate median power of recording for thresholding 
    nc = (10**offset)/(fVec**(exponent)) 

    # Loops through each potential threshold 
    for threshold in thresholds:
        avg_powers = []
        percents = []
        for spectrogram in trial_tfrs:	

            # Calculate average beta power across the spectrogram
            spectrogram_avg_beta_power = np.mean(spectrogram)
            avg_powers.append(spectrogram_avg_beta_power)

            # Check which values surpass X times median power threshold 
            # Note: median power is calculated separately for each frequency band
            spectrogram_thresholded = []
            for i in np.arange(0,16):
                freq_spectrogram_thresholded = np.where(spectrogram[i] > nc[i]*threshold, 1, 0)
                spectrogram_thresholded.append(freq_spectrogram_thresholded)
            spectrogram_thresholded = np.asarray(spectrogram_thresholded)

            # Calculates the percentage of pixels in this epoch for which the power surpasses the threshold
            count = np.count_nonzero(spectrogram_thresholded)
            total = np.size(spectrogram_thresholded)
            percent_of_pixels = count/total
            percents.append(percent_of_pixels)

        # calculate correlation coefficient between spectrogram_avg_beta_power and percent_of_pixels
        x = np.corrcoef(avg_powers, percents)[1,0]      
        values = [threshold, x]
        columns = ['threshold', 'coef']
        dictionary = dict(zip(columns, values))
        df_plot = df_plot.append(dictionary, ignore_index=True)
    
    return df_plot

def spectralevents_find (findMethod, thrFOM, tVec, fVec, TFR, classLabels, neighbourhood_size, threshold, Fs):
    # SPECTRALEVENTS_FIND Algorithm for finding and calculating spectral 
    #   events on a trial-by-trial basis of of a single subject/session. Uses 
    #   one of three methods before further analyzing and organizing event 
    #   features:
    #
    #   1) (Primary event detection method in Shin et al. eLife 2017): Find 
    #      spectral events by first retrieving all local maxima in 
    #      un-normalized TFR using imregionalmax, then selecting suprathreshold
    #      peaks within the frequency band of interest. This method allows for 
    #      multiple, overlapping events to occur in a given suprathreshold 
    #      region and does not guarantee the presence of within-band, 
    #      suprathreshold activity in any given trial will render an event.
    #
    # specEv_struct = spectralevents_find(findMethod,eventBand,thrFOM,tVec,fVec,TFR,classLabels)
    # 
    # Inputs:
    #   findMethod - integer value specifying which event-finding method 
    #       function to run. Note that the method specifies how much overlap 
    #       exists between events.
    #   eventBand - range of frequencies ([Fmin_event Fmax_event]; Hz) over 
    #       which above-threshold spectral power events are classified.
    #   thrFOM - factors of median threshold; positive real number used to
    #       threshold local maxima and classify events (see Shin et al. eLife 
    #       2017 for discussion concerning this value).
    #   tVec - time vector (s) over which the time-frequency response (TFR) is 
    #       calcuated.
    #   fVec - frequency vector (Hz) over which the time-frequency response 
    #       (TFR) is calcuated.
    #   TFR - time-frequency response (TFR) (trial-frequency-time) for a
    #       single subject/session.
    #   classLabels - numeric or logical 1-row array of trial classification 
    #       labels; associates each trial of the given subject/session to an 
    #       experimental condition/outcome/state (e.g., hit or miss, detect or 
    #       non-detect, attend-to or attend away).
    #
    # Outputs:
    #   specEv_struct - event feature structure with three main sub-structures:
    #       TrialSummary (trial-level features), Events (individual event 
    #       characteristics), and IEI (inter-event intervals from all trials 
    #       and those associated with only a given class label).
    #
    # See also SPECTRALEVENTS, SPECTRALEVENTS_FIND, SPECTRALEVENTS_TS2TFR, SPECTRALEVENTS_VIS.

    # Initialize general data parameters
    # Number of elements in discrete frequency spectrum
    flength = TFR.shape[1]
    # Number of point in time
    tlength = TFR.shape[2]
    # Number of trials
    numTrials = TFR.shape[0]
    classes = np.unique(classLabels)

    # Median power at each frequency across all trials
    TFRpermute = np.transpose(TFR, [1, 2, 0]) # freq x time x trial
    TFRreshape = np.reshape(TFRpermute, (flength, tlength*numTrials))
    medianPower = np.median(TFRreshape, axis=1)

    # Spectral event threshold for each frequency value
    eventThresholdByFrequency = thrFOM*medianPower

    # Validate consistency of parameter dimensions
    if flength != len(fVec):
        sys.exit('Mismatch in frequency dimensions!')
    if tlength != len(tVec): 
        sys.exit('Mismatch in time dimensions!')
    if numTrials != len(classLabels): 
        sys.exit('Mismatch in number of trials!')

    # Find spectral events using appropriate method
    #    Implementing one for now
    if findMethod == 1:
        spectralEvents = find_localmax_method_1(TFR, fVec, tVec, eventThresholdByFrequency, classLabels, medianPower, neighbourhood_size, threshold, Fs)

    return spectralEvents

def find_localmax_method_1(TFR, fVec, tVec, eventThresholdByFrequency, classLabels, medianPower, neighbourhood_size,threshold, Fs):
    # 1st event-finding method (primary event detection method in Shin et
    # al. eLife 2017): Find spectral events by first retrieving all local
    # maxima in un-normalized TFR using imregionalmax, then selecting
    # suprathreshold peaks within the frequency band of interest. This
    # method allows for multiple, overlapping events to occur in a given
    # suprathreshold region and does not guarantee the presence of
    # within-band, suprathreshold activity in any given trial will render
    # an event.

    # spectralEvents: 12 column matrix for storing local max event metrics:
    #        trial index,            hit/miss,         maxima frequency,
    #        lowerbound frequency,     upperbound frequency,
    #        frequency span,         maxima timing,     event onset timing,
    #        event offset timing,     event duration, maxima power,
    #        maxima/median power
    # Number of elements in discrete frequency spectrum
    flength = TFR.shape[1]
    # Number of point in time
    tlength = TFR.shape[2]
    # Number of trials
    numTrials = TFR.shape[0]

    spectralEvents = []

    # Retrieve all local maxima in TFR using python equivalent of imregionalmax
    for ti in range(numTrials):

        # Get TFR data for this trial [frequency x time]
        thisTFR = TFR[ti, :, :]

        # Find local maxima in the TFR data
        data = thisTFR
        data_max = filters.maximum_filter(data, neighbourhood_size)
        maxima = (data == data_max)
        data_min = filters.minimum_filter(data, neighbourhood_size)
        diff = ((data_max - data_min) > threshold)
        maxima[diff == 0] = 0
        labeled, num_objects = ndimage.label(maxima)
        xy = np.array(ndimage.center_of_mass(data, labeled, range(1, num_objects + 1)))

        numPeaks = len(xy)

        peakF = []
        peakT = []
        peakPower = []
        for thisXY in xy:
            peakF.append(int(thisXY[0]))
            peakT.append(int(thisXY[1]))
            peakPower.append(thisTFR[peakF[-1], peakT[-1]])

        # Find local maxima lowerbound, upperbound, and full width at half max
        #    for both frequency and time
        Ffwhm = []
        Tfwhm = []
        for lmi in range(numPeaks):
            thisPeakF = peakF[lmi]
            thisPeakT = peakT[lmi]
            thisPeakPower = peakPower[lmi]

            # Indices of TFR frequencies < half max power at the time of a given local peak
            TFRFrequencies = thisTFR[:, thisPeakT]
            lowerInd, upperInd, FWHM = fwhm_lower_upper_bound1(TFRFrequencies,
                                                               thisPeakF, thisPeakPower)
            lowerEdgeFreq = fVec[lowerInd]
            upperEdgeFreq = fVec[upperInd]
            FWHMFreq = FWHM

            # Indices of TFR times < half max power at the frequency of a given local peak
            TFRTimes = thisTFR[thisPeakF, :]
            lowerInd, upperInd, FWHM = fwhm_lower_upper_bound1(TFRTimes,
                                                               thisPeakT, thisPeakPower)
            lowerEdgeTime = tVec[lowerInd]
            upperEdgeTime = tVec[upperInd]
            FWHMTime = FWHM / Fs

            # Put peak characteristics to a dictionary
            #        trial index,            hit/miss,         maxima frequency,
            #        lowerbound frequency,     upperbound frequency,
            #        frequency span,         maxima timing,     event onset timing,
            #        event offset timing,     event duration, maxima power,
            #        maxima/median power
            peakParameters = {
                'Trial': ti,
                'Hit/Miss': classLabels[ti],
                'Peak Frequency': fVec[thisPeakF],
                'Lower Frequency Bound': lowerEdgeFreq,
                'Upper Frequency Bound': upperEdgeFreq,
                'Frequency Span': FWHMFreq,
                'Peak Time': tVec[thisPeakT],
                'Event Onset Time': lowerEdgeTime,
                'Event Offset Time': upperEdgeTime,
                'Event Duration': FWHMTime,
                'Peak Power': thisPeakPower,
                'Normalized Peak Power': thisPeakPower / medianPower[thisPeakF],
                'Outlier Event': thisPeakPower > eventThresholdByFrequency[thisPeakF]
            }

            # Build a list of dictionaries
            spectralEvents.append(peakParameters)

    return spectralEvents
