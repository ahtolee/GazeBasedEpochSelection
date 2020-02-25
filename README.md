# GazeBasedEpochSelection
Gaze data -based epoch selection algorithm for eye tracker assisted visual evoked potential paradigm

Function GazeBasedOptimization.m
Used in steady-state VEP analysis to perform optimal epoch selection based on the eye tracking data (validity codes) collected during the steady-state
stimulus presentation. The output of the algorithm is the epoch division that maximized overall gaze quality for the given number of epochs.
Algorithm also takes into account predefined artefact events, as those segements are excluded from the epochs.

Attached file test_data.mat contains authentic test data from one recording in a form fully compatible with the script. Instructions for the use are provided in the header part of the m-file. 

Developed for the study:
"Use of complex visual stimuli allows controlled recruitment of cortical networks in infants"
Authors: Eero Ahtola1,2, Susanna Stjerna1, Anton Tokariev1,3, Sampsa Vanhatalo1,3

Eero Ahtola. 2020.
eero.ahtola@hus.fi

INPUTS VARIABLES:
EpochN_select     = Give the number epochs of that will be searched (target, e.g. 60)
EpochTime         = Duration of one epoch in sec (e.g. 1)
Reversals         = Matrix array of information of all reversal segments in the EEG data 
                 1st column -> Start index (sample) of each reversal segment
                 2nd column -> Duration (in samples) of each reversal segment
                 3rd column -> The sample where the phase reversal occured in the stimulus. Should be in the middle of the reversal segment.
                  4th column -> Block number of for which each reversal segment belongs to.
blockind          = Separator indices (start, end) for each stimulus block (uninterrupted stimulation sequence). In gaze data samples. 
ValidityLeftEye   = Validity codes from the eye tracker for left eye. Code 0 means that gaze is directed towards the stimulus.
ValidityRightEye	= Validity codes from the eye tracker for right eye. Code 0 means that gaze is directed towards the stimulus.
Artefacts         = Predefined reversal segments that are likely to contain EEG artefacts. These segments will be excluded from the final epoch selection.  
Fs                = Sampling frequency for EEG data (e.g. for samples in Reversal matrix)
FsGaze            = Sampling frequency for gaze data (for ValidityLeftEye and ValidityRightEye; probably 120Hz)

OUTPUT VARIABLES:
Epochs_optimized = Output matrix with information about the optimal epoch division. Size: EpochN_select x 4.
                 1st column -> Start index (sample) of each epoch (using given sampling rate of the EEG data)
                 2nd column -> End index (sample) of each epoch (using given sampling rate of the EEG data)
                 3rd column -> Block number of for which each epoch belongs to.
                 4th column -> Gaze Quality index (e.g mean validity between 0-1) for each epoch. Overall Gaze Quality = mean(orind(:,4))  
