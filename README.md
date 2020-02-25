# GazeBasedEpochSelection
Gaze data -based epoch selection algorithm for eye tracker assisted visual evoked potential paradigm

Function GazeBasedOptimization.m
Used in steady-state VEP analysis to perform optimal epoch selection (for EEG) based on the eye tracking data (validity codes) collected during the steady-state stimulus presentation. The output of the algorithm is the epoch division that maximized overall gaze quality for the given number of epochs.

Algorithm also takes into account predefined EEG artefact events, as those segements are excluded from the epochs.

Attached file GazeBasedOptimization.m includes the code for the algorithm

Attached file test_data.mat contains authentic test data from one recording in a form fully compatible with the script. Instructions for the use are provided in the header part of the m-file. 

Developed for the study:
"Use of complex visual stimuli allows controlled recruitment of cortical networks in infants"
Authors: Eero Ahtola, Susanna Stjerna, Anton Tokariev, Sampsa Vanhatalo

Eero Ahtola. 2020.
eero.ahtola@hus.fi
