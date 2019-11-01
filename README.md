# HeALTH_Tracker

Behavioral Analysis Code for HeALTH System
HeALTH_Tracker
------------------------------------------

Included are the MATLAB files to analyze the behavioral data from the HeALTH platform. After selecting the folder containing the videos of interest, the code will extract out lifespan data and behavioral features from each selected worm of interest within a device longitudinally.
Additionally, after analyzing each individual population, user can then selectively compile trials together (based on experimental conditions) and perform additional, downstream behavioral analysis as needed. 

Note: alphanumericSortFiles created by Stephen Cobeldick (Copyright (c) 2018)

Last Updated: 10/9/19


Files
-------------------------------------------
> generalManager.m
The main file handler ('generalManager.m') calls all of the functions required to obtain lifespan data and to extract behavioral metrics from the videos. 
Depending on user needs, certain sections can be commented out as needed. The user will select the video dataset to analyze from here.

> getChamberLocations
This folder/subfiles identifies worms to analyze in the device and gets chamber center locations for all files of the device in the chosen folder via cross-correlation to perform automated chamber identification.
It contains roomLocator.m (the 'main' function within the folder), createChamberMask.m, getChamberCenters.m, and checkCenters.m
MAIN OUTPUTS: 
- seasonalSnapshots (w/ chamber locations, cropped filtered images for each chamber)
- combCenters
- chamberIDs
- refChamber

> getPixelChanges
This folder/subfiles calculates the pixel difference between the initial and final frame of the video as a rough measure of movement and behavior.
It contains pixelChange.m (the 'main' function within the folder), and motionStabilization.m
MAIN OUTPUTS: 
- pxlChange (number of worms x number of vid matrix with pixel difference values)
- seasonalSnapshots (updated w/ shift amount across frames, image of subtracted frames)
- deathVid (number of worms x 1 array with video of death for each individual)

> getSegmentedWorms
This folder/subfiles segments the worm using a consensus approach across individuals and across sequential video recordings. It iteratively thresholds and segments out individuals based on morphology (ex. size) and positional information
It contains hotelManager.m (the 'main' function within the folder), circleMe.m, evalExpectations.m, findWorm.m, findWormWrapper.m, getExpectedIntensity.m, getExpectedSizes.m, getGradientField.m, getSeasonStats.m, getThreshold.m, getWidthVals.m,
getWormStrel.m, and segWorm.m.
MAIN OUTPUTS:
- seasonalSnapshots (updated with intensity and threshold value used for segmentation, segmented frame, flagged error frames, worm size, and overall 'beliefs' (i.e. consensus size/intensity values))

> getBehavior
This folder/subfiles extracts behavioral information from the segmented frames, giving information about the centroid speed, amplitude, swimming frequency, etc.
It contains getBehavior.m and getMovementReadouts.m (the 'main' functions within the folder), along with getAmplitude.m, calculateTangentsFromBackbone.m, getBehaviorParameters.m, getBodyLineVelocity.m, getLine.m, plotBodyMovement.m, SimpleEndpoints.m
MAIN OUTPUTS:
- combMovData (number of vids x 1 cell array with movement data) 
- freqData (60 x number of vids + 1 with frequency values over time and movement statistics for behavioral decline)
- ampData (60 x number of vids + 1 with amplitude values over time and movement statistics for behavioral decline)
- csData (60 x number of vids + 1 with centroid speed values over time and movement statistics for behavioral decline)
- PCData (60 x number of vids + 1 with pixel change values over time and movement statistics for behavioral decline)

> behavioralAnalysis
This folder/subfiles contain code for subsequent behavioral analysis with metrics obtained from getBehavior/generalManager. When using these files, we compiled experimental trials together 
It contains compareBehaviors.m (the 'main' function within the folder), plottingBehavioralParameters, combinedRelativeBehavioralDecline, and regressionModel.m
MAIN OUTPUTS:
- combinedBP (number of individuals x 6 cell array with compiled longevity and behavioral metrics across all trails of a selected experimental condition)
- meanVals, relMeanVals, etc. (mean population values of raw movement over time)
- PCA
- LASSO regression model parameters and residuals for inputted data 
