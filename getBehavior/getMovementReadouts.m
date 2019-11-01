function [freqData,ampData,csData,dpData,areaData] = getMovementReadouts(combMovData, deathVids_init, chamberIDs)
%% getMovementReadouts.m
% Takes raw behavioral readouts (from combMoveData and getBehavior.m) and
% outputs it into a parseable format with behavioral decline metrics for
% each behavioral metric 

% % Inputs:
% combMovData - cell array of movement data across all videos for a single device over time (ex. segment error, tangents, centroid, amplitude, etc.)
% deathVids_init - time of death for each worm
% chamberIDs - list of chambers to analyze

% % Outputs:
% freqData - cell array (60 x number of videos + 1) of frequency values over time. 
%            columns 1 - num of videos have mean value, max value, and # of segments undergoing sinusoidal motion for each video. 
%            Last column has summary of movement statistics for the
%            different mean, max, and segment number stats describing
%            behavioral decline and duration of high v. low activity
%            movement. (see getBehavioralParamters for more info)
% ampData - cell array (60 x number of videos + 1) of amplitude values over time. 
%           Similar format to freqData (except looking at just mean and max
%           values)
% csData - cell array (60 x number of videos + 1) of centroid speed values over time. 
%           Similar format to freqData (except looking at just mean and max
%           values)
% dpData - cell array (60 x number of videos + 1) of delta pixel values over time. 
%           Similar format to freqData (except looking at just mean and max
%           values)
% areaData - cell array (60 x number of videos + 1) of area values over time. 
%           Similar format to freqData (except looking at just mean and max
%           values)
%%

movementTime = max(find(~cellfun('isempty', combMovData))); 

freqData = cell(length(chamberIDs),movementTime+1); 
ampData = cell(length(chamberIDs),movementTime+1); 
csData = cell(length(chamberIDs),movementTime+1); 
dpData = cell(length(chamberIDs),movementTime+1); 
areaData = cell(length(chamberIDs),movementTime+1); 

minWidthVal = 3; maxWidthVal = 30;
for i = 1:movementTime %for each video
    movementData = combMovData{i,1}.movData;

    for k = 1:length(chamberIDs) %for each individual
        tempTang = zeros(25,size(movementData,2));
        tempAmp = zeros(size(movementData,2),1);
        tempCS = zeros(size(movementData,2),1);
        tempDP = zeros(size(movementData,2),1);
        tempArea = zeros(size(movementData,2),1);
         
        for j = 1:size(movementData,2) %for each frame
            if isfield(movementData{chamberIDs(k),j},'Amplitude') == 1
                tempAmp(j) = NaN; tempCS(j) = NaN; tempDP(j) = NaN;
                tempArea(j) = NaN; tempTang(:,j) = zeros(25,1); 
            elseif isfield(movementData{chamberIDs(k),j},'amplitude') == 0
                tempAmp(j) = NaN; tempCS(j) = NaN; tempDP(j) = NaN;
                tempArea(j) = NaN; tempTang(:,j) = zeros(25,1); 
            else
                tempAmp(j) = movementData{chamberIDs(k),j}.amplitude;
                tempCS(j) = movementData{chamberIDs(k),j}.centroidSpeed;
                tempDP(j) = movementData{chamberIDs(k),j}.deltaPix;
                tempArea(j) = movementData{chamberIDs(k),j}.area;
                tempTang(:,j) = movementData{chamberIDs(k),j}.tangents; 
            end
        end
        
        % checking for bad frames (via area and centroid speed sanity checks/IDing outliers)     
        [TF,L,~,C] = isoutlier(tempArea,'movmedian',14);
        outlierVid = zeros(length(tempArea),1);
        thresh = .05;
        for j = 1:size(movementData,2)
            if ((tempArea(j) <L(j)*(1-thresh) && ~isnan(tempArea(j)))) || tempCS(j) > 7
                tempAmp(j) = NaN; tempCS(j) = NaN; tempDP(j) = NaN; 
                tempArea(j) = NaN;
            end
        end
        
        ampData{k,i}.maxVal = max(tempAmp); ampData{k,i}.meanVal = nanmean(tempAmp);
        csData{k,i}.maxVal = max(tempCS); csData{k,i}.meanVal = nanmean(tempCS);
        dpData{k,i}.maxVal = max(tempDP); dpData{k,i}.meanVal = nanmean(tempDP);
        areaData{k,i}.meanVal = nanmean(tempArea);
        
        % getting frequency values from tang
        iter = 4;
        bodySeg =[2,7,12,18,24];
        freqData{k,i}.freqVals = zeros(length(bodySeg),1);
        freqData{k,i}.numPeaks = zeros(length(bodySeg),1);
        freqData{k,i}.maxPeakHeight = zeros(length(bodySeg),1);

        for j = 1:length(bodySeg)
            tempSegTang = tempTang(bodySeg(j),:);
            % band pass filtering signal (0.25 - 3 Hz)
            fs = 14; N  = length(tempSegTang);
            dF = fs/N; f  = (-fs/2:dF:fs/2-dF)';
            BPF = ((0.3 < abs(f)) & (abs(f) < 7));

            tempSegTang = tempSegTang-nanmean(tempSegTang);
            spektrum = fftshift(fft(tempSegTang))/N;
            spektrum = BPF'.*spektrum;
            fTempSegTang =  ifft(ifftshift(spektrum)); 
            fTempSegTang = real(fTempSegTang);
            
            crossLoc = find(diff(sign(fTempSegTang)));
            diffCL = diff(crossLoc);
            peakWidthInRange = diffCL >= minWidthVal & diffCL <= maxWidthVal;
            crossNum = 0; 
            for m=1:length(diffCL)-1
                if peakWidthInRange(m) == 1 && peakWidthInRange(m+1) ==1
                    if diffCL(m+1)/2 < 3
                        buffer= 3;
                    else
                        buffer = diffCL(m+1)/2;
                    end
                    if (diffCL(m+1) <= diffCL(m) + buffer) && (diffCL(m+1) >= diffCL(m) -buffer)
                        crossNum = crossNum + 1;
                        if crossNum ==3
                            break
                        end
                    else
                        crossNum = 0;
                    end
                else
                    crossNum = 0;
                end
            end
                
            if crossNum == 3
                [pxx,f] = periodogram(fTempSegTang,[],[],fs);
                % find peaks in PSD
                [maxVal,maxIdx] = max(pxx);
                [pkx,pIdx] = findpeaks(pxx,'MinPeakHeight',maxVal/2);
                if numel(pkx) > 2 || maxVal < 0.000005 
                    freqData{k,i}.freqVals(j) = NaN;
                else
                    freqData{k,i}.freqVals(j) = f(maxIdx);
                    freqData{k,i}.numPeaks(j) = numel(pkx);
                    freqData{k,i}.maxPeakHeight(j) = maxVal; 
                end
            else
                freqData{k,i}.freqVals(j) = NaN;
                freqData{k,i}.numPeaks(j) = 0;
                freqData{k,i}.maxPeakHeight(j) = 0; 
            end
            iter = iter+3;
        end
        
        freqData{k,i}.segVals = find(~isnan(freqData{k,i}.freqVals));
        freqData{k,i}.numSeg = numel(find(~isnan(freqData{k,i}.freqVals)));
        if isempty(freqData{k,i}.segVals)
            freqData{k,i}.segVals=NaN;
            freqData{k,i}.numSeg=0;            
        end
        
        if sum(isnan(freqData{k,i}.freqVals)) == 4 %i.e. only 1 segment is moving
            if freqData{k,i}.numPeaks(find(~isnan(freqData{k,i}.freqVals))) > 1 %if there are multiple peaks, censor
                freqData{k,i}.numPeaks(find(~isnan(freqData{k,i}.freqVals))) = 0;
                freqData{k,i}.maxPeakHeight(find(~isnan(freqData{k,i}.freqVals))) = 0;
                freqData{k,i}.freqVals(~isnan(freqData{k,i}.freqVals)) = NaN;
            end
        end
        % going through segment cases
        if sum(isnan(freqData{k,i}.freqVals)) <5 && sum(isnan(freqData{k,i}.freqVals)) >1 
            diffFD = diff(freqData{k,i}.segVals);
            if isempty(find(diffFD<3)) %no consecutive body movements
                freqData{k,i}.numPeaks(find(~isnan(freqData{k,i}.freqVals))) = 0;
                freqData{k,i}.maxPeakHeight(find(~isnan(freqData{k,i}.freqVals))) = 0;
                freqData{k,i}.freqVals(find(~isnan(freqData{k,i}.freqVals))) = NaN;
            end
        end
        freqData{k,i}.meanVal = nanmean(freqData{k,i}.freqVals);
        freqData{k,i}.maxVal = max(freqData{k,i}.freqVals);
    end
end

%% checking for area outliers (indicator of bad segmentation) across videos to censor videos from analysis
tempPopArea = NaN(length(chamberIDs),movementTime); %converting from cell to workable matrix
for i = 1:movementTime
    for k =1:length(chamberIDs)
        tempPopArea(k,i) = areaData{k,i}.meanVal;
    end
end

for k = 1:length(chamberIDs)
    tempArea = tempPopArea(k,:);
    t = 1:1:size(tempArea,2);
    
    [TF,L,U,C] = isoutlier(tempArea,'movmedian',14);
    outlierVid = zeros(length(tempArea),1);
    thresh = .05;
    for i=1:length(tempArea)
        if (tempArea(i) <L(i)*(1-thresh) && ~isnan(tempArea(i))) || (tempArea(i) >=475 && i <=30)
            ampData{k,i}.maxVal = NaN; ampData{k,i}.meanVal = NaN;
            csData{k,i}.maxVal = NaN; csData{k,i}.meanVal = NaN;
            dpData{k,i}.maxVal = NaN; dpData{k,i}.meanVal = NaN;
            areaData{k,i}.meanVal = NaN;
        end
    end
end

% last check for outliers with populational mean
mean2Area = nanmean(tempPopArea);
for i = 1:length(mean2Area)
    [TF,L,U,C] = isoutlier(mean2Area,'movmedian',10);
    outlierVid = zeros(length(mean2Area),1);
    thresh = .05;
    
    if (mean2Area(i) <L(i)*(1-thresh))
            outlierVid(i) = 1;
        for k = 1:length(chamberIDs)
            ampData{k,i}.maxVal = NaN; ampData{k,i}.meanVal = NaN;
            csData{k,i}.maxVal = NaN; csData{k,i}.meanVal = NaN;
            dpData{k,i}.maxVal = NaN; dpData{k,i}.meanVal = NaN;
            areaData{k,i}.meanVal = NaN;
        end
    end
end

% % %% Plot everything here
% % plotBodyMovement(areaData,'meanVal','vid number','area(pxl)','Area ',chamberIDs,0);
% % plotBodyMovement(ampData,'maxVal','vid number','amplitude (pxl)','Max Amp ',chamberIDs,0);
% % plotBodyMovement(ampData,'meanVal','vid number','amplitude (pxl)','Mean Amp ',chamberIDs,0);
% % plotBodyMovement(maxDeltaPixData,'maxVal', 'vid number','max DP(pxl)','Max Delta Pix ',chamberIDs,0);
% % plotBodyMovement(meanDeltaPixData,'meanVal','vid number','mean DP(pxl)','Mean Delta Pix ',chamberIDs,0);
% % plotBodyMovement(maxCentroidSpeedData,'maxVal', 'vid number','max speed','Max Centroid Speed ',chamberIDs,0);
% % plotBodyMovement(meanCentroidSpeedData,'meanVal','vid number','mean speed','Mean Centroid Speed ',chamberIDs,0);
% % plotBodyMovement(freqDataDay, 'time(days)','freq(Hz)','Frequency ',chamberIDs);
% % plotBodyMovement(freqData ,'avgFreqVal','vid number','freq(Hz)','Freq ',chamberIDs,0);
% % plotBodyMovement(freqData,'meanValDay','vid number','freq(Hz)','Freq ',chamberIDs,0);
% % plotBodyMovement(freqData,'maxNumSeg','vid number','freq(Hz)','Freq ',chamberIDs,0);
% % plotBodyMovement(freqData,'meanNumSeg','vid number','freq(Hz)','Freq ',chamberIDs,1);

%% amplitude
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(ampData, 'maxVal', 'movAvg', 0, 'max amp ', deathVids_init);
for k =1:length(chamberIDs)
    ampData{k,movementTime+1}.maxVal.declinePT = declinePT(k); ampData{k,movementTime+1}.maxVal.relativeHighPeriod = relativeHighPeriod(k);
    ampData{k,movementTime+1}.maxVal.declineSlope_all = declineSlope_all(k); ampData{k,movementTime+1}.maxVal.declineInt_all = declineInt_all(k);
    ampData{k,movementTime+1}.maxVal.declineSlope_high = declineSlope_high(k); ampData{k,movementTime+1}.maxVal.declineInt_high = declineInt_high(k);
    ampData{k,movementTime+1}.maxVal.maxValDay = maxValDay(k,:); ampData{k,movementTime+1}.maxVal.avgValDay = avgValDay(k,:);
    ampData{k,movementTime+1}.maxVal.maxValDay_high = maxValDay_high(k,:); ampData{k,movementTime+1}.maxVal.avgValDay_high = avgValDay_high(k,:);
    ampData{k,movementTime+1}.maxVal.threshVal = threshVal(k); ampData{k,movementTime+1}.maxVal.declinePT_75 = declinePT_75(k); 
    ampData{k,movementTime+1}.maxVal.relativeHighPeriod_75 = relativeHighPeriod_75(k); 
    ampData{k,movementTime+1}.maxVal.declinePT_duration = declinePT_duration(k); 
    ampData{k,movementTime+1}.maxVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(ampData, 'meanVal', 'movAvg', 0, 'max amp ',deathVids_init);
for k =1:length(chamberIDs)
    ampData{k,movementTime+1}.meanVal.declinePT = declinePT(k); ampData{k,movementTime+1}.meanVal.relativeHighPeriod = relativeHighPeriod(k);
    ampData{k,movementTime+1}.meanVal.declineSlope_all = declineSlope_all(k); ampData{k,movementTime+1}.meanVal.declineInt_all = declineInt_all(k);
    ampData{k,movementTime+1}.meanVal.declineSlope_high = declineSlope_high(k); ampData{k,movementTime+1}.meanVal.declineInt_high = declineInt_high(k);
    ampData{k,movementTime+1}.meanVal.maxValDay = maxValDay(k,:); ampData{k,movementTime+1}.meanVal.avgValDay = avgValDay(k,:);
    ampData{k,movementTime+1}.meanVal.maxValDay_high = maxValDay_high(k,:); ampData{k,movementTime+1}.meanVal.avgValDay_high = avgValDay_high(k,:);
    ampData{k,movementTime+1}.meanVal.threshVal = threshVal(k);
    ampData{k,movementTime+1}.meanVal.declinePT_75 = declinePT_75(k); ampData{k,movementTime+1}.meanVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    ampData{k,movementTime+1}.meanVal.declinePT_duration = declinePT_duration(k); ampData{k,movementTime+1}.meanVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
%% delta pixel change
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(dpData, 'maxVal', 'kmeans', 0, 'max amp ', deathVids_init);
for k =1:length(chamberIDs)
    dpData{k,movementTime+1}.maxVal.declinePT = declinePT(k); dpData{k,movementTime+1}.maxVal.relativeHighPeriod = relativeHighPeriod(k);
    dpData{k,movementTime+1}.maxVal.declineSlope_all = declineSlope_all(k); dpData{k,movementTime+1}.maxVal.declineInt_all = declineInt_all(k);
    dpData{k,movementTime+1}.maxVal.declineSlope_high = declineSlope_high(k); dpData{k,movementTime+1}.maxVal.declineInt_high = declineInt_high(k);
    dpData{k,movementTime+1}.maxVal.maxValDay = maxValDay(k,:); dpData{k,movementTime+1}.maxVal.avgValDay = avgValDay(k,:);
    dpData{k,movementTime+1}.maxVal.maxValDay_high = maxValDay_high(k,:); dpData{k,movementTime+1}.maxVal.avgValDay_high = avgValDay_high(k,:);
    dpData{k,movementTime+1}.maxVal.threshVal = threshVal;
    dpData{k,movementTime+1}.maxVal.declinePT_75 = declinePT_75(k); dpData{k,movementTime+1}.maxVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    dpData{k,movementTime+1}.maxVal.declinePT_duration = declinePT_duration(k); dpData{k,movementTime+1}.maxVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(dpData, 'meanVal', 'kmeans', 0, 'max amp ',deathVids_init);
for k =1:length(chamberIDs)
    dpData{k,movementTime+1}.meanVal.declinePT = declinePT(k); dpData{k,movementTime+1}.meanVal.relativeHighPeriod = relativeHighPeriod(k);
    dpData{k,movementTime+1}.meanVal.declineSlope_all = declineSlope_all(k); dpData{k,movementTime+1}.meanVal.declineInt_all = declineInt_all(k);
    dpData{k,movementTime+1}.meanVal.declineSlope_high = declineSlope_high(k); dpData{k,movementTime+1}.meanVal.declineInt_high = declineInt_high(k);
    dpData{k,movementTime+1}.meanVal.maxValDay = maxValDay(k,:); dpData{k,movementTime+1}.meanVal.avgValDay = avgValDay(k,:);
    dpData{k,movementTime+1}.meanVal.maxValDay_high = maxValDay_high(k,:); dpData{k,movementTime+1}.meanVal.avgValDay_high = avgValDay_high(k,:);
    dpData{k,movementTime+1}.meanVal.threshVal = threshVal;
    dpData{k,movementTime+1}.meanVal.declinePT_75 = declinePT_75(k); dpData{k,movementTime+1}.meanVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    dpData{k,movementTime+1}.meanVal.declinePT_duration = declinePT_duration(k); dpData{k,movementTime+1}.meanVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
%% centroid speed
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(csData, 'maxVal', 'kmeans', 0, 'max amp ',deathVids_init);
for k =1:length(chamberIDs)
    csData{k,movementTime+1}.maxVal.declinePT = declinePT(k); csData{k,movementTime+1}.maxVal.relativeHighPeriod = relativeHighPeriod(k);
    csData{k,movementTime+1}.maxVal.declineSlope_all = declineSlope_all(k); csData{k,movementTime+1}.maxVal.declineInt_all = declineInt_all(k);
    csData{k,movementTime+1}.maxVal.declineSlope_high = declineSlope_high(k); csData{k,movementTime+1}.maxVal.declineInt_high = declineInt_high(k);
    csData{k,movementTime+1}.maxVal.maxValDay = maxValDay(k,:); csData{k,movementTime+1}.maxVal.avgValDay = avgValDay(k,:);
    csData{k,movementTime+1}.maxVal.maxValDay_high = maxValDay_high(k,:); csData{k,movementTime+1}.maxVal.avgValDay_high = avgValDay_high(k,:);
    csData{k,movementTime+1}.maxVal.threshVal = threshVal;
    csData{k,movementTime+1}.maxVal.declinePT_75 = declinePT_75(k); csData{k,movementTime+1}.maxVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    csData{k,movementTime+1}.maxVal.declinePT_duration = declinePT_duration(k); csData{k,movementTime+1}.maxVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all,declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(csData, 'meanVal', 'kmeans', 0, 'max amp ',deathVids_init);
for k =1:length(chamberIDs)
    csData{k,movementTime+1}.meanVal.declinePT = declinePT(k); csData{k,movementTime+1}.meanVal.relativeHighPeriod = relativeHighPeriod(k);
    csData{k,movementTime+1}.meanVal.declineSlope_all = declineSlope_all(k); csData{k,movementTime+1}.meanVal.declineInt_all = declineInt_all(k);
    csData{k,movementTime+1}.meanVal.declineSlope_high = declineSlope_high(k); csData{k,movementTime+1}.meanVal.declineInt_high = declineInt_high(k);
    csData{k,movementTime+1}.meanVal.maxValDay = maxValDay(k,:); csData{k,movementTime+1}.meanVal.avgValDay = avgValDay(k,:);
    csData{k,movementTime+1}.meanVal.maxValDay_high = maxValDay_high(k,:); csData{k,movementTime+1}.meanVal.avgValDay_high = avgValDay_high(k,:);
    csData{k,movementTime+1}.meanVal.threshVal = threshVal;
    csData{k,movementTime+1}.meanVal.declinePT_75 = declinePT_75(k); csData{k,movementTime+1}.meanVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    csData{k,movementTime+1}.meanVal.declinePT_duration = declinePT_duration(k); csData{k,movementTime+1}.meanVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
%% frequency
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(freqData, 'meanVal', 'none', 1, 'freq(hz)',deathVids_init);
for k =1:length(chamberIDs)
    freqData{k,movementTime+1}.meanVal.declinePT = declinePT(k); freqData{k,movementTime+1}.meanVal.relativeHighPeriod = relativeHighPeriod(k);
    freqData{k,movementTime+1}.meanVal.declineSlope_all = declineSlope_all(k); freqData{k,movementTime+1}.meanVal.declineInt_all = declineInt_all(k);
    freqData{k,movementTime+1}.meanVal.declineSlope_high = declineSlope_high(k); freqData{k,movementTime+1}.meanVal.declineInt_high = declineInt_high(k);
    freqData{k,movementTime+1}.meanVal.maxValDay = maxValDay(k,:); freqData{k,movementTime+1}.meanVal.avgValDay = avgValDay(k,:);
    freqData{k,movementTime+1}.meanVal.maxValDay_high = maxValDay_high(k,:); freqData{k,movementTime+1}.meanVal.avgValDay_high = avgValDay_high(k,:);
    freqData{k,movementTime+1}.meanVal.threshVal = threshVal;
    freqData{k,movementTime+1}.meanVal.declinePT_75 = declinePT_75(k); freqData{k,movementTime+1}.meanVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    freqData{k,movementTime+1}.meanVal.declinePT_duration = declinePT_duration(k); freqData{k,movementTime+1}.meanVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(freqData, 'maxVal', 'none', 1, 'freq(hz)',deathVids_init);
for k =1:length(chamberIDs)
    freqData{k,movementTime+1}.maxVal.declinePT = declinePT(k); freqData{k,movementTime+1}.maxVal.relativeHighPeriod = relativeHighPeriod(k);
    freqData{k,movementTime+1}.maxVal.declineSlope_all = declineSlope_all(k); freqData{k,movementTime+1}.maxVal.declineInt_all = declineInt_all(k);
    freqData{k,movementTime+1}.maxVal.declineSlope_high = declineSlope_high(k); freqData{k,movementTime+1}.maxVal.declineInt_high = declineInt_high(k);
    freqData{k,movementTime+1}.maxVal.maxValDay = maxValDay(k,:); freqData{k,movementTime+1}.maxVal.avgValDay = avgValDay(k,:);
    freqData{k,movementTime+1}.maxVal.maxValDay_high = maxValDay_high(k,:); freqData{k,movementTime+1}.maxVal.avgValDay_high = avgValDay_high(k,:);
    freqData{k,movementTime+1}.maxVal.threshVal = threshVal;
    freqData{k,movementTime+1}.maxVal.declinePT_75 = declinePT_75(k); freqData{k,movementTime+1}.maxVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    freqData{k,movementTime+1}.maxVal.declinePT_duration = declinePT_duration(k); freqData{k,movementTime+1}.maxVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end

[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(freqData, 'numSeg', 'none', 1, 'Num Seg',deathVids_init);
for k =1:length(chamberIDs)
    freqData{k,movementTime+1}.numSeg.declinePT = declinePT(k); freqData{k,movementTime+1}.numSeg.relativeHighPeriod = relativeHighPeriod(k);
    freqData{k,movementTime+1}.numSeg.declineSlope_all = declineSlope_all(k); freqData{k,movementTime+1}.numSeg.declineInt_all = declineInt_all(k);
    freqData{k,movementTime+1}.numSeg.declineSlope_high = declineSlope_high(k); freqData{k,movementTime+1}.numSeg.declineInt_high = declineInt_high(k);
    freqData{k,movementTime+1}.numSeg.maxValDay = maxValDay(k,:); freqData{k,movementTime+1}.numSeg.avgValDay = avgValDay(k,:);
    freqData{k,movementTime+1}.numSeg.maxValDay_high = maxValDay_high(k,:); freqData{k,movementTime+1}.numSeg.avgValDay_high = avgValDay_high(k,:);
    freqData{k,movementTime+1}.numSeg.threshVal = threshVal;
    freqData{k,movementTime+1}.numSeg.declinePT_75 = declinePT_75(k); freqData{k,movementTime+1}.numSeg.relativeHighPeriod_75 = relativeHighPeriod_75(k);
    freqData{k,movementTime+1}.numSeg.declinePT_duration = declinePT_duration(k); freqData{k,movementTime+1}.numSeg.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
end

end
