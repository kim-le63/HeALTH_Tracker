function [declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, ...
    declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75,...
    declinePT_duration, relativeHighPeriod_duration] = ...
    getBehaviorParameters(data, fieldOfInterest, HLmethod, getPlot, yVal, deathVids)
%% getBehavioralParameters.m
% % Gets metrics for characterizing behavioral decline over time

% % Inputs:
% % data - cell array of general behavioral parameter to analyze
% % fieldOfInterest - field name of specific parameter to analyze
% % HLmethod - string input to determine which method to determine high v.
% % low activity threshold (either 'kmeans', 'movAvg', or 'none')
% % getPlot - binary determining whether to generate plots for each parameter
% % yVal - y-axis label for plotting
% % deathVids = array for video of death for each worm

% % Outputs:
% % declinePT - nx1 array of videos where activity transitions from high to low for each individual
% % relativeHighPeriod -  nx1 array of duration of high activity/lifespan for each individual
% % declineSlope_all - nx1 array of slope of linear regression of all non NaN datapoints for the individual
% % declineInt_all - nx1 array of intercept of linear regression of all non NaN datapoints for the individual
% % declineSlope_high - nx1 array of slope of linear regression of high (or low for amp) non NaN datapoints for the individual
% % declineInt_high - nx1 array of intercept of linear regression of high (or low for amp) non NaN datapoints for the individual
% % maxValDay - nx1 array of maximum value for each day for each individual
% % avgValDay - nx1 array of average value for each day for each individual
% % maxValDay_high - nx1 array of maximum value for high activity each day for each individual
% % avgValDay_high - nx1 array of average value for high activity each day for each individual
% % threshVal: threshold value between high and low activity
% % declinePT_75: nx1 array of 75% percentile of video number with high activity
% % relativeHighPeriod_75: nx1 array of relative period of high activity (using
% % declinePT_75 as the period of high activity)
% % declinePT_duration: nx1 array of number of videos with high activity
% % relativeHighPeriod_duration: nx1 array of relative period of high activity 
% % (using declinePT_duration as the period of high activity)

%% Finding threshold between high v. low-intensity movement
% setting up cell to workable matrix array
if ~isempty(fieldOfInterest)
    if isfield(data{1,1},fieldOfInterest) == 0
        disp('field does not exist')
        return
    end
 
    behavioralParameter = zeros(size(data,1),size(data,2));
    for i = 1:size(data,2)-1
        for k = 1:size(data,1)
            behavioralParameter(k,i) = data{k,i}.(fieldOfInterest);
        end
    end
else
    behavioralParameter = data;
end

if nargin < 3 
    getPlot =0; HLmethod = 'kmeans';
    yVal = ''; 
elseif nargin < 4
    getPlot = 0; yVal = ''; 
elseif nargin <5 && getPlot == 1
    yVal = '';
end

declinePT = zeros(size(behavioralParameter,1),1);
relativeHighPeriod = zeros(size(behavioralParameter,1),1);
declinePT_75 = zeros(size(behavioralParameter,1),1);
relativeHighPeriod_75 = zeros(size(behavioralParameter,1),1);
declinePT_duration = zeros(size(behavioralParameter,1),1);
relativeHighPeriod_duration = zeros(size(behavioralParameter,1),1);

declineSlope_all = zeros(size(behavioralParameter,1),1);
declineInt_all = zeros(size(behavioralParameter,1),1);
declineSlope_high = zeros(size(behavioralParameter,1),1);
declineInt_high = zeros(size(behavioralParameter,1),1);

if strcmp(HLmethod,'kmeans') == 1
    reshapedBP = reshape(behavioralParameter,[size(behavioralParameter,1)*size(behavioralParameter,2),1]);
    omitNaNBP= reshapedBP(~isnan(reshapedBP));
    kBP = omitNaNBP(omitNaNBP>0);

    idx = kmeans(kBP,2); 
    oneIdx = find(idx==1); twoIdx = find(idx==2);
    one = kBP(oneIdx); two = kBP(twoIdx);
    threshVal = min(max(min(one),min(two)), min(max(one),max(two)))*.9;

    if getPlot == 1
        figure; hold on; 
        scatter([find(idx==1)],one,'*'); scatter([find(idx==2)],two,'o'); 
        xlim([0 180]); ylim([0 1])
        title([data,' all worms- threshVal: ',num2str(threshVal)])
        xlabel('vid num'); ylabel(yVal)
    end
elseif strcmp(HLmethod,'movAvg') == 1
    threshVal = zeros(size(behavioralParameter,1),1);
end

if getPlot ==1
    figure; hold on; 
end

for i = 1:size(behavioralParameter,1)
    tempVal = behavioralParameter(i,:);
    cMovAvg = tempVal(1:find(~isnan(tempVal), 1,'last'));

    if strcmp(HLmethod,'kmeans') == 1 %getting decline point
        tempkAvg = cMovAvg(~isnan(cMovAvg)); %assuming population level threshold
        tempIdx = find(~isnan(cMovAvg));
        highIdx = []; highVal = []; lowIdx = []; lowVal = [];
        for j = 1:length(tempkAvg)
            if tempkAvg(j)>threshVal
                highIdx = cat(1,highIdx,tempIdx(j));
                highVal = cat(1,highVal,tempkAvg(j));
            else
                lowIdx = cat(1,lowIdx,tempIdx(j));
                lowVal = cat(1,lowVal,tempkAvg(j));
            end
        end

        if isempty(highIdx)
            declinePT(i) = NaN; declinePT_75(i) = NaN; declinePT_duration(i) = NaN;
        else
            declinePT(i) = highIdx(end);
            declinePT_75(i) = quantile(highIdx,0.75);
            declinePT_duration(i) = numel(highIdx);
        end

        if getPlot ==1 
            subplot(5,ceil(size(behavioralParameter,1)/5),i); hold on;
            scatter(highIdx,highVal,'*'); scatter(lowIdx,lowVal,'o'); 
            xlim([0 ceil(size(behavioralParameter,2)/5)*5]); ylim([0 ceil(max(max(behavioralParameter))/5)*5]);
            title(['worm ',num2str(i),' declinePT: ',num2str(declinePT(i))])
            xlabel('vid num'); ylabel(yVal)
        end 

    elseif strcmp(HLmethod,'movAvg') == 1
        tempkAvg = cMovAvg(~isnan(cMovAvg)); %finding non-NaN values (for later plotting)
        tempIdx = find(~isnan(cMovAvg));
        
        movAvg = movmedian(tempkAvg,40,'omitnan');
        meanMovAvg = nanmean(movAvg);
        normMovAvg = movAvg-meanMovAvg;
        threshVal(i) = meanMovAvg;
        
        crossLocations = find(diff(sign(normMovAvg)))+1; % find intercept between moving avg and avg of moving avg        
        
        declinePts = crossLocations(find(normMovAvg(crossLocations-1)>0));
        if isempty(declinePts)
            declinePT(i) = NaN; declinePT_75(i) = NaN; declinePT_duration(i) = NaN;
        elseif declinePts(end) == crossLocations(end) && crossLocations(end)~=length(normMovAvg)
            declinePT(i) = declinePts(end); declinePT_75(i) = quantile(declinePts,0.75);
            declinePT_duration(i) = numel(declinePts);
        else
            declinePT(i) = NaN;
        end

        if getPlot == 1
            subplot(5,ceil(size(behavioralParameter,1)/5),i); hold on;
            scatter([1:1:length(cMovAvg)],cMovAvg,'*');
            plot(meanMovAvg*ones(length(movAvg),1),'LineWidth',1,'Color','m')
            plot(movAvg,'LineWidth',2,'Color','k')
            title(['worm ',num2str(i),' decline pt: ',num2str(declinePT(i))])
            xlabel('vid num'); ylabel(yVal)
            xlim([0 ceil(size(behavioralParameter,2)/5)*5]); ylim([0 ceil(max(max(behavioralParameter))/5)*5]);
        end
        
    elseif strcmp(HLmethod,'none') == 1
        tempkAvg = cMovAvg(~isnan(cMovAvg)); %finding non-NaN values (for later plotting)
        tempIdx_NaN = find(~isnan(cMovAvg));
        tempIdx_zero = find(cMovAvg);
        tempIdx = intersect(tempIdx_NaN,tempIdx_zero);
        if length(tempIdx) < 2
            declinePT(i) = NaN; declinePT_75(i) = NaN; declinePT_duration(i) = NaN;
        else
            declinePT(i) = tempIdx(end); declinePT_75(i) = quantile(tempIdx,0.75);
            declinePT_duration(i) = numel(tempIdx);
        end
        
        threshVal = 0;
    else
        disp('input valid thresholding method');
        return
    end

    if isnan(declinePT(i)) %getting relative period of high activity
        relativeHighPeriod(i) = NaN; relativeHighPeriod_75(i) = NaN; 
        relativeHighPeriod_duration(i) = NaN;
        if numel(tempIdx)>1
            t1 = tempIdx; p = polyfit(tempIdx,tempkAvg,1);
            declineSlope_all(i) = p(1); declineInt_all(i) = p(2);
            declineSlope_high(i) = p(1); declineInt_high(i) = p(2);
        else
            declineSlope_all(i) = NaN; declineInt_all(i) = NaN;
            declineSlope_high(i) = NaN; declineInt_high(i) = NaN;
        end
    else
        relativeHighPeriod(i) = declinePT(i)/deathVids(i); %length(cMovAvg);
        relativeHighPeriod_75(i) = declinePT_75(i)/deathVids(i);% length(cMovAvg);
        relativeHighPeriod_duration(i) = declinePT_duration(i)/deathVids(i); %length(cMovAvg);

        if strcmp(HLmethod,'kmeans') == 1
            p = polyfit(tempIdx,tempkAvg,1);
            pHigh = polyfit(highIdx',cMovAvg(highIdx),1);
        elseif strcmp(HLmethod,'movAvg') == 1 %for movAvg method
            p = polyfit(tempIdx,movAvg,1); 
            pHigh = polyfit(tempIdx(declinePT(i):end),movAvg(declinePT(i):end),1);
        else % for no separation between high v. low activity
            cMovAvg_zero = tempVal(1:find(tempVal~=0, 1,'last'));
%             relativeHighPeriod(i) = declinePT(i)/deathVids(i);
            p = polyfit(tempIdx,cMovAvg_zero(tempIdx),1); 
            pHigh = p;
        end
        declineSlope_all(i) = p(1); declineInt_all(i) = p(2);
        declineSlope_high(i) = pHigh(1); declineInt_high(i) = pHigh(2);
    end    
end

% finding daily avg/maxs for the behavioral parameter
maxDay  = round((size(behavioralParameter,2)-2)/4)+1;
maxValDay = zeros(size(behavioralParameter,1),maxDay);
avgValDay = zeros(size(behavioralParameter,1),maxDay);
maxValDay_high = zeros(size(behavioralParameter,1),maxDay);
avgValDay_high = zeros(size(behavioralParameter,1),maxDay);

for i = 1:size(behavioralParameter,1)
    iter = 3;
    if numel(threshVal) > 1
        threshVal2 = threshVal(i);
    else
        threshVal2 = threshVal;
    end
    
    for j = 1:maxDay
        if j==1
            maxValDay(i,j) = max(behavioralParameter(i,1:2));
            avgValDay(i,j) = nanmean(behavioralParameter(i,1:2));
            highBPVals = behavioralParameter(i,(behavioralParameter(i,1:2)>threshVal2));
        elseif j<maxDay
            maxValDay(i,j) = max(behavioralParameter(i,iter:iter+3));
            avgValDay(i,j) = nanmean(behavioralParameter(i,iter:iter+3));
            highBPVals = behavioralParameter(i,iter-1+find(behavioralParameter(i,iter:iter+3)>threshVal2));%+iter-1);
            iter = iter+4;
        else
            maxValDay(i,j) = max(behavioralParameter(i,iter:end));
            avgValDay(i,j) = nanmean(behavioralParameter(i,iter:end));
            highBPVals = behavioralParameter(i,iter-1+find(behavioralParameter(i,iter:end)>threshVal2));%+iter-1);
        end
        
        if ~isempty(highBPVals) 
            maxValDay_high(i,j) = max(highBPVals);
            avgValDay_high(i,j) = nanmean(highBPVals);
        else
            maxValDay_high(i,j) = NaN; avgValDay_high(i,j) = NaN;
        end
    end
end

end
