function [meanVals, adultPC_rel_indiv, relMeanVals, meanVals_short, relMeanVals_short, meanVals_long, ...
    relMeanVals_long, p_ks, AUC] = combinedRelativeBehavioralDecline(y_Vals,X_Vals)
% finding average movement across a population/subpopulations
% INPUTS:
% y_Vals - time of death for each individual
% X_Vals - raw movement over time for each individual
% OUTPUTS:
% meanVals - (number of time pts x 3) array. 1st column is population mean.
% 2nd column is SEM. 3rd is population size
% adultPC_rel_indiv - (number of invididuals x 11 cell array). Has each
% individual's raw movement value for each relative portion of its lifespan
% relMeanVals - (11 x 3 array). mean vals array for relative lifespan
% meanVals_short - mean vals for shortest-lived cohort of population
% meanVals_long - mean vals for longest-lived cohort of population
% relMeanVals_short - mean vals for shortest-lived cohort normalized to
% each invidual's lifespan
% relMeanVals_long - mean vals for longest-lived cohort normalized to each
% invidiual's lifespan
% p_ks - p-val from Kolmogorov-Smirnov test for shortest and longest-lived
% cohort's relative activity normalized to its lifespan
% AUC - difference between area under the curve of shortest and
% longest-lived cohort's relative activity normalized to its lifespan

%% averaging raw movement across and within populations

popY = y_Vals-6;
popX = X_Vals;

aPopX = popX(:,7:end); 

%% hourly movement average across population
hrNum = size(aPopX,2)/2;
meanVals = zeros(hrNum,3); 
iter = 1;
for i =1:hrNum
    tempVals = reshape(aPopX(:,iter:iter+1),size(aPopX,1)*2,1);
    meanVals(i,1) = nanmean(tempVals); %mean
    meanVals(i,2) = nanstd(tempVals)/sqrt(length(tempVals)); %SEM
    meanVals(i,3) = length(tempVals); %number
    iter = iter + 2;
end

%% relative movement average across population
adultPC_rel_raw = cell(11,1);
adultPC_rel_raw_indiv = cell(length(popY),11);
for j = 1:11
    for i = 1:length(popY)
        tSep = linspace(0,11.001,(popY(i)));
        adultPC_rel_raw{j} = [adultPC_rel_raw{j} aPopX(i,tSep<j & tSep>=j-1)];
        adultPC_rel_raw_indiv{i,j} = [adultPC_rel_raw_indiv{i,j} aPopX(i,tSep<j & tSep>=j-1)];
    end
end

adultPC_rel_indiv = zeros(length(popY),11);
for j = 1:11
    for i = 1:length(popY)
        adultPC_rel_indiv(i,j) = nanmean(adultPC_rel_raw_indiv{i,j});
    end
end

relMeanVals = zeros(11,3);
for i =1:11
    tempVals = adultPC_rel_raw{i};
    tempVals = tempVals; % tempVals(find(tempVals));
    relMeanVals(i,1) = nanmean(tempVals); %mean
    relMeanVals(i,2) = nanstd(tempVals)/sqrt(length(tempVals)); %SEM
    relMeanVals(i,3) = length(tempVals); %number
end

 %% using percentiles
shortVal = prctile(popY,20); % 20th percentile
longVal = prctile(popY,80); % 80th percentile
 
shortIdx = [find(popY<=shortVal)]; longIdx = [find(popY>=longVal)];

% short-lived cohort average
meanVals_short = zeros(hrNum,3);
iter = 1;
for i =1:hrNum
    tempVals_short = reshape(aPopX(shortIdx,iter:iter+1),length(shortIdx)*2,1);
    meanVals_short(i,1) = nanmean(tempVals_short); %mean
    meanVals_short(i,2) = nanstd(tempVals_short)/sqrt(length(tempVals_short)); %SEM
    meanVals_short(i,3) = length(tempVals_short); %number
    iter = iter+2;
end

% short-lived relative movement average 
adultPC_rel_raw_short = cell(11,1);
for j = 1:11
    for i = 1:length(shortIdx)
        tSep = linspace(0,11.001,(popY(shortIdx(i))));
        adultPC_rel_raw_short{j} = [adultPC_rel_raw_short{j} aPopX(shortIdx(i),tSep<j & tSep>=j-1)];
    end
end

relMeanVals_short = zeros(11,3);
for i =1:11
    relMeanVals_short(i,1) = nanmean(adultPC_rel_raw_short{i}); %mean
    relMeanVals_short(i,2) = nanstd(adultPC_rel_raw_short{i})/sqrt(length(adultPC_rel_raw_short{i})); %SEM
    relMeanVals_short(i,3) = length(adultPC_rel_raw_short{i}); %number
end

% long-lived cohort average
meanVals_long = zeros(hrNum,3); %CHECK
iter = 1;
for i =1:hrNum
    tempVals_long = reshape(aPopX(longIdx,iter:iter+1),length(longIdx)*2,1);
    meanVals_long(i,1) = nanmean(tempVals_long); %mean
    meanVals_long(i,2) = nanstd(tempVals_long)/sqrt(length(tempVals_long)); %SEM
    meanVals_long(i,3) = length(tempVals_long); %number
    iter = iter + 2;
end

% long-lived relative movement average 
adultPC_rel_raw_long = cell(11,1);
for j = 1:11
    for i = 1:length(longIdx)
        tSep = linspace(0,11.001,(popY(longIdx(i))));
        adultPC_rel_raw_long{j} = [adultPC_rel_raw_long{j} aPopX(longIdx(i),tSep<j & tSep>=j-1)];
    end
end

relMeanVals_long = zeros(11,3);
for i =1:11
    relMeanVals_long(i,1) = nanmean(adultPC_rel_raw_long{i}); %mean
    relMeanVals_long(i,2) = nanstd(adultPC_rel_raw_long{i})/sqrt(length(adultPC_rel_raw_long{i})); %SEM
    relMeanVals_long(i,3) = length(adultPC_rel_raw_long{i}); %number
end

% KS test
[~,p_ks] = kstest2(relMeanVals_short(:,1) ,relMeanVals_long(:,1));

%difference between AUC
AUC = trapz(relMeanVals_short(:,1)) - trapz(relMeanVals_long(:,1));
end
