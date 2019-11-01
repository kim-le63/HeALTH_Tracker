%% plottingBehavioralParameters.m
function [combinedBP] = plottingBehavioralParameters(expConditions, folderName)
%% getting all behavioral data across all trials for specified conditions (for raw data)
% INPUT:
% expConditions - 1x1 cell array with experimental condition to analyze
% folderName - folder path where files are saved
% OUTPUT:
% combinedBP - nx6 cell array with lifespans and behavioral data across all
% trials for specified experimental conditions. column 1 - sorted death
% videos. column 2 - sorted lifespans (in days). column 3 - number of
% individuals x behavioral decline metrics for each individual
% column 4 - variable names for the behavioral metrics. column 5 - number of
% individuals within an experimental condition x number of vids x 10 cell 
% array with all of the raw behavioral metrics for each video. Individuals 
% are sorted by time of death. Columns within the cell array contain 
% ampData(1-2), csData(3-4), dpData(5-6), freqData(7-9), and PCData(10). 
% column 6- sorted indices of individuals

fileList_lifespan = dir(fullfile(folderName,'lifespanData_*.*')); fileList_lifespan = {fileList_lifespan.name};
fileList_movementData = dir(fullfile(folderName,'movementData_*.*')); fileList_movementData = {fileList_movementData.name};
fileList_chamberData = dir(fullfile(folderName,'chamberData_*.*')); fileList_chamberData = {fileList_chamberData.name};

combinedBP = cell(6,1);

if length(expConditions) == 1
    fileList_lifespan_2 = fileList_lifespan(contains(fileList_lifespan,expConditions{1}));
    fileList_movementData_2 = fileList_movementData(contains(fileList_movementData,expConditions{1}));
    fileList_chamberData_2 = fileList_chamberData(contains(fileList_chamberData,expConditions{1}));
else
    fileList_lifespan_2 = []; fileList_movementData_2 = []; fileList_chamberData_2 = [];
    for i = 1:length(expConditions)
        fileList_lifespan_2 = cat(2,fileList_lifespan_2,fileList_lifespan(contains(fileList_lifespan,expConditions{i})));
        fileList_movementData_2 = cat(2,fileList_movementData_2,fileList_movementData(contains(fileList_movementData,expConditions{i})));
        fileList_chamberData_2 = cat(2,fileList_chamberData_2,fileList_chamberData(contains(fileList_chamberData,expConditions{i})));
    end
    
     fileList_lifespan_2 = []; fileList_movementData_2 = []; fileList_chamberData_2 = [];
    for i = 1:length(expConditions)
        fileList_lifespan_2 = cat(2,fileList_lifespan_2,fileList_lifespan(contains(fileList_lifespan,expConditions{i})));
        fileList_movementData_2 = cat(2,fileList_movementData_2,fileList_movementData(contains(fileList_movementData,expConditions{i})));
        fileList_chamberData_2 = cat(2,fileList_chamberData_2,fileList_chamberData(contains(fileList_chamberData,expConditions{i})));
    end
    
    if any(contains(expConditions,'~')) % set up exclude function (ex. 15C <-> 20C oscillation condition)
        cIdx = contains(fileList_lifespan_2,expConditions{find(contains(expConditions,'~'))}(2:end));
        fileList_lifespan_2 = fileList_lifespan_2(~cIdx);
        fileList_movementData_2 = fileList_movementData_2(~cIdx);
        fileList_chamberData_2 = fileList_chamberData_2(~cIdx);
    end
end
% getting lifespans for all individuals under the same experimental
% condition
y_val = []; y_val_day = [];
for j = 1:length(fileList_lifespan_2)
    fullPath = [folderName '\' fileList_lifespan_2{j}]; tempLife = load(fullPath);
    y_val = cat(1,y_val,tempLife.([fileList_lifespan_2{j}(1:end-4)]).deathVid);
    y_val_day = cat(1,y_val_day,tempLife.([fileList_lifespan_2{j}(1:end-4)]).lifeSpan);
end

% find max lifespan within experimental condition
maxVid = max(y_val_day);

% getting behavioral parameters for all individuals under the same
% experimental conditions
X_val = []; X_decline = [];
for j = 1:length(fileList_movementData_2)
    fullPath = [folderName '\' fileList_movementData_2{j}]; tempMove = load(fullPath);
    ampData = tempMove.([fileList_movementData_2{j}(1:end-4)]).ampData;
    csData = tempMove.([fileList_movementData_2{j}(1:end-4)]).csData;
    dpData = tempMove.([fileList_movementData_2{j}(1:end-4)]).dpData;
    freqData = tempMove.([fileList_movementData_2{j}(1:end-4)]).freqData;
    PCData = tempMove.([fileList_movementData_2{j}(1:end-4)]).PCData;

%         tempDeclineX = getDeclinePts(ampData,csData,dpData,freqData,PCData);
%         X_decline = cat(1,X_decline,tempDeclineX);
    [tempDeclineX,iVarNames] = getBehavioralDecline(ampData,csData,dpData,freqData,PCData); 
    X_decline = cat(1,X_decline,tempDeclineX);


    tempX = getRawData(ampData,csData,dpData,freqData,PCData);
    if j > 1 && size(tempX,2) > size(X_val,2)
        padX = zeros(size(X_val,1),(size(tempX,2)-size(X_val,2)),10);
        X_val = cat(2, X_val,padX);
    elseif j > 1 && size(tempX,2) < size(X_val,2)
        padX = zeros(size(tempX,1),(size(X_val,2)-size(tempX,2)),10);
        tempX = cat(2,tempX,padX);
    end
    X_val = cat(1,X_val,tempX);
end

% checking for censored individuals
y_censor = [];
for j = 1:length(fileList_chamberData_2)
    fullPath = [folderName '\' fileList_chamberData_2{j}]; tempChamber = load(fullPath);
    if isfield(tempChamber.([fileList_chamberData_2{j}(1:end-4)]),'censor')
        tempCensor = tempChamber.([fileList_chamberData_2{j}(1:end-4)]).censor;
    else
        tempCensor = zeros(length(tempChamber.([fileList_chamberData_2{j}(1:end-4)]).chamberIDs),1);
    end
    y_censor = cat(1,y_censor,tempCensor);
end

censorIdx = find(~y_censor);
y_val_censor = y_val(censorIdx);
y_val_day_censor = y_val_day(censorIdx);
X_val_censor = X_val(censorIdx,:,:);
X_decline_censor = X_decline(censorIdx,:);

[s_val,sIdx_val] = sort(y_val_censor);

mask = zeros(size(X_val_censor,1),size(X_val_censor,2));
for j =1:size(y_val_censor,1)
   mask(j,1:y_val_censor(j)) = 1;
end

maskX = X_val_censor.*mask;

for i = 1:size(maskX,1)
    for j = 1:size(maskX,2)
        if mask(i,j)==0
            maskX(i,j,:) = 0;
        end
    end
end

s_maskX = maskX(sIdx_val,:,:);

combinedBP{1} = s_val; combinedBP{2} = y_val_day_censor(sIdx_val);
combinedBP{3} = X_decline_censor(sIdx_val,:); combinedBP{4} = iVarNames;
combinedBP{5} = s_maskX; combinedBP{6} = sIdx_val;

end


%% getting behavioral decline values from movementData cell array into parsable, matrix form
function [X,iVarNames] = getBehavioralDecline(ampData,csData,dpData,freqData,PCData) 
nWorms = size(ampData,1);
X = zeros(nWorms,1);
iVarNames = {};

behavioralParameters = {'ampData';'csData';'dpData';'freqData';'PCData'};
iter = 1;
for i = 1:numel(behavioralParameters)
    bpVar = eval(behavioralParameters{i});

    tempData = bpVar{1,size(bpVar,2)};
    tempFields = fieldnames(tempData); 
    for m = 1:length(tempFields)
        tempData2 = tempData.(tempFields{m});
        if isstruct(tempData2) == 1
            tempFields2 = fieldnames(tempData2); 
            for mm = 1:length(tempFields2)
                tempData3 = tempData2.(tempFields2{mm});
                %name column variable name
                iVarNames{iter} = [behavioralParameters{i},'_',tempFields{m},'_',tempFields2{mm}];
                if numel(tempData3) == 1
                    for k = 1:nWorms
                        X(k,iter) = bpVar{k,size(bpVar,2)}.(tempFields{m}).(tempFields2{mm});
                        if i==4 && isnan(X(k,iter))
                            X(k,iter) = 0;
                        end
                    end
                    iter = iter + 1;
                else %for maxValDay
                    dayCount = 0;
                    for mmm = 1:11 %Up to Day 10 Adult
                        iVarNames{iter} = [behavioralParameters{i},'_',tempFields{m},'_',tempFields2{mm},'_Day',num2str(dayCount)];
                        for k = 1:nWorms
                            X(k,iter) = bpVar{k,size(bpVar,2)}.(tempFields{m}).(tempFields2{mm})(mmm);
                            if i==4 && isnan(X(k,iter))
                                X(k,iter) = 0;
                            end
                        end
                        iter = iter+1; dayCount = dayCount + 1;
                    end
                end
            end
        else
            %name column variable name
            iVarNames{iter}  = [behavioralParameters{i},'_',tempFields{m}];
            if numel(tempData2) == 1
                for k = 1:nWorms
                    X(k,iter) =  bpVar{k,size(bpVar,2)}.(tempFields{m});
                end
                iter = iter + 1;
            else %for maxValDay
                dayCount = 0;
                for mmm = 1:11
                    iVarNames{iter} = [behavioralParameters{i},'_',tempFields{m},'_Day',num2str(dayCount)];
                    for k = 1:nWorms
                        X(k,iter) =  bpVar{k,size(bpVar,2)}.(tempFields{m})(mmm);
                    end
                    iter = iter+1; dayCount = dayCount + 1;
                end
            end
        end
    end
end

end

%% getting raw values from movementData cell array into parsable, matrix form with imputable data
% X - ampData (1-2), csData (3-4), dpData (5-6), freqData (7-9), PCData (10)
function [X_raw] = getRawData(ampData,csData,dpData,freqData,PCData)
behavioralParameters = {'ampData';'csData';'dpData';'freqData'};

X_raw = NaN(size(PCData,1),size(PCData,2)-1,10);
iter = 1;

for i = 1:numel(behavioralParameters)
    bpVar = eval(behavioralParameters{i});
    
    tempData = bpVar{1,1};
    tempFields = fieldnames(tempData); 
    if length(tempFields) <3 % amp, cs, and dp data
        for m = 1:length(tempFields)
            for k = 1:size(bpVar,1)
                for j = 1:size(bpVar,2)-1
                    if j <= size(PCData,2)-1
                        X_raw(k,j,iter) =  bpVar{k,j}.(tempFields{m});
                    end
                end
            end
        iter = iter+1;
        end
    else
        for m = 5:length(tempFields) %freq data
            for k = 1:size(bpVar,1)
                for j = 1:size(bpVar,2)-1 %get rid of processed data column
                    if j <= size(PCData,2)-1
                        X_raw(k,j,iter) =  bpVar{k,j}.(tempFields{m});
                    end
                end
            end
            iter = iter+1;
        end
    end
end

% add PC data
X_raw(:,:,10) = cell2mat(PCData(:,1:end-1));
for i = 1:size(X_raw,1)
    for j = 1:size(X_raw,2)
        if X_raw(i,j,10) >= 1000 %censor outliers
            X_raw(i,j,10) = NaN;
        end
    end
end

% adjusting NaNs for freq values
for i = 1:size(X_raw,1)
    for j = 1:size(X_raw,2)
        if all(~isnan(X_raw(i,j,1:6)))
            if any(isnan(X_raw(i,j,7:9)))
                X_raw(i,j,7:9) = 0;
            end
        end
    end
end


end
