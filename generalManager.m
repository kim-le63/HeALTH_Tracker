%% generalManager.m
% Start here! 
% Handles the data structure holding the output of consecutive videos
% takes in folderName containing videos, outputs data structure seasonalSnapshots
% Responsible for finding chamber locations, maintaining and updating beliefs 
% about the chambers, and revisiting old chambers or reporting errors/conflicting beliefs

%% Load videos of interest
% disp('Select folder with videos of interest')
% folderName = uigetdir(); % select folder with videos of interest
% experimentName = regexp(folderName,filesep,'split'); 
% experimentName = experimentName(end); experimentName = char(experimentName);
% experimentName = experimentName(10:end-15);
% 
% hotelSeasons = dir(fullfile(folderName,'*.avi'));
% hotelSeasons = {hotelSeasons.name};
% [orderedHotelSeasons,ndx] = natsortfiles(hotelSeasons);

seasonalSnapshots  = cell(61,numel(orderedHotelSeasons)); 

%% get chamber center locations for all chambers in device and select chambers to analyze
[seasonalSnapshots, combCenters, chamberIDs, refChamber] = roomLocator(folderName, orderedHotelSeasons, seasonalSnapshots);

%% get lifespan data/ raw movement data
[pxlChange, seasonalSnapshots, deathVid, deathVidName, lifeSpan] = pixelChange(...
    orderedHotelSeasons, seasonalSnapshots, chamberIDs, combCenters, 1, 1, refChamber, 0, [], 1, 50);
% 
% %% get segmented data
seasonalSnapshots = hotelManager(folderName, seasonalSnapshots, orderedHotelSeasons,chamberIDs,37,deathVid);
% 
% %% get behavioral metrics and decline from segmented data
combMovData = getBehavior(seasonalSnapshots,chamberIDs,orderedHotelSeasons);
[freqData,ampData,csData,dpData,areaData] = getMovementReadouts(combMovData, deathVid, chamberIDs);
 
[declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, ...
    declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75,...
    declinePT_duration, relativeHighPeriod_duration] = ...
    getBehaviorParameters(pxlChange, [], 'kmeans', 0, 'pxl change', deathVid);
    PCData = cell(size(pxlChange,1),size(pxlChange,2)+1);
    for k = 1:size(pxlChange,1)
        for j = 1:size(pxlChange,2)
            PCData{k,j} = pxlChange(k,j);
        end
        PCData{k,size(pxlChange,2)+1}.declinePT = declinePT(k); PCData{k,size(pxlChange,2)+1}.relativeHighPeriod = relativeHighPeriod(k);
        PCData{k,size(pxlChange,2)+1}.declineSlope_all = declineSlope_all(k); PCData{k,size(pxlChange,2)+1}.declineInt_all = declineInt_all(k);
        PCData{k,size(pxlChange,2)+1}.declineSlope_high = declineSlope_high(k); PCData{k,size(pxlChange,2)+1}.declineInt_high = declineInt_high(k);
        PCData{k,size(pxlChange,2)+1}.maxValDay = maxValDay(k,:); PCData{k,size(pxlChange,2)+1}.avgValDay = avgValDay(k,:);
        PCData{k,size(pxlChange,2)+1}.maxValDay_high = maxValDay_high(k,:); PCData{k,size(pxlChange,2)+1}.avgValDay_high = avgValDay_high(k,:);
        PCData{k,size(pxlChange,2)+1}.threshVal = threshVal;
        PCData{k,size(pxlChange,2)+1}.declinePT_75 = declinePT_75(k); PCData{k,size(pxlChange,2)+1}.relativeHighPeriod_75 = relativeHighPeriod_75(k);
        PCData{k,size(pxlChange,2)+1}.declinePT_duration = declinePT_duration(k); PCData{k,size(pxlChange,2)+1}.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
    end
    
%% save data 
chamberData.(['chamberData_',experimentName]).combCenters = combCenters;
chamberData.(['chamberData_',experimentName]).refChamber = refChamber;
chamberData.(['chamberData_',experimentName]).chamberIDs = chamberIDs;

lifespanData.(['lifespanData_',experimentName]).deathVid = deathVid;
lifespanData.(['lifespanData_',experimentName]).deathVidName = deathVidName;
lifespanData.(['lifespanData_',experimentName]).lifeSpan = lifeSpan;

movementData.(['movementData_',experimentName]).freqData = freqData;
movementData.(['movementData_',experimentName]).ampData = ampData;
movementData.(['movementData_',experimentName]).csData = csData;
movementData.(['movementData_',experimentName]).dpData = dpData;
movementData.(['movementData_',experimentName]).areaData = areaData;
movementData.(['movementData_',experimentName]).PCData = PCData;

SS.(['seasonalSnapshots_',experimentName]).seasonalSnapshots = seasonalSnapshots;
PC.(['pxlChange_',experimentName]).pxlChange = pxlChange;
CMD.(['combMovData_',experimentName]).combMovData = combMovData;

save(['chamberData_',experimentName,'.mat'],'-struct','chamberData');
save(['lifespanData_',experimentName,'.mat'],'-struct','lifespanData');
save(['movementData_',experimentName,'.mat'],'-struct','movementData');
save(['seasonalSnapshots_',experimentName,'.mat'],'-struct','SS','-v7.3');
save(['pxlChange_',experimentName,'.mat'],'-struct','PC');
save(['combMovData_',experimentName,'.mat'],'-struct','CMD','-v7.3');

