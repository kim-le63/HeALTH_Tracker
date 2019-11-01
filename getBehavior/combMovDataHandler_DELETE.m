%% combMovDataHandler.m
% File handler for combMovData.mat files to get easily parseable behavioral information 
% (ex. max daily frequency, centroid speed, amplitude, delta pixel change, pixel change value, etc.)

folderName = 'F:\Tracker2018_v2_forKim\Tracker2018_v2_forKim';
PCfolderName = 'F:\pxlChange';

fileList = dir(fullfile(folderName,'combMovData_*.*')); %getting all combMovData files
fileList = {fileList.name};

PCfileList = dir(fullfile(PCfolderName,'pxlChange_*.*')); %getting all combMovData files
PCfileList = {PCfileList.name}; 

%% load deathVids_init variable (deathVids_comb.mat)

h = waitbar(0,'Getting Behavior Data...');

for i =1:length(fileList)
    tic
    fullPath = [folderName '\' fileList{i}]; 
    
    fileToExamine = load(fullPath);
    combMovData = fileToExamine.combMovData;
    dv = deathVids_init{i}.deathTimes;
    [freqData,ampData,csData,dpData,areaData] = getMovementReadouts(combMovData, dv);
    
    PCfullPath = [PCfolderName '\' PCfileList{i}]; 
    PCtoExamine = load(PCfullPath);
    PC=PCtoExamine.pxlChange;
    
    [declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, ...
    declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal] = ...
    getBehaviorParameters(PC, [], 'kmeans', 0, 'pxl change', 0, dv);
    PCData = cell(size(PC,1),size(PC,2)+1);
    for k = 1:size(PC,1)
        for j = 1:size(PC,2)
            PCData{k,j} = PC(k,j);
        end
        PCData{k,size(PC,2)+1}.declinePT = declinePT(k); PCData{k,size(PC,2)+1}.relativeHighPeriod = relativeHighPeriod(k);
        PCData{k,size(PC,2)+1}.declineSlope_all = declineSlope_all(k); PCData{k,size(PC,2)+1}.declineInt_all = declineInt_all(k);
        PCData{k,size(PC,2)+1}.declineSlope_high = declineSlope_high(k); PCData{k,size(PC,2)+1}.declineInt_high = declineInt_high(k);
        PCData{k,size(PC,2)+1}.maxValDay = maxValDay(k,:); PCData{k,size(PC,2)+1}.avgValDay = avgValDay(k,:);
        PCData{k,size(PC,2)+1}.maxValDay_high = maxValDay_high(k,:); PCData{k,size(PC,2)+1}.avgValDay_high = avgValDay_high(k,:);
        PCData{k,size(PC,2)+1}.threshVal = threshVal;
    end
    
    newFileName =['movementData_',fileList{i}(12:end)];
    save(newFileName,'freqData','ampData','csData','dpData','areaData','PCData','-v7.3')
    
    fprintf(['time processing: ' num2str(toc) char(10)])
    waitbar(i/length(fileList))
end

disp('done!')
close(h)