%% pixelChange.m
% Calculating the pixel difference between the initial and final frame of
% the video as a rough measure of movement and behavior.
%
% Inputs: 
% orderedHotelSeasons - cell array of ordered video names to examine
% seasonalSnapshots - cell matrix with positional information, frames, etc.
% of all chambers of chosen videos
% chamberIDs - array of chambers to analyze
% combCenters - cell array with chamber center locations
% numVidsinDay - number of videos sampled per day
% backgroundMotion - binary indicating whether to account for subtle device
% movements/motion artifacts (assume 1 if not specified)
% refChamber - Number of empty, reference chamber to use as a reference
% point
% checkRef - binary indicating whether to look at previous video to check
% live/dead times (NOTE: only used when dealing with subsampled data)
% completeSeasons - cell matrix with positional information, frames, etc.
% of all chambers of all vidoes (NOTE: only used when dealing with 
% subsampled data)
% tolerance - movement threshold for separating between worm movement and
% noise (assume 1 if not specified)
% intThresh - intensity threshold for separating between worm v. non-worm
% objects (assume 50 if not specified)
%
% Outputs: 
% pxlChange - matrix (# of worms, # of videos) of summed absolute
% difference of pxls between initial & final frame of a video for each worm
% seasonalSnapshots - updated cell matrix with positional information, 
% frames, shift accounting for motion, etc. of all chambers of chosen videos
% deathVid - array of time of death (in video #) for each individual
% deathVid - cell array of name of video of death for each individual
% lifeSpan - array of lifespan (day) for each worm
% completeSeasons - updated cell matrix w/ positional info, frames, etc.
% of all chambers of all vidoes (NOTE: only used when dealing with 
% subsampled data)
% deathVid_ref - array of time of death (in video #) for each individual
% with previous vid check (NOTE: only used when dealing w/ subsampled data)
% deathVid_ref - cell array of name of video of death for each individual
% with previous vid check (NOTE: only used when dealing w/ subsampled data)
% lifeSpan_ref - array of lifespan (day) for each worm with previous vid 
% check (NOTE: only used when dealing w/ subsampled data)


function [pxlChange, seasonalSnapshots, deathVid, deathVidName, lifeSpan, ...
    completeSeasons, deathVid_ref, deathVidName_ref, lifeSpan_ref] = pixelChange(...
    orderedHotelSeasons, seasonalSnapshots, chamberIDs, combCenters, numVidsinDay, ...
    backgroundMotion, refChamber, checkRef, completeSeasons, tolerance, intThresh)
    

deathVid =[];
deathVidName=[];
lifeSpan=[];
completeSeasons=[];
deathVid_ref=[];
deathVidName_ref=[];
lifeSpan_ref=[];


deathVidName_GND=[]; backgroundMotion = 1;

h=waitbar(0,'Please wait...');

if nargin < 11
    intThresh = 50;
    if nargin < 10
        tolerance = 1;
    end
end

numVidsinDay = 4; %edit based on subsampling/recording frequency scheme
numDaysDead = 2;
deathCountThresh = numVidsinDay*numDaysDead;

deathVid = zeros(length(chamberIDs),1);
deathVidName = cell(length(chamberIDs),1);
pxlChange = zeros(length(chamberIDs),length(orderedHotelSeasons));
% movementStatus = zeros(length(chamberIDs),length(orderedHotelSeasons));
% liveDeadStatus = ones(length(chamberIDs),length(orderedHotelSeasons)); % for creating lifespan curves

% getting shift value for each video from initial to final frame (dealing
% with motion artifacts)
if backgroundMotion == 1
    x_yShift = zeros(length(orderedHotelSeasons),2);
    for j=141:length(orderedHotelSeasons)
        [optimalShift] = motionStabilization(seasonalSnapshots{refChamber,j}.fiImg,seasonalSnapshots{refChamber,j}.ffImg, [-10:10]);
        x_yShift(j,:) = optimalShift;
    end
end

% setting up reference videos to validate death (via hourly posture movement)
if checkRef == 1
    deathVid_ref = zeros(length(chamberIDs),1);
    deathVidName_ref = cell(length(chamberIDs),1);
    liveDeadStatus_ref = ones(length(chamberIDs),length(orderedHotelSeasons)); 

    disp('Select reference folder of interest')
    refFolderName = uigetdir(); % select folder with videos of interest

    refHotelSeasons = dir(fullfile(refFolderName,'*.avi'));
    refHotelSeasons = {refHotelSeasons.name};
    [refOrderedHotelSeasons,~] = natsortfiles(refHotelSeasons);

    if isempty(completeSeasons)
        completeSeasons = cell(61,numel(refOrderedHotelSeasons)); 
    end

    rImg = seasonalSnapshots{61,70}.Image;
    [initCenters,~] = createChamberMask(rImg,37); %create mask for chamber array
    figure(1000);
    [cropFrame,initRect]=imcrop(rImg); % crop out the area of the device to use as a reference pattern 
end

% counts number of pixels above a set intensity threshold (intThresh) in 
% the subtracted filtered frames and determines time of death 
for i=1:length(chamberIDs) 
    for j=1:length(orderedHotelSeasons)
        seasonalSnapshots{chamberIDs(i),j}.motionFsubImg = ...
              imabsdiff(imtranslate(seasonalSnapshots{chamberIDs(i),j}.ffImg,[x_yShift(j,:)]),seasonalSnapshots{chamberIDs(i),j}.fiImg);
        seasonalSnapshots{chamberIDs(i),j}.imgShift = x_yShift(j,:);

       [maskedImg] = circleMe(seasonalSnapshots{chamberIDs(i),j}.motionFsubImg); % adding a circular ROI

        % cropping the image (to get rid of potential shift effects)
        columnSum = nnz(~sum(seasonalSnapshots{chamberIDs(i),j}.ffImg,1));
        rowSum = nnz(~sum(seasonalSnapshots{chamberIDs(i),j}.ffImg,2));

        pxlChange(i,j) = sum(sum(maskedImg((2+rowSum+abs(x_yShift(j,2)):end-abs(x_yShift(j,2))-rowSum-1),(2+columnSum+abs(x_yShift(j,1)):end-abs(x_yShift(j,1))-columnSum-1)) >= intThresh)); %binary count. be aware.
    end
    for j=1:size(pxlChange,2)
        if pxlChange(i,j) > tolerance
            movementStatus(i,j) = 1; 
        else
            movementStatus(i,j) = 0;
        end      
    end
end   

% checking for bacterial flooding in 1st row of videos
floodedChambers = chamberIDs(chamberIDs <=12);
floodVids = zeros(size(pxlChange,2),1); %binary check if need to account for bacteria flow in upper chambers
for i=2:size(pxlChange,2)
    shiftVal = seasonalSnapshots{chamberIDs(1),i-1}.imgShift;
    firstVidDiff = imabsdiff(imtranslate(seasonalSnapshots{1,i}.ffImg,[shiftVal]),seasonalSnapshots{1,i}.fiImg);
    floodVids(i-1) = sum(sum(firstVidDiff >=intThresh));
    if sum(sum(firstVidDiff >=intThresh)) > 50
        for j=1:length(floodedChambers) %looking at the first row of the device
            if pxlChange(j,i)>1 && pxlChange(j,i-1)<2 && pxlChange(j,i+1) <2
                pxlChange(j,i)=NaN;
                if movementStatus(j,i) ==1
                    movementStatus(j,i) = 0;
                end
            end
        end
    end
end

stopPoint =size(pxlChange,2);%length(orderedHotelSeasons);

for i=1:length(chamberIDs) 
    movementCount = 0;
    tic
    positionChange = 1;
     for j=stopPoint:-1:deathCountThresh+1
        if movementCount == 0
            if movementStatus(i,j) == 1 && sum(movementStatus(i,j:-1:j-deathCountThresh+1)) >= 2 %looking for movement
                deathTime = j+1;
                deathVid(i) = deathTime;
                deathVidName{i} = orderedHotelSeasons(deathVid(i));
                if checkRef == 1
                    %% validating death (via hourly posture movement)
                    if deathTime >= length(orderedHotelSeasons)-1
    %                         figure(i); subplot(1,3,1); imshow(imtranslate(seasonalSnapshots{chamberIDs(i),deathTime}.ffImg,[seasonalSnapshots{chamberIDs(i),deathTime}.imgShift]));
    %                         subplot(1,3,2); imshow(seasonalSnapshots{chamberIDs(i),deathTime}.fiImg); subplot(1,3,3); imshow(seasonalSnapshots{chamberIDs(i),deathTime}.motionFsubImg);
                        positionChange = 0;
                    end

                    iter = 1;
                    while positionChange == 1 && deathTime < stopPoint %checking if position of the worm is changing from hour to hour 
                        name = orderedHotelSeasons(deathTime);             
                        if contains(name,'Pulse') 
                            vidNum = str2double(name{1,1}(8:end-9));
                        else
                            vidNum = str2double(name{1,1}(8:end-10));
                        end

                        fiImg = seasonalSnapshots{chamberIDs(i),deathTime}.fiImg; 
                        rImg = seasonalSnapshots{refChamber,deathTime}.fiImg;
                        [fiImg2,rImg2,completeSeasons] = checkHourlyVids(refFolderName, refOrderedHotelSeasons, chamberIDs(i), ...
                            refChamber, orderedHotelSeasons{deathTime},completeSeasons, deathTime, combCenters, cropFrame, initRect, initCenters);

                        columnSum = nnz(~sum(fiImg,1)); rowSum = nnz(~sum(fiImg,2));
                        rImg=rImg(1:size(rImg2,1),1:size(rImg2,2));
                        % account for potential video shifts
                        [optimalShift] = motionStabilization(rImg,rImg2,[-5:5]);

                        if isequal(size(fiImg),size(fiImg2))
                            hSubImg = imabsdiff(imtranslate(fiImg2,optimalShift),fiImg); 
                        else
                            hSubImg = imabsdiff(uint8(zeros(size(fiImg,1),size(fiImg,2))),fiImg);
                        end

                        [maskedHSubImg] = circleMe(hSubImg);
                        
                        iter = iter+1;
                        if sum(sum(maskedHSubImg(2+rowSum+abs(optimalShift(2)):end-abs(optimalShift(2))-1-rowSum,...
                                2+columnSum+abs(optimalShift(1)):end-columnSum-1-abs(optimalShift(1))) >= intThresh)) < 50 % worm is dead, no movement
                                positionChange = 0;
                        else
                            deathTime = deathTime + 1;
                        end
                    end
                    deathVid_ref(i) = deathTime;
                    deathVidName_ref{i} = orderedHotelSeasons(deathVid_ref(i));
                    
                    liveDeadStatus_ref(i,deathVid_ref(i):end) = 0;
                end
%                 disp('death reached')
                liveDeadStatus(i,deathVid(i):end) = 0;
                waitbar(i/length(chamberIDs))
                movementCount = 1;
            end
        end
     end
    toc
end

if any(deathVid==0) % worm dying in last video
    deathVid(find(~deathVid))=length(orderedHotelSeasons);
end

lifeCurve = ones(round(length(orderedHotelSeasons)/numVidsinDay),1)*length(chamberIDs);
lifeSpan = zeros(length(chamberIDs),1);
dayCount = 1;

for i=1:numVidsinDay:length(orderedHotelSeasons)+(numVidsinDay-1)
   if sum(any(deathVid==i:i+(numVidsinDay-1)))>0
        lifeCurve(dayCount:end) = lifeCurve(dayCount) - sum(sum(deathVid==i:i+3));
   end
   dayCount = dayCount + 1;
end

for i = 1:length(deathVid)
    lifeSpan(i) = floor((deathVid(i)-2)/numVidsinDay)+1;
end
lifeSpan(lifeSpan==0)=[];

disp(['mean lifespan: ',num2str(mean(lifeSpan)),' stderr: ',num2str(std(lifeSpan)/sqrt(length(deathVid)))])

if checkRef == 1
    if any(deathVid_ref==0) % worm dying in last video
        deathVid_ref(find(~deathVid_ref))=length(orderedHotelSeasons);
    end
    
    lifeCurve_ref = ones(round(length(orderedHotelSeasons)/numVidsinDay),1)*length(chamberIDs);
    lifeSpan_ref = zeros(length(chamberIDs),1);
    dayCount = 1;

    for i=1:numVidsinDay:length(orderedHotelSeasons)+(numVidsinDay-1)
       if sum(any(deathVid_ref==i:i+(numVidsinDay-1)))>0
            lifeCurve_ref(dayCount:end) = lifeCurve(dayCount) - sum(sum(deathVid_ref==i:i+3));
       end
       dayCount = dayCount + 1;
    end

    for i = 1:length(deathVid_ref)
        lifeSpan_ref(i) = floor((deathVid_ref(i)-2)/numVidsinDay)+1;
    end
    lifeSpan_ref(lifeSpan_ref==0)=[];

    disp(['mean ref lifespan: ',num2str(mean(lifeSpan_ref)),' stderr: ',num2str(std(lifeSpan_ref)/sqrt(length(deathVid_ref)))])
end

%% plotting lifespan curve
lifeCurve = ones(round(length(orderedHotelSeasons)/numVidsinDay),1)*length(chamberIDs);

dayCount = 1;
for i=1:numVidsinDay:length(orderedHotelSeasons)+(numVidsinDay-1)
   if sum(any(deathVid==i:i+(numVidsinDay-1)))>0
        lifeCurve(dayCount:end) = lifeCurve(dayCount) - sum(sum(deathVid==i:i+3));
   end
   dayCount = dayCount + 1;
end

figure; plot(lifeCurve/length(chamberIDs)); 
title('Lifespan Curve'); xlabel('Time(days)'); ylabel('Proportion Alive');

if checkRef == 1
    hold on; 
    lifeCurve_ref= ones(round(length(orderedHotelSeasons)/numVidsinDay),1)*length(chamberIDs);

    dayCount = 1;
    for i=1:numVidsinDay:length(orderedHotelSeasons)+(numVidsinDay-1)
       if sum(any(deathVid_ref==i:i+(numVidsinDay-1)))>0
            lifeCurve_ref(dayCount:end) = lifeCurve_ref(dayCount) - sum(sum(deathVid_ref==i:i+3));
       end
       dayCount = dayCount + 1;
    end
    plot(lifeCurve_ref/length(chamberIDs));

end

%% sorting chambers based on time of death
[sortedVid,sortedIdx] = sort(deathVid);

disp('done with getting lifespan and pixel change info!')

end

%% getting initial frame from subsequent video to check for posture changes at a smaller time interval
function [fiImg,rImg,completeSeasons] = checkHourlyVids(refFolderName, refOrderedHotelSeasons, chamberID, refChamber, suspectedDeathVidName, completeSeasons, deathTime, combCenters, cropFrame, initRect, initCenters) 
% for switching in platforms
if contains(suspectedDeathVidName,'Pulse') 
    vidNum = str2double(suspectedDeathVidName(8:end-9));
else
    vidNum = str2double(suspectedDeathVidName(8:end-10));
end
modVidNum = vidNum;

if contains(suspectedDeathVidName,'Pulse') %look at 1 vid after
    deathVidIdx =  find(contains(refOrderedHotelSeasons,strcat('AutoVid',num2str(modVidNum),'Pulse'))); %suspectedDeathVidName));%find index of supposed death video
    refFullPath = [refFolderName '\' refOrderedHotelSeasons{deathVidIdx + 1}];
    nIdx = deathVidIdx + 1;
else
    deathVidIdx =  find(contains(refOrderedHotelSeasons,strcat('AutoVid',num2str(modVidNum),'Steady'))); %suspectedDeathVidName));%find index of supposed death video
    refFullPath = [refFolderName '\' refOrderedHotelSeasons{deathVidIdx + 2}]; % look at next hour video
    nIdx = deathVidIdx + 2;
end
% chamber of interest
[fiImg,rImg, completeSeasons] = getiFrame(refFullPath, chamberID, refChamber, nIdx, completeSeasons, combCenters,deathTime,cropFrame, initRect, initCenters);

end


function [fiImg, rImg, completeSeasons] = getiFrame(fullPath, chamber, refChamber, nIdx, completeSeasons, combCenters, deathTime, cropFrame, initRect, initCenters) %returns initial frame of chambers of interest from inputed video (via full path)
rad = 37; 
video = importVideo(fullPath);

if ~isfield(completeSeasons{61,nIdx},'FullReport') == 1 
    completeSeasons{61,nIdx}.FullPath = fullPath;
    if sum(sum(video(:,:,1))) > 0 && sum(sum((imabsdiff(video(:,:,end),video(:,:,1))))) < 20000000 
        completeSeasons{61,nIdx}.FullReport = 1;
        rCenters = getChamberCenters(video(:,:,1),rad, cropFrame, initRect, initCenters,1);
        [completeSeasons, centerIssue] = getImages (completeSeasons, rCenters, nIdx, video(:,:,1), rad);
        if centerIssue == 1
            % try centers from previous video 
            centerIssue2 = checkCenters(video(:,:,1), combCenters{deathTime}, 250); %CHECK IF CELL WILL WORK!
            if centerIssue2 == 1
                [rCenters,~] = createChamberMask(video(:,:,1),rad); 
            else
                rCenters = combCenters{deathTime};
            end            
        end
        [completeSeasons] = getImages (completeSeasons, rCenters, nIdx, video, rad);
    else
        completeSeasons{61,nIdx}.FullReport = 0;
        for j = 1:60
            completeSeasons{j,nIdx}.Image = uint8(zeros(rad*2+1,rad*2+1));
        end
    end
end

if completeSeasons{61,nIdx}.FullReport == 1
   iImg = completeSeasons{chamber,nIdx}.Image;
   fiImg = imgaussfilt(iImg);
   rImg = completeSeasons{refChamber,nIdx}.Image;
elseif completeSeasons{61,nIdx}.FullReport == 0
   fiImg = uint8(zeros(rad*2+1,rad*2+1));
   rImg = uint8(zeros(rad*2+1,rad*2+1));
end

end


function [completeSeasons] = getImages (completeSeasons, rCenters, nIdx, video, rad)

 for j = 1:60 % for each chamber getting season stats (center + image of frame) 
    completeSeasons{j,nIdx}.ChamberCenter = rCenters(j,2:-1:1);
    center = completeSeasons{j,nIdx}.ChamberCenter;

    if center(1)-rad < 1 || center(1)+rad > size(video,1) || ...
            center(2)-rad < 1 || center(2) + rad > size(video,2)
        padArray = padMe(center, rad, video, 1);
        completeSeasons{j,nIdx}.Image = padArray;
    else
        completeSeasons{j,nIdx}.Image = video(center(1)-rad:center(1)+rad,...
        center(2)-rad:center(2)+rad); 
    end
 end
 
end

% adding a circular ROI
function [maskedImg] = circleMe(image)

chamberMat = false(size(image,1), size(image,2));
[r2,c2]=meshgrid(1:size(image,2),1:size(image,1)); 
sMask =(((r2-round(size(image,2)/2)).^2+(c2-round(size(image,1)/2)).^2)<=(35)^2); 
chamberMat = chamberMat | sMask;
maskedImg = image.*uint8(chamberMat);         
end