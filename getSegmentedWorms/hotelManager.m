S%% hotelManager.m
% Updates seasonalSnapshots (61 x number of videos cell array) with
% segmented frames, worm size across frames, intensity threshold value for
% segmentation, flagged frames (for potential error cases), and beliefs of
% the worm skeleton size and std for future reference
% 
% Inputs:
% folderName - folder path of videos to analyze
% seasonalSnapshots - (61x number of videos cell array) with chamber locations, 
% initial and final frames
% chamberIDs - array of chambers to analyze
% orderedHotelSeasons - cell array with names of videos to analyze
% rad - radius of chamber
% 
% Outputs:
% seasonalSnapshots - updated seasonalSnapshots cell array with intensity 
% threshold value and std for segmentation, segmented image for each frame, 
% flagged frames (for potential errors), worm size for each frame, and overall
% beliefs from the video (i.e. skeleton size and std) for future reference

function [seasonalSnapshots] = hotelManager(folderName, seasonalSnapshots, orderedHotelSeasons,chamberIDs,rad,deathVids)

h = waitbar(0,'Segmenting Videos...');

for i=1:numel(orderedHotelSeasons)
    fullPath = [folderName '\' orderedHotelSeasons{i}];
    seasonVid = importVideo(fullPath);
    [~,name] = fileparts(fullPath);
    frames = {};

    if i == 1 %manual selection of chambers of interest for the first video of the folder (i.e. single loaded chambers)
        tic
        % attempting 1st pass segmentation w/ background subtraction
        [subVideo, ~] = getBackgroundSubVid (seasonVid, rad);
        [chambersToClean, seasonalSnapshots(:,i)] = concierge (subVideo, {seasonalSnapshots{:,i}}, ...
                chamberIDs,[], name, rad, frames, seasonVid); 
         iter = 0;
            while ~isempty(chambersToClean) && iter < 3         
                for j = 1:length(chambersToClean)
                        flaggedFrames = find([seasonalSnapshots{chambersToClean(j),i}.Flagged.Flagged] ==1);
                        frames(j,1) = {flaggedFrames};
                end
                % segmentation w/o background subtraction for flagged chambers and frames
                [chambersToClean, seasonalSnapshots(:,i)] = concierge (seasonVid, {seasonalSnapshots{:,i}}, ...
                chambersToClean,[], name, rad, frames, seasonVid);
            
                frames={};
                iter = iter+1;
            end
             % adjusting beliefs as needed for next video's reference
             seasonalSnapshots(:,i) = lobbyBoy (seasonVid, {seasonalSnapshots{:,i}}, chamberIDs);

             disp(['Video ' num2str(i) ' done (' num2str(toc) 's)'])   
    else
        % new chamber ID list based on whether video occurs when worm is
        % still alive
        newChamberIDs = [];
        for j=1:length(deathVids)
            if i < deathVids(j)
                newChamberIDs = horzcat(newChamberIDs, chamberIDs(j));
            end
        end
        if isempty(newChamberIDs)
            disp('DONE PROCESSING VIDEOS WITH MOTION')
            return
        end
        
        tic
        % attempting 1st pass segmentation w/ background subtraction
        [subVideo, ~] = getBackgroundSubVid (seasonVid, rad);
        [chambersToClean, seasonalSnapshots(:,i)] = concierge (subVideo, {seasonalSnapshots{:,i}}, ...
                newChamberIDs,{seasonalSnapshots{:,i-1}}, name, rad, frames, seasonVid);
            iter = 0;
            while ~isempty(chambersToClean) && iter < 3 
                for j = 1:length(chambersToClean)
                        flaggedFrames = find([seasonalSnapshots{chambersToClean(j),i}.Flagged.Flagged] ==1);
                        frames(j,1) = {flaggedFrames};
                end
                % segmentation w/o background subtraction for flagged chambers and frames
                [chambersToClean, seasonalSnapshots(:,i)] = concierge (seasonVid, {seasonalSnapshots{:,i}}, ...
                chambersToClean,[], name, rad, frames, seasonVid);
            
                frames={};
                iter = iter+1;
            end
             % adjusting beliefs as needed for next video's reference
             seasonalSnapshots(:,i) = lobbyBoy (seasonVid, {seasonalSnapshots{:,i}}, newChamberIDs);
        disp(['Video ' num2str(i) ' done (' num2str(toc) 's)'])   
    end
    waitbar(i/numel(orderedHotelSeasons))
end                    
end

function [subVideo, background] = getBackgroundSubVid (seasonVid,vidCorrection)
% getting the background subtraction frames of the video
if vidCorrection == 0
    background = backgroundSubtraction(seasonVid);
else
    background = min(seasonVid(:,:,round(size(seasonVid,3)/4)*3:5:end), [], 3);
end

subVideo = seasonVid - background;
edgeMask = edge(background);
se = strel('disk',1);
edgeMask = imdilate(edgeMask,se);
edgeMask = uint8(imcomplement(edgeMask));

subVideo = edgeMask.*subVideo;
end

function [chambersToClean, seasonStats] = concierge (video, seasonStats, ...
     chamberIDs, adjacentSeasonStats, name, rad, frames, initVideo)
% gets segmented images and updated expectations/beliefs (via
% getSeasonStats) and goes through to ID any chambers to resegment them as
% needed (via housekeeping)

seasonStats = getSeasonStats(video, seasonStats, chamberIDs, adjacentSeasonStats, name, rad, frames);
[chambersToClean, seasonStats] = housekeeping(seasonStats, chamberIDs, initVideo);

end

function [seasonStats] = lobbyBoy (video, seasonStats, chamberIDs)
% cleans up after concierge function - modifies beliefs as needed to
% prepare for next video by clearing beliefs for chambers with bad
% segmentation and recalculating overall beliefs (throughout device)

okChambers = chamberIDs;
for j = 1:length(chamberIDs)
    numFlaggedFrames = sum([seasonStats{chamberIDs(j)}.Flagged.Flagged]);
    if numFlaggedFrames > round(0.5*size(video,3)) %marking videos with greater than 50% flagged frames
        seasonStats{chamberIDs(j)}.Beliefs =[];
        okChambers(j)=0;
    end
end
okChambers(okChambers==0)=[];

skelSize = zeros(length(okChambers),1); %wormLength = zeros(length(okChambers),1);  
for k = 1:length(okChambers)
    if ~isempty(seasonStats{okChambers(k)}.Beliefs)
        skelSize(k) = seasonStats{okChambers(k)}.Beliefs.SkelSize;
    end
end
skelSize(skelSize==0)=[];

seasonStats{61}.Beliefs = struct('SkelSize',median([skelSize]),'SkelSizeStd', ...
    2*std([skelSize]));
end


% Checking reasonableness of segmentation via consensus across chambers/frames
function [chambersToClean, seasonStats] = housekeeping( ...
    seasonStats, chamberIDs, initVideo)

numChambers = length(chamberIDs);
numFrames = size(initVideo,3);

numFlaggedFrames = zeros(1,numChambers);
for i = 1:numChambers %for each chamber selected, summing number of flagged frames and removing flagged individual chamber beliefs
    numFlaggedFrames(i) = sum([seasonStats{chamberIDs(i)}.Flagged.Flagged]);
    if numFlaggedFrames(i) > 0.5*size(initVideo,3)
        seasonStats{chamberIDs(i)}.Beliefs =[];
    end
end
if sum(numFlaggedFrames) == numChambers*numFrames
    chambersToClean = chamberIDs;
else
    chambersToClean=chamberIDs(find(numFlaggedFrames));
end  
end
