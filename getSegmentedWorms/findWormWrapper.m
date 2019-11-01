function [newBeliefs] = findWormWrapper(data, thresh, expectation, optOut, oldBeliefs, frames)
% wrapper function to handle passing each frame into findWorm, then
%    integrating the outputs of each layer 
%    analyzing the integrity/stability of the result
% Unlike findWorm, this takes 3D input
% Inputs:
%    data is a (probably 71x71xN) 3D uint8 of a single chamber
%    thresh is a threshold to binarize data with
%    expectation is a struct with fields that interact with findWorm
%        (see findWorm > parseOptions), that fsindWorm uses to decide
%         which binary objects to assess in a given frame

assert(ndims(data) <= 3)
assert(numel(thresh) <= 1)
nFrames = size(data, 3);

if isempty(thresh)
    thresh = oldBeliefs.Intensity;
end

if isempty(frames) %assume analyze all frames
    [wormSize(nFrames), coloredData(:,:,nFrames), ...
    orientationData(:,:,1:2,nFrames), newBeliefs(nFrames)] =...
    findWorm(data(:,:,nFrames), thresh, expectation, optOut, oldBeliefs);

    for i = 1:nFrames
        if i==1
            [wormSize(i), coloredData(:,:,i), orientationData(:,:,1:2,i), newBeliefs(i)]...
            = findWorm(data(:,:,i), thresh, expectation, optOut, oldBeliefs);
        else
            % adding position metric to expectation
            [wormSize(i), coloredData(:,:,i), orientationData(:,:,1:2,i), newBeliefs(i)]...
            = findWorm(data(:,:,i), thresh, expectation, optOut, newBeliefs(i-1));
        end
    end
else
    for i = 1:length(frames)
            if frames(i)==1
                [wormSize(frames(i)), coloredData(:,:,frames(i)), orientationData(:,:,1:2,frames(i)), newBeliefs(frames(i))]...
                = findWorm(data(:,:,frames(i)), thresh, expectation, optOut, oldBeliefs);
            elseif (i~= 1) && (frames(i)-1 == frames(i-1))
                [wormSize(frames(i)), coloredData(:,:,frames(i)), orientationData(:,:,1:2,frames(i)), newBeliefs(frames(i))]...
                = findWorm(data(:,:,frames(i)), thresh, expectation, optOut, newBeliefs(frames(i)-1));
            else
                [wormSize(frames(i)), coloredData(:,:,frames(i)), orientationData(:,:,1:2,frames(i)), newBeliefs(frames(i))]...
                = findWorm(data(:,:,frames(i)), thresh, expectation, optOut, oldBeliefs);
            end
    end
end
% Consolidate each frame's beliefs into one worm
[newBeliefs] = consolidateWorm(wormSize, expectation, newBeliefs);

end


function [modifiedBeliefs] = consolidateWorm(sizes, expects, bels)
% detects any outlying segmented objects in a frame through
% consensus across frames in single video by comparing 1) area and 2)
% position overlap between frames
% outputs new beliefs (intensity, intensitystd, etc.) 
%% Flag Legend
% 0 - no issue
% 1 - too big
% 2 - too small
% 3 - switched positions based on position info from findWorm.m 

[checkFrames, bels] = checkOutliers(sizes, expects, bels);

%flagging potential problem frames
for i = 1:length(sizes)
    if bels(i).PositionError == 1
        bels(i).Flagged = 3;
    else
        bels(i).Flagged = checkFrames(i);
    end
end

modifiedBeliefs = bels;

end


function [checkFrames, beliefs] = checkOutliers(sizes,expectations,beliefs)
% checking for outlying segmented objects across all frames looking at area 
% outputs array listing frames with outlying segmented objects

medianSkelSize = median([beliefs.SkelSize]);

checkFrames = zeros(length(sizes),1); 
% if general area/length of worm is outside of expected worm size range, mark
% all frames in that video for that specific chamber as potentially
% problematic (will check/compare vs. all other chambers)
if (medianSkelSize > (expectations.SkelSize + 1.75*expectations.SkelSizeStd) || ...
        medianSkelSize < (expectations.SkelSize - 1.75*expectations.SkelSizeStd)) 
    if (medianSkelSize > (expectations.SkelSize + 1.75*expectations.SkelSizeStd)) 
        checkFrames = ones(length(sizes),1); % too large
    else
        checkFrames = ones(length(sizes),1)*2; % too small
    end
else
    % checking for issues within individual frames 
    % 'edit' each frame via morphological operations
    se = strel('square',2);  
    
    for ii = 1:length(sizes) %for each frame
        initBW  = beliefs(ii).FullWorm;

        if sum(sum(initBW)) >  medianSkelSize
            eBW = imerode(initBW,se); fBW = bwareafilt(logical(eBW),1);
            BW2 = imdilate(fBW,se);
            newSkelSize = sum(sum(BW2));
            if (newSkelSize >  1.5*medianSkelSize)
                checkFrames(ii) = 1;
            end 
        elseif sum(sum(initBW)) < medianSkelSize
            dBW = imdilate(initBW,se); fBW = bwareafilt(logical(dBW),1);
            BW2 = imerode(fBW,se);
            newSkelSize = sum(sum(BW2));
            if (newSkelSize < 0.5*medianSkelSize)
                checkFrames(ii) = 2; 
            end
        else
            newSkelSize = sum(sum(initBW));
            BW2 = initBW;
        end
        
        beliefs(ii).SkelSize = newSkelSize;
        beliefs(ii).FullWorm = BW2;
    end
end

end

