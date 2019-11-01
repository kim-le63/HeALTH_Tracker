function seasonStats = getSeasonStats(video, seasonStats, chamberIDs,  adjacentSeasonStats, fileName, rad, frames)
% gets segmented image (via findWormWrapper) and marks frames with
% potential errors
% Inputs:
% video - video to analyze
% seasonStats - cell array with chamber location, initial and final frames
% in video
% chamberIDs - array of chambers to analyze
% adjacentSeasonStats - seasonStats info from previous video
% fileName - video name
% rad - radius of chamber
% frames - frames to (re)analyze
% Outputs:
% seasonStats - array with segmented images, intensity threshold values,
% and potential error cases for each frame in the video

supportedExpectations = {'SkelSize', 'Length', 'Width', 'Position', 'WidthA', 'WidthB'};

for jj = 1:length(chamberIDs)
    cCenter = seasonStats{chamberIDs(jj)}.ChamberCenter;
    
    if cCenter(1)-rad < 1 || cCenter(1)+rad > size(video,1) || ...
                cCenter(2)-rad < 1 || cCenter(2) + rad > size(video,2)
        padArray = padMe(cCenter, rad, video, 0);
        vid = padArray;
    else
        vid = video(cCenter(1)-rad:cCenter(1)+rad,...
                cCenter(2)-rad:cCenter(2)+rad,:);
    end
    
    % create circular mask in ROI
    vid = circleMe(vid);

    %% Get threshold details
    intensityStd = 4; %power of 2
    
    %% Get expectations for worm shapes
    if isfield(seasonStats{chamberIDs(jj)}, 'Beliefs') && ~isempty(seasonStats{chamberIDs(jj)}.Beliefs)...
        && ~isempty(intersect(fieldnames(seasonStats{chamberIDs(jj)}.Beliefs), supportedExpectations))
        % There are already specific expectations about this chamber
        expectations = seasonStats{chamberIDs(jj)}.Beliefs;
    elseif isfield(seasonStats{61}, 'Beliefs') && ~isempty(seasonStats{61}.Beliefs)...
            && ~isempty(intersect(fieldnames(seasonStats{61}.Beliefs),...
            supportedExpectations))
        % There are already specific expectations about this video
        expectations = seasonStats{61}.Beliefs;
    else
        % the best threshold must be found
        if nargin >= 4 && ~isempty(adjacentSeasonStats) 
            % There may be expectations about previous videos, but we must check
            if isfield(adjacentSeasonStats{chamberIDs(jj)}, 'Beliefs') && ~isempty(adjacentSeasonStats{chamberIDs(jj)}.Beliefs)...
                && ~isempty(intersect(fieldnames(adjacentSeasonStats{chamberIDs(jj)}.Beliefs), supportedExpectations))
                % There are already specific, reliable expectations about this
                % chamber from the last video
                expectations = adjacentSeasonStats{chamberIDs(jj)}.Beliefs;
            elseif isfield(adjacentSeasonStats{61}, 'Beliefs') && ~isempty(adjacentSeasonStats{61}.Beliefs)...
                    && ~isempty(intersect(fieldnames(adjacentSeasonStats{61}.Beliefs), supportedExpectations))
                        % There are already specific expectations about
                        % this video from the last video
                        expectations = adjacentSeasonStats{61}.Beliefs;
            end
        end
    end
    initVidExpectations = getExpectedSizes(fileName);
    if ~exist('expectations', 'var') || isempty(expectations) || isnan(expectations.SkelSize)
    % Using averaged worm sizes at different video time points [INITIAL GUESS]
       expectations = initVidExpectations;
    end
   
    
    % Checking if expectedd beliefs/sizes make sense; if not, replace with initial guesses
    if expectations.SkelSize > initVidExpectations.SkelSize + 1.5*initVidExpectations.SkelSizeStd ...
            || expectations.SkelSize < initVidExpectations.SkelSize - 1.5*initVidExpectations.SkelSizeStd
        expectations = initVidExpectations;
        disp('change in video expectations')
    end

    intensity = getExpectedIntensity(vid);
    intensityStd = estimateIStd(vid(:,:,30), intensity, expectations);
    
    expectationFields = fieldnames(expectations);
    assert(~isempty(expectationFields));
    for k = 1:numel(expectationFields)
        % Make sure that everything that isn't a request has a Std
        if ~contains(expectationFields{k}, 'Report') && ...
                ~contains(expectationFields{k}, 'Std') && ...
                ~contains(expectationFields{k}, 'Confidence')
            if ~any(cellfun(@(cell) strcmp(cell, [expectationFields{k} ...
                    'Std']), expectationFields))
                if any(cellfun(@(cell) strcmp(cell, [expectationFields{k}...
                        'Confidence']), expectationFields))
                    conf = expectations.([expectationFields{k} 'Confidence']);
                else
                    conf = 1;
                end
                % Normalize by expectation / confidence
                expectations.([expectationFields{k} 'Std']) = ...
                    expectations.(expectationFields{k}) ./ conf;
            end
        end
    end

    if intensityStd > 0
        % Find intensity by finding the threshold that has an object
        %    with a length similar to expected length
        for k = 1:20:size(vid,3)
            warning('off', 'WormHotel:FieldNotReported')
            intensity(k) = getThreshold(vid(:,:,k), @(data, threshold)...
                evalExpectations(findWorm(data, threshold,...
                expectations, []), expectations), intensity, intensityStd, intensityStd-2);
            warning('on', 'WormHotel:FieldNotReported')
        end

        adjInt = intensity(intensity~=0);
        intensity = mean(adjInt);
        intensityStd = estimateIStd(vid(:,:,30), intensity, expectations);
    end
    
    seasonStats{chamberIDs(jj)}.Intensity = intensity;
    seasonStats{chamberIDs(jj)}.IntensityStd = intensityStd;
    
    beliefs.Intensity = intensity;
    beliefs.IntensityStd = estimateIStd(vid(:,:,30), intensity, expectations);
    beliefs.ReportWormObjects = true;
    beliefs.ReportFullWorm = true;
    beliefs.ReportPosition = true;
    beliefs.ReportSize = true;

    frameForChamber = [];
    
    if ~isempty(frames)
        [~,idx]=find(chamberIDs==chamberIDs(jj));
        frameForChamber = frames{idx,1};
    end
    
    [newBeliefs] = findWormWrapper(vid, intensity, expectations, [], beliefs, frameForChamber);

    if isempty(frameForChamber) %first pass
        seasonStats{chamberIDs(jj)}.SegWorm = struct('SegWorm',{newBeliefs.FullWorm});
        seasonStats{chamberIDs(jj)}.Flagged = struct('Flagged',{newBeliefs.Flagged}); 
        seasonStats{chamberIDs(jj)}.wormSize = struct('wormSize',{newBeliefs.SkelSize});%{wormSize});

        if median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]) > initVidExpectations.SkelSize + 1.5*initVidExpectations.SkelSizeStd || ...
        median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]) < initVidExpectations.SkelSize - 1.5*initVidExpectations.SkelSizeStd
            seasonStats{chamberIDs(jj)}.Beliefs = [];
        else
            seasonStats{chamberIDs(jj)}.Beliefs = struct('SkelSize',median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]),'SkelSizeStd',2*std([newBeliefs.SkelSize]));%, ...
        end
        
    else %going through flagged frames individually
        for ii=1:length(frameForChamber)
            seasonStats{chamberIDs(jj)}.SegWorm(frameForChamber(ii)).SegWorm = newBeliefs(frameForChamber(ii)).FullWorm;
            seasonStats{chamberIDs(jj)}.Flagged(frameForChamber(ii)).Flagged = newBeliefs(frameForChamber(ii)).Flagged;
            seasonStats{chamberIDs(jj)}.wormSize(frameForChamber(ii)).wormSize = newBeliefs(frameForChamber(ii)).SkelSize;
        end
        if median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]) > initVidExpectations.SkelSize + 1.5*initVidExpectations.SkelSizeStd || ...
        median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]) < initVidExpectations.SkelSize - 1.5*initVidExpectations.SkelSizeStd
            seasonStats{chamberIDs(jj)}.Beliefs = [];
        else
            seasonStats{chamberIDs(jj)}.Beliefs = struct('SkelSize',median([seasonStats{chamberIDs(jj)}.wormSize.wormSize]),'SkelSizeStd', ...
                2*std([seasonStats{chamberIDs(jj)}.wormSize.wormSize]));
        end
    end
    
end 

end

function intensityStd = estimateIStd(image, intensity, expectation)
% Estimates the uncertainty in a given intensity accurately reflecting 
%    the locations of wormy objects
beliefs.ReportWormObjects = true;
[~,~,~, wormScore1] = findWorm(image, intensity, expectation, [], beliefs);
[~,~,~, wormScore2(1)] = findWorm(image, intensity-4, expectation, [], beliefs);
[~,~,~, wormScore2(2)] = findWorm(image, intensity+4, expectation, [], beliefs);
assert(isfield(wormScore1, 'Worminess'), 'failure in findWorm reporting')
if any(isempty([wormScore1.Worminess, wormScore2(:).Worminess]))
    intensityStd = 255;
else
    intensityStd = 4.*abs(wormScore1.Worminess./min([wormScore2(:).Worminess]));
end
if isempty(intensityStd)
    intensityStd = 255;
end
end
