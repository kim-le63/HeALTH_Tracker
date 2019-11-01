function [combMovData] = getBehavior(seasonalSnapshots,chamberIDs,orderedHotelSeasons)
%% Output Data
%combMovData is a cell array of structures (# of vids x 1) with each cell 
%in format movData{chamberIndex,frameIndex}.PROPERTY
%where Chamber=61 corresponds to the entire frame
%and Frame=max(Frame)+1 corresponds to the chamber across all frames
%Valid properties of movData FOR VALID CHAMBER/FRAME are:
    %segmentError: logical of whether there was an error in processing
    %tangents: 25x1 double of tangent values
    %rgLine: 100x5 double of different data about the worm bodyline
        %1st/2nd col: x/y data
        %3rd col: distance in pixels along line
        %4th/5th col: smoother x/y data
    %centroid: position of worm centroid w/respect to the chamber circle
    %image: image of the chamber circle
    %centroidVelocity/blVelocity: two ways of measuring worm velocity
        %[change in rows; change in columns]
    %centroidSpeed/blSpeed: scalars of aforementioned velocities
    %amplitude: amplitude of the worm
    %area: area of the worm
    %deltaPix: (pixels gained since last frame) + (pixels lost since last)
%%

numVids = size(seasonalSnapshots,2);
 
combMovData = cell(numVids,1);
numSubsections = 25;

h = waitbar(0,'Getting Movement Data...');

for n =1:numVids
    tic
    findNumFrames = 0;
    for ii = 1:length(chamberIDs)
       if isfield(seasonalSnapshots{chamberIDs(ii),n},'SegWorm')
            findNumFrames = 1;
            numFrames = size(seasonalSnapshots{chamberIDs(ii),n}.SegWorm,2);
       end
    end
    if findNumFrames == 0
       return % no more worms to segment/get behavioral info
    end
    
    name = orderedHotelSeasons(n);             
    if contains(name,'Pulse') 
        vidNum = str2double(name{1,1}(8:end-9));
    else
        vidNum = str2double(name{1,1}(8:end-10));
    end
    % setting size expectations for final segmentation check
    if vidNum <= 25
        wormMinSize = 80; wormMaxSize = 250;
    elseif vidNum <=40
        wormMinSize = 150; wormMaxSize = 300;
    elseif vidNum <=80
        wormMinSize = 220; wormMaxSize = 400;
    else
        wormMinSize = 250; wormMaxSize = 500;
    end
    for j = 1:length(chamberIDs)
        iChamber = chamberIDs(j);
 
        for iFrame = 1:numFrames
            % Set initial parameters    
            combMovData{n}.movData{iChamber,iFrame}.segmentError = false;
            combMovData{n}.movData{iChamber,iFrame}.tangents = zeros(numSubsections,...
                1,'double');
            combMovData{n}.movData{iChamber,iFrame}.rgLine = zeros(100,6);
            combMovData{n}.movData{iChamber,iFrame}.centroid = [0,0];
            combMovData{n}.movData{iChamber,iFrame}.area = NaN;
            combMovData{n}.movData{iChamber,iFrame}.eccentricity = NaN;
            if isfield(seasonalSnapshots{iChamber,n},'SegWorm')
                image = seasonalSnapshots{iChamber,n}.SegWorm(iFrame).SegWorm;
                stats = regionprops(image,'Centroid','Area','Eccentricity');
                combMovData{n}.movData{iChamber,iFrame}.image = image;

                sizeCheck = 0;
                if sum(sum(image)) < wormMinSize || sum(sum(image)) > wormMaxSize
                    sizeCheck = 1;
                end

                if ~isempty(stats) && seasonalSnapshots{iChamber,n}.Flagged(iFrame).Flagged ~= 3 && seasonalSnapshots{iChamber,n}.Flagged(iFrame).Flagged ~= 4 && sizeCheck == 0

                    combMovData{n}.movData{iChamber,iFrame}.centroid = stats(1).Centroid;
                    combMovData{n}.movData{iChamber,iFrame}.area = stats(1).Area;
                    combMovData{n}.movData{iChamber,iFrame}.eccentricity = stats(1).Eccentricity;

                    [rgLine, tangents, dStatus, Amp] = getLine(image,...
                    numSubsections);

                    combMovData{n}.movData{iChamber,iFrame}.centroidVelocity = [NaN;NaN];
                    combMovData{n}.movData{iChamber,iFrame}.centroidSpeed = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.blVelocity = [NaN;NaN];
                    combMovData{n}.movData{iChamber,iFrame}.blSpeed = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.deltaPix = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.pxlChange = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.amplitude = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.tangents = NaN(numSubsections,1);

                     %if getLine fails on frame, skip ahead to next frame
                    if dStatus == -1 
                        combMovData{n}.movData{iChamber,iFrame}.segmentError = true;
                    else
                        padsize = 100 - length(rgLine); %% ADJUST BASED ON VIDEO NUMBER!
                        if padsize <0
                            combMovData{n}.movData{iChamber,iFrame}.segmentError = true;
                        else
                            rgLine = padarray(rgLine, [padsize 0], 'post');
                            combMovData{n}.movData{iChamber,iFrame}.rgLine = rgLine;
                            combMovData{n}.movData{iChamber,iFrame}.tangents = tangents;
                            combMovData{n}.movData{iChamber,iFrame}.amplitude = Amp;
                            cent1 = [0,0];
                            line1 = [0,0];
                            if iFrame ~= 1 && ~combMovData{n}.movData{iChamber,iFrame-1}.segmentError...
                                    && ~isequal(sum(sum(combMovData{n}.movData{iChamber,iFrame-1}.rgLine)),0)
                                cent1 = combMovData{n}.movData{iChamber,iFrame-1}.centroid;
                                line1 = combMovData{n}.movData{iChamber,iFrame-1}.rgLine(:,4:5);
                            end
                            cent2 = combMovData{n}.movData{iChamber,iFrame}.centroid;
                            line2 = combMovData{n}.movData{iChamber,iFrame}.rgLine(:,4:5);
                            if ~isequal(cent1,[0,0]) && ~isequal(cent2,[0,0])
                                combMovData{n}.movData{iChamber,iFrame}.centroidSpeed =...
                                    sqrt((cent1(1)-cent2(1))^2+(cent1(2)-cent2(2))^2);
                                combMovData{n}.movData{iChamber,iFrame}.centroidVelocity =...
                                    [cent2(1)-cent1(1);cent2(2)-cent1(1)];
                                thisFrame = combMovData{n}.movData{iChamber,iFrame}.image;
                                lastFrame = combMovData{n}.movData{iChamber,iFrame-1}.image;
                                combMovData{n}.movData{iChamber,iFrame}.deltaPix =...
                                    (sum(sum(thisFrame | lastFrame)) -...
                                    sum(sum(thisFrame & lastFrame))) /...
                                    (sum(sum(thisFrame | lastFrame)));
                                combMovData{n}.movData{iChamber,iFrame}.pxlChange =...
                                    sum(sum(imabsdiff(thisFrame,lastFrame)));
                            end
                            if ~isempty(line1) && ~isempty(line2) &&...
                                    ~isequal(sum(sum(line1)),0) &&...
                                    ~isequal(sum(sum(line2)),0)
                                combMovData{n}.movData{iChamber,iFrame}.blVelocity=...
                                    getBodyLineVelocity(line1,line2);
                                if abs(combMovData{n}.movData{iChamber,iFrame}.blVelocity(1)) > 5 ||...
                                    abs(combMovData{n}.movData{iChamber,iFrame}.blVelocity(2)) > 5
                                    combMovData{n}.movData{iChamber,iFrame}.blVelocity = [-10;-10];
                                end
                                combMovData{n}.movData{iChamber,iFrame}.blSpeed=...
                                    sqrt(combMovData{n}.movData{iChamber,iFrame}.blVelocity(1)^2+...
                                    combMovData{n}.movData{iChamber,iFrame}.blVelocity(2)^2);

                            end
                        end
                    end
                else
                    combMovData{n}.movData{iChamber,iFrame}.segmentError = true;
                    combMovData{n}.movData{iChamber,iFrame}.centroid = [NaN;NaN];
                    combMovData{n}.movData{iChamber,iFrame}.area = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.eccentricity = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.centroidVelocity = [NaN;NaN];
                    combMovData{n}.movData{iChamber,iFrame}.centroidSpeed = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.blVelocity = [NaN;NaN];
                    combMovData{n}.movData{iChamber,iFrame}.blSpeed = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.deltaPix = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.amplitude = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.pxlChange = NaN;
                    combMovData{n}.movData{iChamber,iFrame}.tangents = NaN(numSubsections,1);
                end
            else
                combMovData{n}.movData{iChamber,iFrame}.segmentError = NaN;
                combMovData{n}.movData{iChamber,iFrame}.centroid = [NaN;NaN];
                combMovData{n}.movData{iChamber,iFrame}.area = NaN;
                combMovData{n}.movData{iChamber,iFrame}.eccentricity = NaN;
                combMovData{n}.movData{iChamber,iFrame}.centroidVelocity = [NaN;NaN];
                combMovData{n}.movData{iChamber,iFrame}.centroidSpeed = NaN;
                combMovData{n}.movData{iChamber,iFrame}.blVelocity = [NaN;NaN];
                combMovData{n}.movData{iChamber,iFrame}.blSpeed = NaN;
                combMovData{n}.movData{iChamber,iFrame}.deltaPix = NaN;
                combMovData{n}.movData{iChamber,iFrame}.amplitude = NaN;
                combMovData{n}.movData{iChamber,iFrame}.pxlChange = NaN;
                combMovData{n}.movData{iChamber,iFrame}.tangents = NaN(numSubsections,1);
            end
        end
 
    end
    fprintf(['time processing: ' num2str(toc) char(10)])
    waitbar(n/numVids)

end
disp('Captured all movement!')
close(h)

for n= 1:numVids
    numFrames = size(seasonalSnapshots{chamberIDs(1),n}.SegWorm,2);
    for j = 1:length(chamberIDs)
        iChamber = chamberIDs(j);
        %Index out all the centroid velocities of each chamber for graphing
        for i = 1:numFrames
            combMovData{n}.movData{iChamber,numFrames+1}.centroid(:,i) =...
                combMovData{n}.movData{iChamber,i}.centroid;
            combMovData{n}.movData{iChamber,numFrames+1}.centroidVelocity(:,i) =...
                combMovData{n}.movData{iChamber,i}.centroidVelocity;
            combMovData{n}.movData{iChamber,numFrames+1}.centroidSpeed(i) =...
                combMovData{n}.movData{iChamber,i}.centroidSpeed;
            combMovData{n}.movData{iChamber,numFrames+1}.blVelocity(:,i) =...
                combMovData{n}.movData{iChamber,i}.blVelocity;
            combMovData{n}.movData{iChamber,numFrames+1}.blSpeed(i) =...
                combMovData{n}.movData{iChamber,i}.blSpeed;
            combMovData{n}.movData{iChamber,numFrames+1}.deltaPix(i) =...
                combMovData{n}.movData{iChamber,i}.deltaPix;
            combMovData{n}.movData{iChamber,numFrames+1}.Amplitude(:,i) =...
                combMovData{n}.movData{iChamber,i}.amplitude;
            combMovData{n}.movData{iChamber,numFrames+1}.tangents(:,i) =...
                combMovData{n}.movData{iChamber,i}.tangents;
        end
    end
end

end