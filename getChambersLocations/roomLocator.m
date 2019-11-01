%% roomLocator.m
% Identifies worms to analyze in the device and gets chamber center
% locations for all files of the device in the chosen folder via
% cross-correlation to perform automated chamber identification. User
% selects chamber locations for initial reference frame, along with
% chambers to analyze. The user also selects an empty, reference chamber,
% and crops a 'static' pattern of interest (ex. a device feature that will
% not change drastically over time) for subsequent pattern matching across
% time.
%
% INPUTS: 
% folderName - file path (selected in generalManager.m) 
% orderedHotelSeasons - cell array (1 x number of vids in folder) with
% names of video files (obtained in generalManager.m) 
% seasonalSnapshots - empty cell array (61 x number of vids in folder)
% (initialized in generalManager.m)
% OUTPUTS:
% seasonalSnapshots - updated seasonalSnapshots array with cropped images
% of the initial and final frames of each video (initImg, finImg), 
% gaussian filtered images (fiImg, ffImg), and chamber center
% locations (ChamberCenter). 61 row has information for the reference
% chamber, the background, and a snapshot of the video
% combCenters - cell array (number of vids x 1) with 60x2 matrix of center
% coordinates for each chamber for each video 
% chamberIDs - array of chamber IDs to analyze
% refChamber - user selected empty chamber to use for subsequent
% registration across videos

function [seasonalSnapshots, combCenters, chamberIDs, refChamber] = roomLocator(...
    folderName, orderedHotelSeasons, seasonalSnapshots)

rad = 37; %chamber radius

% Initializing main variables
combCenters = cell(numel(orderedHotelSeasons),1);
h = waitbar(0,'Getting Chamber Centers...');

% Getting reference image for pattern matching
fullPath = [folderName '\' orderedHotelSeasons{5}]; 
seasonVid = importVideo(fullPath); 
background = backgroundSubtraction(seasonVid);

rImg = background(:,:,end);
[initCenters,~] = createChamberMask(rImg,rad);
close;

disp('select worms to analyze')
[chamberIDs] = getOccupiedRooms(seasonVid(:,:,11), initCenters);
disp('select empty reference chamber')
[refChamber] = getOccupiedRooms(seasonVid(:,:,11), initCenters);

figure(1000); title('crop reference pattern')
[cropFrame,initRect]=imcrop(rImg); % crop out the area of the device to use as a reference pattern
close;

errorCounter = 0;

for i= 1:numel(orderedHotelSeasons)
    fullPath = [folderName '\' orderedHotelSeasons{i}];
    [~,name] = fileparts(fullPath);

    if contains(name,'Pulse') 
        vidNum = str2double(name(8:end-5));
    else
        vidNum = str2double(name(8:end-6));
    end

    seasonalSnapshots{61,i}.FullReport = true;
    seasonalSnapshots{61,i}.FullPath = fullPath;
    try
        seasonVid = importVideo(fullPath);
        seasonalSnapshots{61,i}.VideoExists = true;
    catch
        disp(['Video load error in ' fullPath]);
        seasonalSnapshots{61,i}.VideoExists = false;
        continue
    end

    disp(['Processing video ' name])

    background = backgroundSubtraction(seasonVid);
%       Using cross-correlation to track video shifts 
    [centers, centerIssue] = getChamberCenters(seasonVid(:,:,end), rad, cropFrame,initRect, initCenters, 1);

    if i>1 && centerIssue ==1 % try with centers from prior video
        errorCounter = errorCounter+1;
        centerIssue2 = checkCenters(seasonVid(:,:,end), combCenters{i-1}, 250);
        if centerIssue2 == 1
            if errorCounter >=2
                figure(1000); title('crop reference pattern')
                [cropFrame,initRect]=imcrop(background); % crop out the area of the device to use as a reference pattern
                close;
                [centers,~] = createChamberMask(background,rad); %create mask for chamber array
                close;
                errorCounter = 0;
            else
            [centers,~] = createChamberMask(background,rad); %create mask for chamber array
            close;
            end
        else
            figure (1001); imshow(seasonVid(:,:,end)); title(['Estimated Chamber Location']);
            viscircles(gca,combCenters{i-1},rad*ones(60,1))
            centers = combCenters{i-1};
        end
    elseif centerIssue == 1
        errorCounter = errorCounter+1;
        if errorCounter >=4
            figure(1000); title('crop reference pattern')
            [cropFrame,initRect]=imcrop(background); % crop out the area of the device to use as a reference pattern
            close;
            [centers,~] = createChamberMask(background,rad); %create mask for chamber array
            close;
            errorCounter = 0;
        else
            [centers,~] = createChamberMask(background,rad); %create mask for chamber array
            close;
        end
    end
    
    if centerIssue == 0
        errorCounter = 0;
    end

    for j = 1:60
        seasonalSnapshots{j,i}.ChamberCenter = centers(j,2:-1:1);
        seasonalSnapshots{j,i}.initImg = getCroppedImg(seasonalSnapshots{j,i}.ChamberCenter, rad, seasonVid(:,:,round(size(seasonVid,3)/2)));
        seasonalSnapshots{j,i}.finImg = getCroppedImg(seasonalSnapshots{j,i}.ChamberCenter, rad, seasonVid(:,:,end));
        seasonalSnapshots{j,i}.fiImg = imgaussfilt(seasonalSnapshots{j,i}.initImg);
        seasonalSnapshots{j,i}.ffImg = imgaussfilt(seasonalSnapshots{j,i}.finImg);
    end
    seasonalSnapshots{61,i}.refChamber = refChamber;
    seasonalSnapshots{61,i}.refChamberCenter = centers(refChamber,2:-1:1);
    seasonalSnapshots{61,i}.rInitImg = getCroppedImg(seasonalSnapshots{refChamber,i}.ChamberCenter, rad, seasonVid(:,:,round(size(seasonVid,3)/2)));
    seasonalSnapshots{61,i}.rFinImg = getCroppedImg(seasonalSnapshots{refChamber,i}.ChamberCenter, rad, seasonVid(:,:,end));

    seasonalSnapshots{61,i}.Image = seasonVid(:,:,end);
    seasonalSnapshots{61,i}.Background = background;

    combCenters{i} = centers;
    waitbar(i/numel(orderedHotelSeasons))
end

disp('Found all chambers')
close(h)

end

function [chamberIDs] = getOccupiedRooms(Img, centers)
position = centers; %annotation check for chamber number
position(1,:) = centers(1,:)-10;
position(:,2) = centers(:,2)-50;
value = 1:60;
Image = insertText(Img, position, value);

figure(1000), imshow(Image), title('Click on the chamber to analyze')
[r,c] = ginput;
hold on;

tempCenters=zeros(12,2,5);
chamberIDs=zeros(1,length(r));
for j=1:length(r)
    k=1;
    for i=1:12:60
        tempCenters(:,:,k) = centers(i:i+11,:);
        k=k+1;
    end
    avgR=mean(tempCenters(:,2,:));
    [~,idx]=min(abs(avgR-c(j))); %finding matching row of selected worm
    [~,idx2]=min(abs(tempCenters(:,1,idx)-r(j)));
    chamberIDs(j)=(idx-1)*12+idx2;
end

close;
end