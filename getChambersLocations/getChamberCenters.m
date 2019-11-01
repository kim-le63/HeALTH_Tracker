function [centers, centerIssue] = getChamberCenters(image, rad, cropFrame,initRect, initCenters, showCenters)
% get center locations via cross-correlation pattern matching across videos

% INPUTS:
% image - the mxn frame of the video of interest
% rad - radius of chambers
% cropFrame - cropped reference image
% initRect - coordinate info of the cropped reference image
% initCenters - chamber locations in the initial reference image
% showCenters - visualization option for viewing calculated chamber
% locations
% OUTPUTS:
% centers - shifted center locations of chambers
% centerIssue - binary error flag

    c=normxcorr2(cropFrame,image);
    cornerBlocking=c;

    [ypeak, xpeak] = find(cornerBlocking==max(cornerBlocking(:))); 
    
    if isscalar(ypeak)
        yoffSet = ypeak-size(cropFrame,1);
        xoffSet = xpeak-size(cropFrame,2);

        totOff = [xoffSet-initRect(1), yoffSet-initRect(2)]; 
        centers = round(initCenters + totOff);

        if showCenters ==1
            figure (1001); imshow(image); title(['Estimated Chamber Location']);
            viscircles(gca,centers,rad*ones(60,1))
        end

        centerIssue = checkCenters(image, centers, 250);
    else
       [centers,~] = createChamberMask(image,rad); %create mask for chamber array
       centerIssue = 0; 
    end
end

