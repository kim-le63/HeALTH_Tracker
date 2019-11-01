function [optimalShift] = motionStabilization(imgA,imgB, shiftRange)
% to offset motion artifacts within a video (due to shaking tubing)

% INPUTS:
% imgA - cropped, filtered initial frame of the reference chamber 
% imgB - cropped, filtered final frame of the reference chamber of the vid
% shiftRange - matrix of offset values of searching for shifted
% frames for both x and y axis
% OUTPUTS: 
% optimalShift - x and y-value offset value for imgB

xShiftVal = shiftRange; 
yShiftVal = shiftRange;

optimalShift = [0, 0];
for i=1:length(xShiftVal)
    for j=1:length(yShiftVal)
        imgDiff = imtranslate(imgB,[xShiftVal(i), yShiftVal(j)])-imgA; %subtracted image
        sumDiff = sum(sum(imgDiff));
        if i==1 && j==1
            optimalDiff = sumDiff;
            optimalShift = [xShiftVal(i),yShiftVal(j)];
        else
            if sumDiff < optimalDiff
                optimalShift = [xShiftVal(i) yShiftVal(j)];
                optimalDiff = sumDiff;
            end
        end
    end
end

end

