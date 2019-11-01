function [centerIssue] = checkCenters(image, centers, btol)
% checks whether detected chamber centers are within range (based on image
% intensity values) 
%
% INPUTS:
% image - frame of video of interest
% centers - center location values to test
% btol - tolerance
% OUTPUTS:
% centerIssue - binary error flag

    sql = 37; %10;

    testingChambers = [1,12,49,60];
    lowerLims = zeros(length(testingChambers),2);
    upperLims = zeros(length(testingChambers),2);

    for k=1:length(testingChambers)
        for j=1:2
            if centers(testingChambers(k),j)-sql <1
                lowerLims(k,j) = 1;
            else
                lowerLims(k,j) = round(centers(testingChambers(k),j)-sql);
            end
        end
        if centers(testingChambers(k),1)+sql > 1280
            upperLims(k,1) = 1280;
        else
            upperLims(k,1) = round(centers(testingChambers(k),1)+sql);
        end
        if centers(testingChambers(k),2)+sql > 1024
            upperLims(k,2) = 1024;
        else
            upperLims(k,2) = round(centers(testingChambers(k),2)+sql);
        end
    end

    cChamb1 = image(lowerLims(1,2):upperLims(1,2), lowerLims(1,1):upperLims(1,1));
    cChamb12 = image(lowerLims(2,2):upperLims(2,2),lowerLims(2,1):upperLims(2,1));
    cChamb49 = image(lowerLims(3,2):upperLims(3,2), lowerLims(3,1):upperLims(3,1));
    cChamb60 = image(lowerLims(4,2):upperLims(4,2), lowerLims(4,1):upperLims(4,1));
    
%   checking corner chambers for misalignments in finding centers
    if mean(mean(cChamb1)) <=15 || mean(mean(cChamb12)) <=15 || mean(mean(cChamb49)) <=15 ... 
            || mean(mean(cChamb60)) <=15 || numel(find(cChamb1<15))>=btol || numel(find(cChamb12<15))>=btol ...
            || numel(find(cChamb49<15))>=btol || numel(find(cChamb60<15))>=btol || isempty(cChamb1) || ...
            isempty(cChamb12) || isempty(cChamb49) || isempty(cChamb60) % && i >2
        disp('Issue with automated room locator')
        centerIssue = 1;
    else
        centerIssue = 0;
    end
end