function [rgLine, Tangent, dStatus, Amplitude] = getLine(rgBody, numSubsections, hPts)
%% Get vector of points in line in image
% input:    rgBody:         binarized image of worm  
%           hPts:           the position of the head in  the previous frame
%           numSubsections: number of even segments to split skeleton into
% output: Tangent: column vector of normalized tangent directions for numSubsections
%           evenly spaced pts along the skeleton
% rgLine:
%           1st column: x
%           2nd column: y
%           3rd column: distance in pixels along the line (2 norm)
%           4th column: smoothed x
%           5rd column: smoothed y
% Tangent: vector with tangent angles, sampled at numSubsections (the mean 
%          has been subtracted so this vector has zero mean)
% Amplitude: amplitude of body bend
% dStatus: error indicator

Amplitude = []; Tangent = NaN(numSubsections,1);
if nargin ==2
    hPts = [-1, -1];
else
    fprintf('Too many input arguments')
end

rgBody = padarray(rgBody, [2, 2]);
I = bwmorph(rgBody,'thin',Inf);

%% smoothing
iConn = bwconncomp(I);
numPix = cellfun(@numel,iConn.PixelIdxList);
[~,maxInd] = max(numPix);
rgThin = false(size(I));
if sum(sum(I)) > 0
    rgThin(iConn.PixelIdxList{maxInd}) = true;
end
rgThin = SimpleEndpoints(rgThin);

try
      
    iPoints = sum(sum(I));
    if iPoints< numSubsections %ignore the image when the size of the object is too small --by zym
        rgLine = NaN(numSubsections +1, 5);
        dStatus = -1;
        Amplitude = [];
        return;
    end
    
    rgEnds = uint8(rgThin);    
    rgEnds = imfilter(rgEnds, [1,1,1;1,10,1;1,1,1],'same');
    [yEnds, xEnds] = find(rgEnds==11);
    
    if isempty(yEnds)
        dStatus = -1;
        rgLine = NaN(numSubsections +1, 5);
        return
    else
            %Add first point
        rgLine = zeros(iPoints,5);    
        rgLine(1,1:2) = [xEnds(1), yEnds(1)];
        rgLine(1,3) = 0;

        for i = 2:iPoints
    
            lastPoint = rgLine(i-1,1:2);
            rgThin(lastPoint(2),lastPoint(1)) = 0;

            %Get points around last point   
            rgDirection = rgThin(lastPoint(2)-1:lastPoint(2)+1, ...
                lastPoint(1)-1:lastPoint(1)+1).* ...
                [1,2,1;2,0,2;1,2,1];     

            %Find next point 
            v = max(rgDirection(:));
            [yLoc, xLoc] = find(rgDirection==v);
            y = yLoc + rgLine(i-1,2) - 2 ;
            x = xLoc + rgLine(i-1,1) - 2;
            
            if length(yLoc)==1
                idx = 1;
            else % when there is a fork, we walk along the previous direction --by zym
                prev = [rgLine(i-1,1)-rgLine(i-2,1) rgLine(i-1,2)-rgLine(i-2,2)];
                cur = [x-rgLine(i-1,1) y-rgLine(i-1,2)];
                len = sqrt(sum(cur.*cur,2));
                cur = cur./(len*ones(1,2));
                projection = cur*prev';
                [~, idx] = max(projection);           
            end

            %Add point to line 
            rgLine(i,1) = x(idx);
            rgLine(i,2) = y(idx);
        end
    end

catch
    try
        
        [yPix,xPix]=find(rgThin);
        rgP=[xPix, yPix];
        ciDist = [];
        
        %choose new first point of line based on where head was previously 
        if (isempty(yEnds) && hPts(1) ~= -1) || rgThin(yEnds(1),xEnds(1)) == false
            for ci=1:length(yPix)
                ciDist=[ciDist;pdist([hPts;rgP(ci,:)])];
            end
            [~,headIdx]=min(ciDist);
            yEnds(1)=yPix(headIdx);
            xEnds(1)=xPix(headIdx);
        end
      
    catch
        rgLine = NaN(numSubsections +1, 5);
        Amplitude = [];
        dStatus = -1;
        return;
    end
end

rgLine(:,4) = smooth(rgLine(:,1), 8); %smoothed x vals
rgLine(:,5) = smooth(rgLine(:,2), 8); %smoothed y vals

% calculate distance from beginning of curve at each point 

rgLine(1,3) = 0;

for i = 2:iPoints
    %approximate curve length
    try
    rgLine(i,3) = rgLine(i-1,3) + pdist([rgLine(i, 4:5); rgLine(i-1, 4:5)]);
    catch
    end
end

Amplitude = getAmplitude(rgLine(:,4),rgLine(:,5),.09);
backbone = rgLine(:,1:2); 
if ~isempty(find(~backbone))
    backbone = backbone(1:find(backbone(:,1)==0,1)-1,:);
end
Tangent = calculateTangentsFromBackbone(backbone,numSubsections);

dStatus = 1;
    
end


