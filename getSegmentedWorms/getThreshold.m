function [threshold, image] = getThreshold(image, criteria, startI, startP, effort)
% Finds a threshold that fits 'criteria' or the minimum value of
% 'criteria', in a binary search pattern
% e.g. start at 128 with a power of 4, check 112 and 144, find that 144 
%    is better, check 152 and 136, ...
% image can be 3D, must be uint8
% criteria must be in one of the following formats:
%    {'fillFraction', lowerBound, upperBound} (0< lb < ub < 1)
%         Finds a threshold above which between lb and ub pixels reside
%    {'wormSize', lowerBound, upperBound} (0< lb < ub)
%         Finds a threshold above which a worm is between length lb and ub
%    An anonymous function that takes in (data, thresh) and outputs a 
%         number that is to be minimized
% startI is the starting intensity
% startP is the starting power (of 2), e.g. a startP of 4 means that the
%     function analyzes startI +/- 16
% effort is an integer, higher means a more exhaustive search (up to startP)
%     It's only helpful if the function is smooth-ish
% TODO: add option to min/max a param?

assert(isa(image, 'uint8'))
assert(isa(criteria, 'function_handle') && nargin(criteria) == 2 || ...
    iscell(criteria) && ischar(criteria{1}) && isa([criteria{2:end}],'double'))
if nargin < 5 || isempty(effort)
    effort = startP-1;
end
if nargin < 4 || isempty(startP)
    startP = 32;
end
if nargin < 3 || isempty(startI)
    startI = 128;
end

if isa(criteria, 'function_handle')
    threshold = getThresholdHelper(image, criteria, []);
else
    if strcmpi(criteria{1},'FillFraction')
        objFunc = @(data, thresh) sum(data(:) > thresh) ./ numel(data);
    elseif strcmpi(criteria{1},'WormSize')
        objFunc = @(data, thresh) -1.*getWormLength(data, thresh);
    end
    threshold = getThresholdHelper(image, objFunc, [criteria{2:end}],...
            startI, startP, 1, effort);
end
end

function approximateThresh = getThresholdHelper(frames, evalFunc, bounds,...
    startPoint, power, pixelation, backtrack)
% Finds an answer to the question: "are there any intensities within
%    2^(power) distance of (startPoint) at which evalFunc is within bounds?"
% pixelation can speed things up: only every 2^(pixelation) distances
%    are ever checked, so it can stop earlier 
%    (it will perform ~~2^(power-pixelation+1)-1 evaluations of evalFunc
%       on data, so it may be worthwhile to downsample data beforehand)
% evalFunc takes two parameters, 'data' and 'threshold', and returns 1
%    Make sure it returns numbers that can be inside of [bounds] or be sad
%    Make sure that it increases as the threshold is lowered
% bounds are [lowerBound upperBound] if finding a window, [tolerance]
%    to minimize
if nargin < 7
    backtrack = 4;
end
if nargin < 6
    pixelation = 0;
end
if nargin >= 5
    power = floor(power+2.*eps);
end
if nargin  == 4
    power = min(floor(log2(startPoint)),floor(log2(255-startPoint)));
elseif nargin < 4
    power = 6;
    startPoint = 128;
end
if nargin < 3
    bounds = [];
end
if numel(bounds) == 2
    protocol = 'Window';
elseif numel(bounds) == 1
    protocol = 'Minimize';
elseif isempty(bounds)
    protocol = 'Minimize';
    bounds = 0;
end

if power < pixelation
   approximateThresh = startPoint;
   return
end

thresh = startPoint;
shifts = 2.^(power:-1:pixelation);
history = ones(1, numel(shifts));
approximateThresh = [];
if size(frames,3) > 40
    frames = frames(:,:,1:5:end);
end

for i = 1:numel(shifts)
    funcVal = evalFunc(frames, thresh);
    assert(numel(funcVal) == 1)
    if strcmp(protocol, 'Window')
        if funcVal > bounds(2)
            history(i) = 1;
            thresh = thresh + history(i).*shifts(i);
        elseif funcVal < bounds(1)
            history(i) = -1;
            thresh = thresh + history(i).*shifts(i);
        else
            approximateThresh = thresh;
            return
        end
    elseif strcmp(protocol, 'Minimize')
        funcValL = evalFunc(frames, thresh - shifts(i));
        funcValH = evalFunc(frames, thresh + shifts(i));
        [~, idx] = min([funcValL, funcVal, funcValH]);
        history(i) = idx - 2;
        thresh = thresh + (idx - 2) .* shifts(i);
        if i == numel(shifts)
            approximateThresh = thresh;
        end
    end
    %disp(thresh)
end

if isempty(approximateThresh) || (strcmp(protocol,'Minimize') &&...
        ~isempty(history))
    % starting from the "best-case" scenario, recursively undo the smallest
    %   step (lowest power of two), take the alternate fork, and optimize
    for i = 1:backtrack-pixelation
        newStart = startPoint + sum(history(1:end-i) .* shifts(1:end-i))...
            -history(end-i+1).*shifts(end-i+1);
        newPower = round(log2(shifts(end-i+1)))-1;
        approximateThresh2 = getThresholdHelper(frames, evalFunc,...
            bounds, newStart, newPower, pixelation, min(backtrack, newPower));
        if strcmp(protocol,'Window') && ~isempty(approximateThresh2)
            return
        elseif strcmp(protocol,'Minimize') && evalFunc(frames, ...
                approximateThresh2) < evalFunc(...
                frames,approximateThresh) && ~isempty(approximateThresh)
            approximateThresh = approximateThresh2;
        end
    end
end
end

function wormSize = getWormLength(data, threshold)
vals = getWidthVals(21);
threshedImage = data > threshold;
coloredImage = segWorm(threshedImage,{[2 20]});
connComp = bwconncomp(coloredImage);
wormSize = median(cellfun(@numel,connComp.PixelIdxList));
end
