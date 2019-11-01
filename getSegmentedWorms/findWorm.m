function [wormSize, coloredData, orientationData, beliefs] = findWorm(data, thresh,...
    expectation, optOut, beliefs)
% thresholds 'data' by 'thresh', then evaluates which of the binary
%    objects most closely matches the criteria of expectedSize or
%    expectedPos
% Inputs:
% data - image of a single frame of a chamber
% thresh - uint8 to threshold the image
% optOut - cell with any of 'Length', 'Width', or 'SkelSize' that determines
%    whether worm length (=skeleton size),
%    width (=median width) or size (=skeleton size * skeleton width) 
%    are reported as wormSize
% expectation  - handles how the function chooses what binary objects are worms:
%    findWorms chooses objects closer to it (it's formatted as a structure 
%    with the fieldnames as e.g. 'Length' and values as numbers)
% expectedPos is a list of linear pixel indices whose overlap is maximised
% Outputs:
% wormSize - size of worm
% coloredData - color coded worm to indicate thickness
% beliefs - updated beliefs structure

assert(ismatrix(data) && min(size(data)) > 1, ['Image input required:'...
    ' Pass video data to findWormWrapper']);
assert(numel(thresh) == 1)

if nargin < 4 || isempty(optOut)
    if nargin > 2 && ~isempty(expectation)
        expectedFields = fieldnames(expectation);
        optOut = expectedFields((cellfun(@(cell) ~contains(cell, 'Std') && ...
            ~contains(cell, 'Confidence'), expectedFields)) == 1);
    else
        optOut = {'Length'};
    end
    if isempty(optOut)
        optOut = {'Length'};
    end
end
if nargin < 3 || isempty(expectation)
    warning('No expectations given')
end
if nargin < 2 || isempty(thresh)
    assert(islogical(data), 'Data must be binarized')
else
    data = data > thresh;
end
if isfield(expectation, 'WidthRange') && ~isempty(expectation.WidthRange)
    if numel(expectation.WidthRange) == 2
        widths = expectation.WidthRange;
    elseif numel(expectation.WidthRange) == 1
        widths = [expectation.WidthRange expectation.WidthRange];
    end
else
    if isfield(expectation, 'SkelSize') && expectation.SkelSize < 200
            widths = [2 4];
    elseif isfield(expectation, 'SkelSize') && expectation.SkelSize < 300
        widths = [2 6];
    else
        % The range of widths of the worm body, in pixels, that segWorm looks for
        widths = [2 8];
    end
end
fcnList = parseOptions(optOut, expectation);
[coloredData, orientationData] = segWorm(data,{widths});
anySkel = bwconncomp(coloredData > 0);
anySkel.PixelIdxList([cellfun(@numel, anySkel.PixelIdxList) < 10]) = [];
anySkel.NumObjects = numel(anySkel.PixelIdxList);

% size based evaluation of objects
[idx, val] = getExpect(fcnList, expectation, coloredData,...
    anySkel.PixelIdxList);

% checking w/ position info 
if nargin==5 && isfield(beliefs, 'Position') && ~isempty(beliefs.Position) && ~isempty(idx)
    positionOverlap = zeros(length(idx),1);
    
    checkPosition = @(cell, expectedIdx)...
        intersect([cell(:)],expectedIdx)./numel(expectedIdx);
    
    for i=1:length(idx)
        positionOverlap(i) = size(checkPosition(anySkel.PixelIdxList{1,i},beliefs.Position),1);
    end
    if positionOverlap(idx(1)) == 0 % making sure worm has some overlap with previous frame
        numOverlap = find(positionOverlap~=0);
        if isempty(numOverlap)
            positionError = 1;
            selectedObjID = 1;
        elseif isempty(find(numOverlap == idx(2),1))
            positionError = 1;
            selectedObjID = 1; 
        else
            selectedObjID = find(numOverlap==idx(2),1);
            positionError = 0;
        end
    else
        selectedObjID = 1;
        positionError = 0;
    end
        
else
    selectedObjID = 1; % assume no position error, going with object with the 'wormiest' score
    positionError = 0;
end


%% Setting wormSize
nonFunctional = setdiff(optOut, fieldnames(fcnList));
if iscell(nonFunctional) && (numel(nonFunctional) ~= 0 || ~isempty(nonFunctional))
    warning(['Options ' strjoin(nonFunctional) ' had no implementation'])
end

optOut = intersect(fieldnames(fcnList), optOut);
if ~isempty(anySkel.PixelIdxList) 
    wormSkel = anySkel.PixelIdxList{idx(selectedObjID)}; %only wormiest obj.
    % Report sizes
    for i = 1:numel(optOut)
        opVal = fcnList.(optOut{i})(wormSkel, coloredData);
        wormSize.(optOut{i}) = opVal; 
    end
else
    wormSkel = 0;
    for i = 1:numel(optOut)
        wormSize.(optOut{i}) = 0;
    end
end

% Update beliefs if requested in the call to findWorm
if nargout > 3
    if nargin < 5
        beliefs = struct();
    end
    if isfield(beliefs, 'ReportWormObjects') && beliefs.ReportWormObjects
        % Report anything at least half as wormy as the best-case
        if isempty(val)
            beliefs.Worminess = [];
            beliefs.WormObjects = [];
        else
            beliefs.Worminess = val(val <= val(selectedObjID) .* 2);
            beliefs.WormObjects = anySkel.PixelIdxList(idx(val <= (val(selectedObjID) .* 2)));
        end
    end
    if isfield(beliefs, 'WormThreshold') && ~isempty(beliefs.WormThreshold)
        % Potential alternative method of determining worm presence 
        beliefs.NumberOfWorms = sum(val > beliefs.WormThreshold);
    end
    beliefs.Intensity = thresh;
    if isfield(beliefs,'ReportColoredData') && beliefs.ReportColoredWorm
        beliefs.ColoredWorm = coloredData;
    end
    if isfield(beliefs, 'ReportFullWorm') && beliefs.ReportFullWorm %KL edit
        [pointR, pointC] = ind2sub(size(data), wormSkel(1));
        warning('off', 'images:bwselect:outOfRange')
        beliefs.FullWorm = bwselect(data, pointC, pointR);
        warning('on', 'images:bwselect:outOfRange')
    end
    if isfield(beliefs, 'ReportSkelInd') && beliefs.ReportSkelInd
        beliefs.SkelInd = wormSkel;
    end
    if isfield(beliefs, 'ReportPosition')
        beliefs.Position = wormSkel;
        if positionError == 0
            beliefs.PositionError = 0;
        else
            beliefs.PositionError = 1;
        end
    end
    if isfield(beliefs, 'ReportSize') && beliefs.ReportSize
        beliefs.SkelSize = sum(sum(beliefs.FullWorm));
    end
end
end

function funcList = parseOptions(outputOptions, expectedVals)
assert((iscell(outputOptions) && isstruct(expectedVals)))
chosenFcns = unique([outputOptions; fieldnames(expectedVals)]);
funcList = cell2struct(cell(numel(chosenFcns),1), chosenFcns);
% functions take in a cell and a colored image (except Position)
% NOTE: whenever the functions on this list are updated, the variable
%    'supportedExpectations' in getSeasonStats must be updated

weightFunction = @(numbers) .5.^(1-mod(double(numbers),2));

if isfield(funcList, 'SkelSize') || isfield(funcList, 'WidthA')
    funcList.SkelSize = @(cell, image)...
        floor(sum(double(image([cell(:)])).*weightFunction(image([cell(:)]))));
end
if isfield(funcList, 'Length') || isfield(funcList, 'WidthA')
    funcList.Length = @(cell, image)...
        floor(sum(weightFunction(image([cell(:)]))));
end
if isfield(funcList, 'Width')
    funcList.Width = @(cell, image) mean((double(image([cell(:)])).*...
        weightFunction(image([cell(:)]))));
end
if isfield(funcList, 'WidthGradient')
    funcList.WidthGradient = getWidthGrad;
end

if isfield(funcList, 'WidthA')
    funcList.WidthA = @(cell, image) funcList.SkelSize(cell, image) ./...
        funcList.Length(cell, image);
end

if isfield(funcList, 'WidthB') 
    widthVals = getWidthVals();
end
if isfield(funcList, 'WidthB')
    funcList.WidthB = @(cell, image) sum(widthVals(image([cell(:)])));
end
end

function widthGrad = getWidthGrad(cell, image)
        locs = false(size(image));
        locs([cell(:)]) = true;
        widthGrad = @(cell, image) getGradientield(image, locs);
end

function [indices, values] = getExpect(funcList, expectList, coloredImage, objList)
% Takes in functions in funcList and evaluates them over objects in objList, 
%    returning each object's 'score'
% The output of this function can likely be stored and reused to improve
%    performance when a worm's "fit with expectations" is needed, if
%    the expectations stay the same
% funcList and expectList are structs with matching fieldnames
%    (funcList can have additional fieldnames)
%    funcList's fields are anonymous functions, expectList's are 
%    numeric (1x1 for 'Length', 'Width', 'SkelSize'; 1xN for 'Position')
% coloredImage is an RxC int16 colored worm skeleton

if ~isa(coloredImage, 'double')
    coloredImage = double(coloredImage);
end
if isempty(objList)
    indices = [];
    values = [];
    return
end

criteria = intersect(fieldnames(expectList), fieldnames(funcList));
minFuncAcc = {};
% Assess the sizes of worms using different criteria
% Criteria must have either a 'CriteriaConfidence' or 'CriteriaStd'
%    field to be considered- 'Confidence' means to normalize by the 
%    expected value (it's probably not a great feature)
for i = 1:numel(criteria)
    criterion = criteria{i};
    assert(isa(expectList.(criterion), 'double'))
    diffFun = @(cell) 0;
    if isfield(expectList, [criterion 'Weight']) && ~isempty(expectList.Weight)
        wt = expectList.Weight;
    else
        wt = 1;
    end
    if isfield(expectList, [criterion 'Confidence']) && ...
            isfield(expectList, [criterion 'Std']) && ...
            ~isempty(expectList.([criterion 'Confidence'])) && ...
            ~isempty(expectList.([criterion 'Std']))
        error(['Specifying an expectation-normalized confidence level ' ...
            'is incompatible with specifying a standard deviation'])
    elseif isfield(expectList, [criterion 'Confidence']) &&...
            ~isempty(expectList.([criterion 'Confidence']))
        diffFun = @(cell) expectList.([criterion 'Confidence']) .*...
            (abs(funcList.(criterion)(cell,...
            coloredImage) - expectList.(criterion)) ./ expectList.(criterion));
    elseif isfield(expectList, [criterion 'Std']) && ...
            ~isempty(expectList.([criterion 'Std']))
        diffFun = @(cell) abs(funcList.(criterion)(cell,...
            coloredImage) - expectList.(criterion)) ./ ...
            (expectList.([criterion 'Std']));
%     elseif isequal(criterion, 'Position')
%         diffFun = @(cell) funcList.Position(cell,expectList.Position);
%         
    elseif ~contains(criterion, 'Confidence') && ...
            ~contains(criterion, 'Std') && ~contains(criterion, 'Weight')
        warning(['Option ' criterion ' could not be used, because '...
            'there was no way to assess proximity'])
    else
        diffFun = @(cell) 0;
    end
    minFuncAcc{end+1} = @(cell) wt .* diffFun(cell);
end

% For each cell in objCells, evaluates each of funcCells and sums them
% May be possible to save time processing by optionally reporting fcn 
%    values here instead of in relativeFrameEvaluation
minFunc = @(objCells, funcCells) cellfun(@(obj) sum(cellfun(@(func)...
    func(obj), funcCells)), objCells);

% Smaller (better) values are earlier
[values, indices] = sort(minFunc(objList, minFuncAcc));
end