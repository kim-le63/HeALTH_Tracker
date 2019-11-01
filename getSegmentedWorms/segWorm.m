 function [coloredSkel, skelOrient] = segWorm(image, opts, threshold)
% Inputs:
%   image : MxN binary output from some type of thresholding
%   opts: cell containing a double array of a width or range of widths
%   threshold: if image isn't logical, it can be thresholded with this
% Output: MxN uint8 skeleton that is "colored" by the diameter of the worm
%   and the orientations of each diameter found

if nargin < 3
    assert(islogical(image))
else
    assert(~islogical(image))
    image = int16(image > threshold);
end

if nargin < 2
    opts = {'all'};
end

coloredSkel = zeros(size(image),'int16');
skelOrient = zeros([size(image)],'int16');

if any(cellfun(@(cell) strcmp(cell,'all'), opts))
    vec = [1 7];
else
    assert(sum(cellfun(@(cell) isnumeric(cell), opts)) == 1, ...
        'only one numeric option, width/width range, is allowed')
    numIdx = find(cellfun(@(cell) isnumeric(cell), opts));
    assert(numel(opts{numIdx}) > 0 && numel(opts{numIdx}) <= 2,...
        'width must be a single number or range')
    vec = opts{numIdx};
end

if numel(vec) == 1
    vec = [vec vec];
end

if mod(vec(1),2) == 0
    [skel, orient] = findWormHelper(image, vec(1)-1, false, true);
    skelOrient = skelOrient + cat(3, -1+2.*mod(orient,2), ...
        -1+2.*(orient <= 2)) .* int16(coloredSkel == 0 & skel ~= 0); 
    coloredSkel = coloredSkel + skel;
    vec(1) = vec(1)+1;
end

widths = vec(1):2:vec(2);
for i = widths
    if i ~= vec(2)
        [skel, orient] = findWormHelper(image, i, true, true);
    else
        [skel, orient] = findWormHelper(image, i, true, false);
    end
    skelOrient = skelOrient + int16(coloredSkel == 0 & skel ~= 0) .* ...
        cat(3, -1+2.*mod(orient, 2), -1+2.*(orient <= 2));
    coloredSkel = coloredSkel + int16(coloredSkel == 0) .* skel;
end

end

function [skelOut, orient] = findWormHelper(imgIn, width, inclOdd, inclEven)
% Takes in binary image imgIn containing a single object
% Returns 2*width, where width is a multiple of .5
if nargin < 3
    inclEven = true;
    inclOdd = true;
end
assert(islogical(imgIn))
imgIn = int16(imgIn);
width = int16(width);
wormStrel = int16(4*getWormStrel(width)); 
wormStrelD = int16(4*getWormStrel(width, 'diagonal'));

skelOut(:,:,1) = int16(conv2(imgIn, wormStrelD, 'same'));
skelOut(:,:,2) = int16(conv2(imgIn, wormStrel, 'same'));
skelOut(:,:,3) = int16(conv2(imgIn, wormStrel', 'same'));
skelOut(:,:,4) = int16(conv2(imgIn, wormStrelD(:, end:-1:1), 'same'));
%skelOut(skelOut < 4.*width-1) = intmax('int16');
[skelOut, orient] = max(skelOut,[],3,'omitnan');
orient(isnan(skelOut) | skelOut == intmax('int16')) = 0;
skelOut(isnan(skelOut) | skelOut == intmax('int16')) = 0;

orient = int16(inclOdd) .* int16(orient .* (skelOut == (4.*width))) + ...
    int16(inclEven) .* int16(orient .* (skelOut == (4.*width-1)));
skelOut = int16(inclOdd) .* width .* int16(skelOut == (4.*width)) + ...
    int16(inclEven) .* (width+1) .* int16(skelOut == (4.*width-1));
end
