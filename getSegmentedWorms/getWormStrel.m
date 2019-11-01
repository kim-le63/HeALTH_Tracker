function wormStrel = getWormStrel(wormWidth,opt)
% convolution output:
%   wormWidth where there is a worm,
%   wormWidth-.25 where the worm is 1 pixel thicker (skeleton is 2 px thicc)
%   wormWidth-.5 where the worm is 2 pixels thicker (skeleton is 1 px thick)
if nargin < 2 || ~exist('opt','var') || isempty(opt)
    opt = 'straight';
end
assert(any(cellfun(@(cell) strcmp(cell,opt),{'straight', 'diagonal'})))
assert(mod(wormWidth,2) == 1, 'filter must be centered on one pixel')

wormStrel = ones(1,wormWidth+4);
wormStrel([1,end]) = -2;
wormStrel([2,end-1]) = -.25;
if strcmp(opt,'diagonal')
    wormStrel = diag(wormStrel);
end
end