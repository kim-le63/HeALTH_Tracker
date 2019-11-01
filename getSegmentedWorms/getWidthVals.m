function widthVal = getWidthVals(idx)
% Note that index is 1 to 9: index with widthVals(number + 1)
if nargin < 1
    idx = 9;
end
for dbl = 0:idx-1
    widthVal(dbl+1) = log10(dbl+2).*(1-exp(-(dbl-.75).^2)).^2;
end