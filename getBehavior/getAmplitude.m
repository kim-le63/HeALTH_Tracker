function Amplitude = getAmplitude(xVals,yVals,frxn)

%frxn is the outer end (usually .1) that will be ignored in getting the
%amplitude
% lBound = ceil(numel(xVals).*frxn);
% uBound = floor(numel(xVals).*(1-frxn));

xlist = find(xVals);
ylist = find(yVals);

lBound = ceil(numel(xlist).*frxn);
uBound = floor(numel(ylist).*(1-frxn));

x1 = xVals(lBound);
x2 = xVals(uBound);
y1 = yVals(lBound);
y2 = yVals(uBound);

Amplitude = 0;
for i = lBound+1:1:uBound-1
    x0 = xVals(i);
    y0 = yVals(i);
    Distance = abs( (y2-y1).*x0-(x2-x1)*y0+x2*y1-y2*x1 ) ./ ...
        sqrt( (y2-y1)^2+(x2-x1)^2 );
    if Distance > Amplitude
        Amplitude = Distance;
    end
end
end