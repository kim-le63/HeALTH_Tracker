function bodyLineVelocity = getBodyLineVelocity(line1,line2)
%input lines in column format: [Xvals YVals]
%output a column of difference vectors
%General strategy: if one line is 5 points long and the other is 8, 
%   use points 2 thru 6 of the second line
%Establish a correspondence between the central members of the lines, 
%   (using the middle 80% of worm) and calculate the difference between
%   corresponding members
a1 = line1(:,1);
b1 = line1(:,2);
a1 = a1(a1 ~= 0);
b1 = b1(b1 ~= 0);
a2 = line2(:,1);
b2 = line2(:,2);
a2 = a2(a2 ~= 0);
b2 = b2(b2 ~= 0);
A1 = a1(1+floor(numel(a1)*.1):end-floor(numel(a1)*.1));
B1 = b1(1+floor(numel(a1)*.1):end-floor(numel(a1)*.1));
A2 = a2(1+floor(numel(a2)*.1):end-floor(numel(a2)*.1));
B2 = b2(1+floor(numel(a2)*.1):end-floor(numel(a2)*.1));
if ~mod(numel(A1),2) %is even
    A1(end) = [];
end
if ~mod(numel(B1),2)
    B1(end) = [];
end
if ~mod(numel(A2),2)
    A2(end) = [];
end
if ~mod(numel(B2),2)
    B2(end) = [];
end
Size = min(numel(A1),numel(A2));
Buffer1 = floor((numel(A1)-Size)/2);
Buffer2 = floor((numel(A2)-Size)/2);
A1 = A1(1+Buffer1:end-Buffer1);
B1 = B1(1+Buffer1:end-Buffer1);
A2 = A2(1+Buffer2:end-Buffer2);
B2 = B2(1+Buffer2:end-Buffer2);
Line1 = cat(2,A1,B1);
Line2 = cat(2,A2,B2);
try
velocities = Line2 - Line1;
catch
end
bodyLineVelocity = mean(velocities,1);
bodyLineVelocity = bodyLineVelocity';
%below code only calculates speed
%{
half1 = numel(line1)/4;
half2 = numel(line2)/4;
bodyLineVelocity = 0;
for i = 1:min(floor(half1),floor(half2))
    i = i - .01;
    point1L = ceil(half1 - i);
    point1U = floor(half1 + i);
    point2L = ceil(half2 - i);
    point2U = floor(half2 + i);
    dist = sqrt((line1(point1L,1)-line2(point2L,1))^2+...
        (line1(point1U,2)-line2(point2U,2))^2);
    bodyLineVelocity = bodyLineVelocity + dist;
end
bodyLineVelocity = bodyLineVelocity / i;
%}
end