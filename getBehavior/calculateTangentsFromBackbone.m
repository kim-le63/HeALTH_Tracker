function [theta] = calculateTangentsFromBackbone(backbone, nBodyPoints)
% CALCULATETANGENTSFROMBACKBONE Calculate tangent angles etc. from the backbone
%
% Given the worm backbone, calculates the tangent angles, sampled at a given
% number of points (equal arclengths along the body).
% (equal arclengths along the body)
%
% INPUT:
% backbone = Nx2 matrix, with each row the (x,y) coordinates of a point along
%       the backbone. Ordered from head (first row) to tail (last row).
% nBodyPoints = sample this number of points along the backbone.
%
% OUTPUT:
% theta = vector with tangent angles, sampled at "nBodyPoints" equally-spaced
%       points along the backbone. The mean has been subtracted (see
%       "thetaMean"), so this vector has zero mean.


x = backbone(:,1)';
y = backbone(:,2)';
n = size(backbone,1);

%% Calculate worm length.
s = zeros(1,n);
for i = 2:n
    s(i) = s(i-1) + sqrt((x(i)-x(i-1))^2+(y(i)-y(i-1))^2);
end
wormLength = s(end);


%% Resample worm at equal arc-lengths.
% Resample using cubic spline interpolation.
curve_spline = spline(s, [x;y]);
dS = max(s)/(nBodyPoints+1);
scale = 0:dS:max(s);
splinePosition = ppval(curve_spline, scale);
splineX = splinePosition(1,:);
splineY = splinePosition(2,:);

% M = diag(3:-1:1, 1);
% 
% %first derivative
% d1 = curve_spline;
% d1.coefs = d1.coefs*M;
% 
% %second derivative
% d2 = d1;
% d2.coefs = d2.coefs*M;
% 
% %find index of maximum curvature
% K = ppval(d2, scale);
% [~, maxKI] = max(K);


%% Calculate tangent angle at each interpolated body point.
theta = zeros(1,nBodyPoints);
for i = 2:length(scale)-1
    % Use a symmetric derivative.
    theta(i-1) = atan2((splineY(i+1)-splineY(i-1))/2, (splineX(i+1)-splineX(i-1))/2);
end
theta = unwrap(theta);
thetaMean = mean(theta);
theta = theta - thetaMean;


% %% Calculate worm thickness at each interpolated body point.
% wormThickness = zeros(1,nBodyPoints);
% for i = 2:length(scale)-1
%     b_dist = zeros(1,size(boundaryXY,1));
%     for j = 1:size(boundaryXY,1)
%         b_dist(j) = sqrt( (splineX(i)-boundaryXY(j,1))^2 + (splineY(i)-boundaryXY(j,2))^2 );
%     end
%     wormThickness(i-1) = min(b_dist);
% end

end
