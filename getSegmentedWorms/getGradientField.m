function gradientField = getGradientField(binaryImage, obsField, wIdx)
% Given a MxNx2 set of cardinal + ordinal orientations, returns gradient
% (first input is the second output of segWorm)
% Note that distance must be conserved e.g. the distance between N-S and
%    E-W must be larger than the distance between N-S and NE-SW in orientationField input
%    e.g. (note that example also has a set, low distance to 0): 
%        N-S: [-1;; 1]
%        NE-SW: [1;; 1]
%        E-W: [1;; -1]
%        NW-SE: [-1;; -1]
%        (;; as these are pairs of 2D gradient values)
% Also note that it may be useful to raise the floor of the gradient by
%    taking the maximum of (gradientField-2, 0)
obsField = double(obsField);
obsField(:,:,sum(sum(obsField,1),2)==0) = [];
if ~islogical(binaryImage)
    binaryImage = binaryImage ~= 0;
end
gradientField = zeros(size(obsField));
for i = 1:size(obsField, 3)
    gradientField(:,:,i) = binaryImage.*imgradient(obsField(:,:,i), 'intermediate');
    if i == wIdx
        gradientField(:,:,i) = max(gradientField(:,:,i)...
            - sqrt(2).*obsField(:,:,i), 0);
    else
        gradientField(:,:,i) = max(gradientField(:,:,i) - sqrt(2), 0);
    end
end
end