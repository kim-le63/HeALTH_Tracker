function [chamberCenters, chamberMask] = createChamberMask(rgImage,radius)
% based on user input, define chamber locations (60 chamber lifespan device)

% INPUTS:
% rgImage - the mxn image
% radius (optional)- the radius of chambers (default 25)
% OUTPUTS:
% chamberCenters - 60x2 double coordinates of the centers of chambers
% chamberMask - mxn mask of the areas inside chambers

if nargin < 2
    radius = 25;
end
chamFig = figure(1000);
chamberAxis = axes();
vidHeight = size(rgImage,1);
vidWidth = size(rgImage,2);
k = false;
dStatus = false;
while ~k
radii = ones(60,1) .* radius;

axis image
imshow(rgImage)
title('Click on the centers of the corner chambers')
[r, c] = ginput(4); % grab four corner points
hold on; 

locMat = zeros(4,2);
chamberCenters = zeros(60, 2);
chamberMat = false(vidHeight, vidWidth);

%% Determine which corner is which 
[~, i] = sort(r); %sorts in ascending order
[~, j] = sort(c);



for k= 1:2
    locMat(i(k),1) = 1; %1 for 2 lower values
    locMat(j(k),2) = 1; %1 for 2 lower values
end

checkMat = [1 1;0 1; 1 0; 0 0];
for f= 1:4
    b = ismember(locMat, checkMat(f,:), 'rows');
    idx = find(b);
    if sum(b) == 1
        if f ==1
            topleft = [r(idx), c(idx)];
        elseif f==2
            topright = [r(idx), c(idx)];
        elseif f== 3
            bottomleft = [r(idx), c(idx)];
        elseif f== 4
            bottomright = [r(idx), c(idx)];
        end
    else
        dStatus = true;
        disp('manual chamber identification error')
        createChamberMask(chamberAxis, vidHeight, vidWidth, rgImage);
    end
    
end


%% Interpolate centers of all chambers

if ~dStatus

    leftHeightComp = bottomleft-topleft; % x_distance and y_distance between top left and bottom left corners
    rightHeightComp = bottomright-topright; 
    heightComp = [leftHeightComp; rightHeightComp];

    topCorners = [topleft; topright];
    leftAndRightCenters = zeros(5, 4); %first two cols are (x,y) centers of far left column, 3rd and 4th cols are (x,y) centers of far right column

    for h=1:2 %h = 1 corresponds to left col centers, h=2 corresponds to right col centers
        height = sqrt(sum(heightComp(h, :).^2));                                      % distance between top and bottom left corners
        step = height/4;                                                    % step depends how many rows there are. Step = height /(nRow-1)
        angle = atan2(heightComp(h,2), heightComp(h,1));                            % return the angle between the x-axis and left-corners line

        %calculate x and y coordinates for chambers
        leftAndRightCenters(:,1+(h-1)*2)=topCorners(h,1)+(0:4)*step*cos(angle);   %x            
        leftAndRightCenters(:,2+(h-1)*2)=topCorners(h,2)+(0:4)*step*sin(angle);   %y  
    end


    % Calculate the coordinates of each chamber in each row (same process used above, applied to rows)
    for row=1:5
        widthComp=leftAndRightCenters(row,3:4)-leftAndRightCenters(row,1:2);                                        % x_distance and y_distance between far right and far left chambers of the row
        width=sqrt(sum(widthComp.^2));                                             % distance between the far right and far left chambers of the row
        step=width/11;                                                        % step depends how many chambers there are in a row. Step = width /(nColumn-1)
        angle=atan2(widthComp(2), widthComp(1));                              % return the angle between the y-axis and the line of the row
        chamberCenters(12*(row-1)+(1:12),1)=leftAndRightCenters(row,1)+(0:11)*step*cos(angle);  % calculate the x coordinates for all the chambers of the same row
        chamberCenters(12*(row-1)+(1:12),2)=leftAndRightCenters(row,2)+(0:11)*step*sin(angle);  % calculate the y coordinates for all the chambers of the same row
    end
    chamberCenters=round(chamberCenters);

    %% Create mask
    [r,c]=meshgrid(1:vidWidth,1:vidHeight); 
    % create a mask for one chamber (1's where chamber is)
    %and add mask to mask matrix
    for i=1:60
        sMask =(((r-chamberCenters(i,1)).^2+(c-chamberCenters(i,2)).^2)<=(radii(1))^2);           
        chamberMat = chamberMat | sMask;                                                        
    end
end

title('click mouse to redo chamber ID; otherwise press any key')
viscircles(chamberCenters, radii); 

k = waitforbuttonpress;
hold off
end
chamberMask = chamberMat;
close(chamFig)