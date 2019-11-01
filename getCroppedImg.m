function [cImg] = getCroppedImg(center, rad, video)
  
if center(1)-rad < 1 || center(1)+rad > size(video,1) || ...
        center(2)-rad < 1 || center(2) + rad > size(video,2)
    %go to padding function
    padArray = padMe(center, rad, video, 1);
    cImg = padArray;
else
    cImg = video(center(1)-rad:center(1)+rad,...
    center(2)-rad:center(2)+rad,1); % single-frame image, NOT VIDEO
end

end

function padArray = padMe(center, rad, video, image)
padAmounts=zeros(4,1);

if center(1)-rad < 1
    padAmounts(1) = abs(center(1)-rad)+1;
elseif center(1)+rad > size(video,1)
    padAmounts(2) = (center(1)+rad)-size(video,1);
end

if center(2)-rad <1
    padAmounts(3) = abs(center(2)-rad)+1;
elseif center(2)+rad > size(video,2)
    padAmounts(4) = (center(2)+rad)-size(video,2);
end

if image == 1
    cImg = video(center(1)-rad+padAmounts(1):center(1)+rad-padAmounts(2),...
            center(2)-rad+padAmounts(3):center(2)+rad-padAmounts(4)); % single-frame image, NOT VIDEO
    %padding array
    padArray = padarray(cImg,[padAmounts(1) padAmounts(3)],'pre');
    padArray = padarray(padArray,[padAmounts(2) padAmounts(4)],'post');
else
    cVid = video(center(1)-rad+padAmounts(1):center(1)+rad-padAmounts(2),...
            center(2)-rad+padAmounts(3):center(2)+rad-padAmounts(4),:); % VIDEO
    %padding array
    padArray = padarray(cVid,[padAmounts(1) padAmounts(3)],'pre');
    padArray = padarray(padArray,[padAmounts(2) padAmounts(4)],'post');
end

end