function [maskedImg] = circleMe(image)
 % adding a circular ROI
chamberMat = false(size(image,1), size(image,2));
[r2,c2]=meshgrid(1:size(image,2),1:size(image,1)); 
sMask =(((r2-round(size(image,2)/2)).^2+(c2-round(size(image,1)/2)).^2)<=(35)^2); 
chamberMat = chamberMat | sMask;
maskedImg = image.*uint8(chamberMat);         
end