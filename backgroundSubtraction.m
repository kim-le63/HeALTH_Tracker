function [background] = backgroundSubtraction(video)
background = min(video(:,:,1:15:end), [], 3); % median if moving objects weren't brighter
end