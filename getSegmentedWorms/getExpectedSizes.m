function expectations = getExpectedSizes(fileName)
% adjust as needed based on video recording scheme (i.e. what videos occur
% when the worm is L4, Day 1 adult, etc.)

numName=sscanf(fileName,'AutoVid%d');
if numName <= 10
    expectations = struct('SkelSize', 140, 'SkelSizeStd', 25);
elseif numName <= 30
	expectations = struct('SkelSize', 260, 'SkelSizeStd', 25);
elseif numName <= 50
	expectations = struct('SkelSize', 335, 'SkelSizeStd', 40);
elseif numName <= 100
	expectations = struct('SkelSize', 400, 'SkelSizeStd',45);
elseif numName <= 300
	expectations = struct('SkelSize', 450, 'SkelSizeStd', 50);
else
    expectations = struct('SkelSize', 500, 'SkelSizeStd', 55);
end

end
