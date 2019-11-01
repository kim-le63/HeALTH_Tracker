function vid = importVideo(fileName)
try
    fileName = char(fileName(fileName ~= '"'));
    vRO = VideoReader(fileName);
    vid = zeros(1024,1280,64);
    vid = uint8(vid);
    count = 1;
    while vRO.hasFrame()
        frame = vRO.readFrame();
        vid(:,:,count) = frame(:,:,1);
        count = count+1;
    end
    vid(:,:,count+1:end) = [];

catch ME
    if strcmp(ME.identifier,'MATLAB:audiovideo:VideoReader:CodecNotFound')
        error(['Missing codec for ' fileName])
    elseif strcmp(ME.identifier, 'MATLAB:audiovideo:VideoReader:FileNotFound')
        error(['File ' fileName ' not found'])
    else
        rethrow(ME);
    end
end
end