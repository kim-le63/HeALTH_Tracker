function intensity = getExpectedIntensity(vid)
    avgInt = mean2(vid(:,:,30)); 
    if avgInt >= 50 %rough estimate for device OD
        intensity = 100;
    else
        intensity = 90;
    end
end