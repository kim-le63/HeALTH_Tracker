%% freq fixer
% fix max freq value
folderName = 'C:\Users\kim63\Dropbox (GaTech)\MATLAB\healthspanPlatform_Tracker\healthspanBehavioralAnalysis';

fileList_movementData = dir(fullfile(folderName,'movementData_*.*')); fileList_movementData = {fileList_movementData.name};
fileList_lifespanData = dir(fullfile(folderName,'lifespanData_*.*')); fileList_lifespanData = {fileList_lifespanData.name};

for j = 1:length(fileList_movementData)
    fullPath_movement = [folderName '\' fileList_movementData{j}]; tempMove = load(fullPath_movement);
    fullPath_lifespan = [folderName '\' fileList_lifespanData{j}]; tempLifespan = load(fullPath_lifespan);

    freqData = tempMove.([fileList_movementData{j}(1:end-4)]).freqData;
    for i = 1:size(freqData,1) %for each chamber of interest
        for ii = 1:size(freqData,2)-1 %for each vid
            if ~isnan(freqData{i,ii}.meanVal) %find actual max freq value
                freqData{i,ii}.maxVal = max(freqData{i,ii}.freqVals);
            end
        end
    end
    
    deathVid = tempLifespan.([fileList_lifespanData{j}(1:end-4)]).deathVid;
    
    [declinePT, relativeHighPeriod, declineSlope_all, declineInt_all, declineSlope_high, declineInt_high, maxValDay, avgValDay, ...
    maxValDay_high, avgValDay_high, threshVal, declinePT_75, relativeHighPeriod_75, declinePT_duration, ...
    relativeHighPeriod_duration] = getBehaviorParameters(freqData, 'maxVal', 'none', 0, 'freq(hz)',deathVid);

    for k =1:size(freqData,1)
        freqData{k,size(freqData,2)}.maxVal.declinePT = declinePT(k); freqData{k,size(freqData,2)}.maxVal.relativeHighPeriod = relativeHighPeriod(k);
        freqData{k,size(freqData,2)}.maxVal.declineSlope_all = declineSlope_all(k); freqData{k,size(freqData,2)}.maxVal.declineInt_all = declineInt_all(k);
        freqData{k,size(freqData,2)}.maxVal.declineSlope_high = declineSlope_high(k); freqData{k,size(freqData,2)}.maxVal.declineInt_high = declineInt_high(k);
        freqData{k,size(freqData,2)}.maxVal.maxValDay = maxValDay(k,:); freqData{k,size(freqData,2)}.maxVal.avgValDay = avgValDay(k,:);
        freqData{k,size(freqData,2)}.maxVal.maxValDay_high = maxValDay_high(k,:); freqData{k,size(freqData,2)}.maxVal.avgValDay_high = avgValDay_high(k,:);
        freqData{k,size(freqData,2)}.maxVal.threshVal = threshVal;
        freqData{k,size(freqData,2)}.maxVal.declinePT_75 = declinePT_75(k); freqData{k,size(freqData,2)}.maxVal.relativeHighPeriod_75 = relativeHighPeriod_75(k);
        freqData{k,size(freqData,2)}.maxVal.declinePT_duration = declinePT_duration(k); freqData{k,size(freqData,2)}.maxVal.relativeHighPeriod_duration = relativeHighPeriod_duration(k);
    end

    ampData = tempMove.([fileList_movementData{j}(1:end-4)]).ampData;
    csData = tempMove.([fileList_movementData{j}(1:end-4)]).csData;
    dpData = tempMove.([fileList_movementData{j}(1:end-4)]).dpData;
    areaData = tempMove.([fileList_movementData{j}(1:end-4)]).areaData;
    PCData = tempMove.([fileList_movementData{j}(1:end-4)]).PCData;
    
    movementData.([fileList_movementData{j}(1:end-4)]).freqData = freqData;
    movementData.([fileList_movementData{j}(1:end-4)]).ampData = ampData;
    movementData.([fileList_movementData{j}(1:end-4)]).csData = csData;
    movementData.([fileList_movementData{j}(1:end-4)]).dpData = dpData;
    movementData.([fileList_movementData{j}(1:end-4)]).areaData = areaData;
    movementData.([fileList_movementData{j}(1:end-4)]).PCData = PCData;
    
    save(fileList_movementData{j},'-struct','movementData');
end
% 
% 
% save('movementData_DRTrial_Exp45_OD10.mat','movementData_DRTrial_Exp45_OD10');
% save('movementData_DRTrial_Exp51_OD10.mat','movementData_DRTrial_Exp51_OD10');
% save('movementData_DRTrial_ExpA_OD10.mat','movementData_DRTrial_ExpA_OD10');
% save('movementData_DRTrial_ExpB_OD10.mat','movementData_DRTrial_ExpB_OD10');
% 
% save('movementData_DRTrial_Exp47_OD2_5.mat','movementData_DRTrial_Exp47_OD2_5');
% save('movementData_DRTrial_Exp48_OD2_5.mat','movementData_DRTrial_Exp48_OD2_5');
% save('movementData_DRTrial_Exp53_OD2_5.mat','movementData_DRTrial_Exp53_OD2_5');
% save('movementData_DRTrial_Exp54_OD2_5.mat','movementData_DRTrial_Exp54_OD2_5');
% save('movementData_DRTrial_ExpE_OD2_5.mat','movementData_DRTrial_ExpE_OD2_5');
% save('movementData_DRTrial_ExpF_OD2_5.mat','movementData_DRTrial_ExpF_OD2_5');
% 
% save('movementData_DRTrial_ExpC_OD5_v2.mat','movementData_DRTrial_ExpC_OD5_v2');
% save('movementData_DRTrial_ExpD_OD5_v2.mat','movementData_DRTrial_ExpD_OD5_v2');
% save('movementData_DRTrial_Exp49_OD5_v2.mat','movementData_DRTrial_Exp49_OD5_v2');
% save('movementData_DRTrial_Exp50_OD5_v2.mat','movementData_DRTrial_Exp50_OD5_v2');
% 
% save('movementData_mutantTrial_Exp20_QLN2_v2.mat','movementData_mutantTrial_Exp20_QLN2_v2');
% save('movementData_mutantTrial_Exp21_QLN2_v2.mat','movementData_mutantTrial_Exp21_QLN2_v2');
% save('movementData_mutantTrial_Exp29_QLN2_v2.mat','movementData_mutantTrial_Exp29_QLN2_v2');
% save('movementData_mutantTrial_Exp60_QLN2_v2.mat','movementData_mutantTrial_Exp60_QLN2_v2');
% save('movementData_mutantTrial_Exp65_QLN2_v2.mat','movementData_mutantTrial_Exp65_QLN2_v2');
% 
% save('movementData_mutantTrial_Exp32_QL12.mat','movementData_mutantTrial_Exp32_QL12');
% save('movementData_mutantTrial_Exp61_QL12_v2.mat','movementData_mutantTrial_Exp61_QL12_v2');
% save('movementData_mutantTrial_Exp62_QL12.mat','movementData_mutantTrial_Exp62_QL12');
% save('movementData_mutantTrial_Exp66_QL12.mat','movementData_mutantTrial_Exp66_QL12');
% 
% save('movementData_mutantTrial_Exp18_QL390.mat','movementData_mutantTrial_Exp18_QL390');
% save('movementData_mutantTrial_Exp31_QL390.mat','movementData_mutantTrial_Exp31_QL390');
% 
% save('movementData_tempTrial_Exp0_20C.mat','movementData_tempTrial_Exp0_20C');
% save('movementData_tempTrial_Exp15_15C_20C.mat','movementData_tempTrial_Exp15_15C_20C');
% save('movementData_tempTrial_Exp6_20C.mat','movementData_tempTrial_Exp6_20C');
% save('movementData_tempTrial_Exp9_15C.mat','movementData_tempTrial_Exp9_15C');
% 
