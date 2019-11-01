%% Combining and comparing longevity/behavioral information across trials 

% change expCondition depending on what experimental condition to examine/compare
expConditions = {'QL390'}; %example daf-16 conditions. If want to include multiple conditions together, list them within the array. To exclude a specific phrase, use '~'

disp('Select folder with files of interest')
folderName = uigetdir(); % select folder with videos of interest

[combinedBP] = plottingBehavioralParameters(expConditions, folderName); %getting compiled behavioral information across all trials

%% examining avg populational decline and intra-populational differences (short v. long-lived cohorts)
xVals = combinedBP{5,1};
[meanVals, adultPC_rel_indiv,relMeanVals, meanVals_short, relMeanVals_short, meanVals_long, ...
    relMeanVals_long, p_ks, AUC] = combinedRelativeBehavioralDecline(combinedBP{1,1},xVals(:,:,10));

%% PCA
% compiling across all conditions combined individual relative decline
% (in order): N2, daf-16, daf-2, OD10, OD2.5, 20C, 17.5C, 15C, 15C <-> 20C
combAdultPC_rel_indiv = cat(1,adultPC_rel_indiv_QLN2_OD5,adultPC_rel_indiv_QL12,adultPC_rel_indiv_QL390,...
    adultPC_rel_indiv_OD10,adultPC_rel_indiv_OD2_5,adultPC_rel_indiv_20C, adultPC_rel_indiv_17_5C,...
    adultPC_rel_indiv_15C,adultPC_rel_indiv_15C_20C);

[wcoeff,score,latent,tsquared,explained] = pca(combAdultPC_rel_indiv,'centered',false);

%% LASSO regression
[B_LASSO_model,fitInfo_LASSO_model,varToKeep_LASSO_model,residuals_train_model,residuals_test_model, ...
    yhat_test,residuals_test] = regressionModel(combinedBP_QLN2_OD5{2,1},combinedBP_QLN2_OD5{3,1},combinedBP_QLN2_OD5{4,1});
