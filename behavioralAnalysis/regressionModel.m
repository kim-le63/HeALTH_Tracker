function [B_LASSO_model,fitInfo_LASSO_model,varToKeep_LASSO_model,residuals_train_model, ...
    residuals_test_model,yhat_test,residuals_test] = regressionModel(y_day_censor, X_censor, varNames, X_pop, Y_pop)
% getting LASSO regression model from population
% INPUTS:
% y_day_censor - time of death for each individual excluding censored
% individuals
% X_censor - behavioral information for each individual over time
% varNames - names of behavioral metric being examined
% X_pop - population of individuals to be fitted onto the model
% Y_pop - time of death of individuals to be fitted onto the model
% OUTPUTS:
% B_LASSO_model - fitted coefficient values for LASSO model
% fitInfo_LASSO_model - fit information for LASSO model
% varToKeep_LASSO - names of behavioral metrics selected for the model
% residuals_train_model - residuals for training data set
% residuals_test_model - residuals for testing data set
% yhat_test - predicted lifespans for individuals inputted into the model
% residuals_test - residuals for individuals inputted into the model

% normalization of dataset
[cX_censor, varNames, imputedX, varIdx] = normNaN_X (X_censor,y_day_censor,varNames); % normalize predictor values and ID columns with many NaNs

% splitting dataset into testing and training data for model
n = length(y_day_censor); c = cvpartition(n,'HoldOut',0.3); %training w/ 70% of data
idxTrain = training(c,1); idxTest = ~idxTrain;
XTrain = cX_censor(idxTrain,:); yTrain= y_day_censor(idxTrain);
XTest = cX_censor(idxTest,:); yTest = y_day_censor(idxTest);

[B_LASSO_model,fitInfo_LASSO_model,varToKeep_LASSO_model,residuals_train_model, ...
    residuals_test_model] = getLASSO(XTrain,yTrain,XTest,yTest,varNames,'');

% testing (different) experimental groups on the model
if nargin > 4
    [yhat_test,residuals_test] =  getFittedResponse(fitInfo_LASSO_model,B_LASSO_model, ...
        varToKeep_LASSO_model,varNames, X_pop, Y_pop);
else
    [yhat_test,residuals_test] =  getFittedResponse(fitInfo_LASSO_model,B_LASSO_model, ...
        varToKeep_LASSO_model,varNames, XTest, yTest);
end

end

function[yhat, residuals] = getFittedResponse(fitInfo,B,varToKeep_LASSO,varNames,X,y)
idxLambda1SE = fitInfo.Index1SE;
coef_LASSO = B(:,idxLambda1SE); 
coef0_LASSO = fitInfo.Intercept(idxLambda1SE);

modX = []; mod_LASSO_coef = [];
coeffVal = find(coef_LASSO);

for i = 1:length(varToKeep_LASSO)
    for j = 1:length(varNames)
        if contains(varNames{j},varToKeep_LASSO{i,1}) 
            modX = cat(2,modX,X(:,j));
            mod_LASSO_coef = cat(1,mod_LASSO_coef,coef_LASSO(coeffVal(i)));
        end
    end
end

yhat = modX*mod_LASSO_coef + coef0_LASSO;
residuals = y - yhat;
end


function [cX, varNames, imputedX, varIdx] = normNaN_X (X,y,iVarNames) % normalize predictor values and ID columns with many NaNs
%% v2 normalization code
cX = []; imputedX = []; varNames =[]; varIdx = [];
iter = 1;
for i = 1:size(X,2)
    if ~contains(iVarNames{i},'threshVal') && ~contains(iVarNames{i},'declinePT') &&  ~contains(iVarNames{i},'relativeHighPeriod')  && ~contains(iVarNames{i},'declineSlope') &&  ~contains(iVarNames{i},'declineInt')
        if sum(isnan(X(:,i))) < round(size(X,1)/2) && sum(X(:,i)) >0 % if variable not too sparse
            [idxNaN] = find(isnan(X(:,i))); 
            for j = 1:length(idxNaN)
               X(idxNaN(j),i) = nanmedian(X(:,i));
            end
            imputedX(:,iter) = X(:,i);
            cX(:,iter) = normalize(X(:,i));
            varNames{iter} = iVarNames{i};
            varIdx{iter} = i;
            iter = iter+1;
        end
    elseif contains(iVarNames{i},'declinePT_duration') || contains(iVarNames{i},'relativeHighPeriod_duration') 
        if sum(isnan(X(:,i))) < round(size(X,1)/2) && sum(X(:,i)) >0 % if variable not too sparse
            [idxNaN] = find(isnan(X(:,i))); 
            for j = 1:length(idxNaN)
               X(idxNaN(j),i) = nanmedian(X(:,i));
            end
            imputedX(:,iter) = X(:,i);
            cX(:,iter) = normalize(X(:,i));
            varNames{iter} = iVarNames{i};
            varIdx{iter} = i;
            iter = iter+1;
        end
    end
end

end

%% LASSO
function [B,fitInfo,varToKeep, residuals_train,residuals_test] = getLASSO(xTrain,yTrain,XTest,yTest,varNames,condition)
[B,fitInfo] = lasso(xTrain,yTrain,'Alpha',1,'CV',10,'Standardize',true);
lassoPlot(B,fitInfo,'PlotType','CV'); title ([condition,' LASSO - standardized'])

idxLambda1SE = fitInfo.Index1SE;
coef = B(:,idxLambda1SE); coef0 = fitInfo.Intercept(idxLambda1SE);

yhat_train = xTrain*coef + coef0;
figure; subplot(2,2,1); hold on; scatter(yTrain,yhat_train,'bo'); plot(yTrain, yTrain)
xlabel('Observed Response'); ylabel('Fitted Response'); title([condition,' LASSO - fitted v. trained for training data'])
residuals_train = yTrain - yhat_train;
subplot(2,2,3); hold on; stem(residuals_train)
xlabel('Train Observation'); ylabel('Residual');

yhat_test = XTest*coef + coef0;
subplot(2,2,2); hold on; scatter(yTest,yhat_test,'bo'); plot(yTest,yTest)
xlabel('Observed Response'); ylabel('Fitted Response'); title([condition,' LASSO - fitted v. trained for test data'])
residuals_test = yTest - yhat_test;
subplot(2,2,4); hold on; stem(residuals_test)
xlabel('Test Observation'); ylabel('R esidual');

[rows,~] = find(B(:,idxLambda1SE)); uRows = unique(rows);
varToKeep = [];
for i = 1:length(uRows)
    varToKeep{i,1} = varNames{uRows(i)};
    varToKeep{i,2} = B(uRows(i),idxLambda1SE);
end
end



