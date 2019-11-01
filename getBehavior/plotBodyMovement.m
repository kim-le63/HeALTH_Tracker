%% plotting function  - plots both individuals separately and within a combined plot

function plotBodyMovement(data, fieldOfInterest, xVal, yVal, titleVal, chamberIDs, day)
if day == 0 % plotting by individual videos
    % setting up cell to workable matrix array
    if isfield(data{1,1},fieldOfInterest) == 0
        disp('field does not exist')
        return
    end

    tempData = zeros(size(data,1),size(data,2)-1);
    for i = 1:size(data,2)-1
        for k = 1:size(data,1)
            tempData(k,i) = data{k,i}.(fieldOfInterest);
        end
    end

    % plotting individuals together
    tempMeanVal = nanmean(tempData);
    Markers = {'*','o','x'};

    figure; hold on; title(['Averaged ',titleVal])
    iter = 1;
    t = 1:1:size(tempData,2);
    for k = 1:length(chamberIDs)
        scatter(t,tempData(k,:),'Marker',Markers{iter})
        iter = iter+1;
        if iter==4
            iter = 1;
        end
    end
    plot(tempMeanVal,'LineWidth',2,'Color','k')
    xlabel(xVal); ylabel(yVal)

    % plotting individuals separately 
    figure; hold on; 
    t = 1:1:size(tempData,2);
    for k = 1:length(chamberIDs)
        subplot(5,ceil(length(chamberIDs)/5),k); 
        hold on; title([titleVal, ': Worm ',num2str(chamberIDs(k))])
        scatter(t,tempData(k,:),'*'); xlabel(xVal); ylabel(yVal)
    end
else %plotting by days
    % setting up cell to workable matrix array
    if isfield(data,fieldOfInterest) == 0
        disp('field does not exist')
        return
    end

    tempData = data.(fieldOfInterest);

    % plotting individuals together
    tempMeanVal = nanmean(tempData);
    Markers = {'*','o','x'};

    figure; hold on; title(['Averaged ',titleVal])
    iter = 1;
    t = 1:1:length(tempData);
    for k = 1:length(chamberIDs)
        scatter(t,tempData(k),'Marker',Markers{iter})
        iter = iter+1;
        if iter==4
            iter = 1;
        end
    end
    plot(tempMeanVal,'LineWidth',2,'Color','k')
    xlabel(xVal); ylabel(yVal)

    % plotting individuals separately 
    figure; hold on; 
    t = 1:1:length(tempData);
    for k = 1:length(chamberIDs)
        movAvg = movmedian(tempData(k,:),20,'omitnan');
        subplot(5,ceil(length(chamberIDs)/5),k); 
        hold on; title([titleVal, ': Worm ',num2str(chamberIDs(k))])
        scatter(t,tempData(k,:),'*'); xlabel(xVal); ylabel(yVal)
        plot(movAvg,'LineWidth',2,'Color','k')
    end    
end

end