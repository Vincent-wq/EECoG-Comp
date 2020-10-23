function [ind, cnt] = countCatNum1d(Data,TITLE)
% Count the number of categorial data values, and plot bar diagram.
% by Vincent 2019.4.26
[ind, jj, kk]=unique(Data);
cnt=histc(kk,1:numel(ind));
if ischar(TITLE)
    figure
    bar(cnt)
    hold on
    plot(ones(length(cnt))*mean(cnt),'LineStyle', ':', 'Color','r','LineWidth',1.5)
    xticks(1:length(ind))
    if iscell(ind)
        xticklabels(ind)
    else
        xticklabels({num2str(ind)})
    end
    title([TITLE ' Range: ' num2str(min(cnt)) '-' num2str(max(cnt)) ' with mean ' num2str(mean(cnt))])
end
end