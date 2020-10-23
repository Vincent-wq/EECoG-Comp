function [catA, cntA, listA, cntB, listB] = cntByCat(A,B,TITLE)
nA = length(A);
if nA==length(B)
    N = nA;
else
    catA =[]; cntA = [];
    disp('Length should be the same!')
    return
end
catA = unique(A); nCatA = length(catA);
listA = cell(nCatA,1); cntA = zeros(nCatA,1);
listB = cell(nCatA,1); cntB = zeros(nCatA,1);

for i =1:N
    for j = 1: nCatA
        if A(i) == catA(j)
            listA{j}=[listA{j}, B(i)];
            listB{j}=[listB{j}, i];
        end
    end
end

for i = 1:nCatA
    cntA(i) = length(unique(listA{i}));
    cntB(i) = length(unique(listB{i}));
end
if ischar(TITLE)
    figure
    bar(cntA)
    hold on
    plot(ones(length(cntA))*mean(cntA),'LineStyle', ':', 'Color','r','LineWidth',1.5)
    xticks(1:length(catA))
    if iscell(catA)
        xticklabels(catA)
    else
        xticklabels({num2str(catA)})
    end
    title([TITLE ' Range: ' num2str(min(cntA)) '-' num2str(max(cntA)) ' with mean ' num2str(mean(cntA))])
end
end