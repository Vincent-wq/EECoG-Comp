function [Xx,Yx,AUCx] = calPartialROC(X,Y,p)
[dim1, dim2] = size(X); 
Xx = cell(dim1, dim2); Yx = cell(dim1, dim2); AUCx = cell(dim1, dim2);
for idim1 = 1:dim1
    for idim2 = 1:dim2
        tmpX = X{idim1,idim2}; 
        [~,N]= size(tmpX); 
        for iN = 1:N
             tmpx = tmpX{iN}; Diff = abs(tmpx-p);
             [~, indMinDiff] = min(Diff);
             Xx{idim1,idim2}{iN} = X{idim1,idim2}{iN}(1:indMinDiff,1);
             Yx{idim1,idim2}{iN} = Y{idim1,idim2}{iN}(1:indMinDiff,1);
             a = p*p/2;
             AUCx{idim1,idim2}{iN} = (1+(trapz(Xx{idim1,idim2}{iN},Yx{idim1,idim2}{iN})-a)/(p-a))/2;
        end
    end
end

end