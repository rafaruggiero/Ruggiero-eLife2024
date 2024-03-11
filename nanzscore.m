% NANZSCORE - Calculate the z-score per dimension, omitting NaN values.
% 
% Authors: Rafael Naime Ruggiero, 2021 

function [zscored] = nanzscore(x,dim)

zscored=zeros(size(x));
if dim == 1
    for irow = 1:size(x,1)
        zscored(irow,:) = (x(irow,:) - mean(x(irow,:),'omitnan'))./std(x(irow,:),'omitnan');
    end
elseif dim == 2
    for icol = 1:size(x,2)
        zscored(:,icol) = (x(:,icol) - mean(x(:,icol),'omitnan'))./std(x(:,icol),'omitnan');
    end
end
end