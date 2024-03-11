% STATEMAP2D() - Projects spectral power estimates on 2D map that separates
%       sleep-wake states
%    
%   Usage
%       [ratio1,ratio2] = statemap2D(data,f)
% 
%   Inputs 
%       data = time-frequency power spectra per region (frequency,epochs,region)
%           Original article used 2-s epochs in 1-s steps
%       f = frequency vector
% 
%   Outputs 
%       ratio1 = PC1 of ratio 1 (0.5-20/0.5-55) of regions;
%       ratio2 = PC1 of ratio 2 (0.5-4.5/0.5-9) of regions;
% 
% References:
%   Gervasoni et al. Global forebrain dynamics predict rat behavioral 
%       states and their transitions. J Neurosci. 2004. 
%       doi: 10.1523/JNEUROSCI.3524-04.2004.
%   Dzirasa et al. Dopaminergic control of sleep-wake states. J Neurosci.
%       2006. doi: 10.1523/JNEUROSCI.1767-06.2006.
% 
% Author: Danilo Benette Marques, 2020

function [ratio1, ratio2] = statemap2D(data,f)
% 2D state map

%Define frequencies of ratios (based on Dzirasa et al., 2006)
f_ratio1n = find(f>=2 & f<=20);
f_ratio1d = find(f>=2 & f<=55);

f_ratio2n = find(f>=2 & f<=4.5);
f_ratio2d = find(f>=2 & f<=9);

%Calculate the ratios per regions
for region = 1:size(data,3)
ratio1(region,:) = sum(data(f_ratio1n,:,region),1) ./ sum(data(f_ratio1d,:,region),1);
ratio2(region,:) = sum(data(f_ratio2n,:,region),1) ./ sum(data(f_ratio2d,:,region),1);
end

%PCA of ratios of all regions
[~,pcscore] = pca(ratio1');
ratio1 = pcscore(:,1);

[~,pcscore] = pca(ratio2');
ratio2 = pcscore(:,1);

%Hann window smoothing
% original articles uses 20 s (adjust according to length)
ratio1 = conv(ratio1,hann(10),'same');
ratio2 = conv(ratio2,hann(10),'same');

plotmap = 0; %1 to plot statemap
if plotmap
figure,
subplot(2,1,1)
    hist(ratio1,100)
    ylabel('Ratio 1')
    title('Ratio 1 distribution')
subplot(2,1,2);
    hist(ratio2,100);
    ylabel('Ratio 2')
    title('Ratio 2 distribution')

figure,scatter(ratio2,ratio1,'.k','sizedata',10)
xlabel('Ratio 2')
ylabel('Ratio 1')

title('2D state map')
end

end