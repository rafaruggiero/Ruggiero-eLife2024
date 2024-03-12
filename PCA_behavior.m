% PCA OF BEHAVIORAL DATA
% 
% 
% 
% Authors: Rafael Naime Ruggiero, 2021 

%Load data matrix
load('Matriz_corr.mat')

%Zscore data matrix
Xz = nanzscore(X,2);

%Reordering the matrix.
Xz = Xz(:,[2,1,3,4,5,6]);

%PCA
[coeff, score, latent, tsquared, explained, mu] = pca(Xz);

coeff(:,1)=coeff(:,1)*-1;
score(:,1)=score(:,1)*-1;

[rho,~] = corr([score(:,1) Xz]);


%Plot the projections and loadings of the PCA
f5=figure;
set(f5,'color','w');
subplot(3,2,[1,5])
scatter(score(1:11,1),score(1:11,2), 60,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[1 1 1])
hold on
scatter(score(12:end,1),score(12:end,2),60,'filled','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[219 78 78]./256)
xlabel('PC1')
ylabel('PC2')

subplot(3,2,2)
plot(explained,'k','LineWidth',2)
xlabel('PCs')
ylabel('Explained Variance')

names={'Distance'; 'Errors'; 'Startle';'PPI71';'PPI777'; 'PPI83' };
empty={};
subplot(3,2,4)

bar(rho(2:end,1),'FaceColor',[226 110 110]./256,'EdgeColor',[0 0 0]), hold on
set(gca,'xticklabel',empty)
line([0 7],[0.5 0.5],'Color','black','LineStyle','--')
line([0 7],[-0.5 -0.5],'Color','black','LineStyle','--')
ylabel('Loadings (PC1)')


[rho,~] = corr([score(:,2) Xz]);

subplot(3,2,6)
bar(rho(2:end,1), 'FaceColor',[226 110 110]./256,'EdgeColor',[0 0 0]), hold on
line([0 7],[0.5 0.5],'Color','black','LineStyle','--')
%line([0 7],[-0.5 -0.5],'Color','black','LineStyle','--')
ylabel('Loadings (PC2)')
set(gca,'xticklabel',names)

% figure
% biplot(coeff(:,1:2), 'Scores', score(1:11,1:2), 'Color','b','Marker','o', 'varlabels',names); hold on
% biplot(coeff(:,1:2), 'Scores', score(12:end,1:2), 'Color','r','Marker','o');
%
figure

plot([mean(score(1:11,1)), mean(score(1:11,2))],'bo'), hold on
plot([mean(score(12:end,1)), mean(score(12:end,2))],'ro')
hold on
errorbar([mean(score(1:11,1)), mean(score(1:11,2))],[std(score(1:11,1),[],1)./sqrt(11), std(score(1:11,2),[],1)./sqrt(11)],'b')
xlim([.5 2.5])
errorbar([mean(score(12:end,1)), mean(score(12:end,2))],[std(score(12:end,1),[],1)./sqrt(14), std(score(12:end,2),[],1)./sqrt(14)],'r')
xlim([.5 2.5])
ylabel('PC SCore')
xlabel('Component number')

%Statistical test of the scores of the PCs for each group
[H,P,stats,ci] = ttest2(score(1:11,1),score(12:end,1))
[H,P,stats,ci] = ttest2(score(1:11,2),score(12:end,2))

