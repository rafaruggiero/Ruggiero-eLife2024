% GROUP COMPARISONS ON SPECTRAL ESTIMATES AND SLEEP/WAKE STATES
% 
% Organization and selection
%   specific periods
%   concatenate states
% 
% Power and Coherence
%   per frequency
%   per band
%   per band peak
% 
% 2D statemap
%   states and transitions
%   automatic clustering
%   
% Similarity between states
%   average euclidean distance
% 
% Load all tf_data extraced from analysis_lfpstates_rat
% 
% Authors: Danilo Benette Marques and Rafael Naime Ruggiero

clearvars -except *tf *_emg *_states* *_periods*
variaveis=who;
for rato=1:numel(variaveis)
%rat name
    idx_=findstr(variaveis{rato},'_');
    idx_=idx_(1);
ratoName=variaveis{rato}(1:idx_-1); clear idx_  

%concatenate in cells
if ~isempty(findstr(variaveis{rato},'PFC'))
    PFC_tfdata{rato,1}=eval(variaveis{rato}); 
elseif ~isempty(findstr(variaveis{rato},'HPC'))
    HPC_tfdata{rato,1}=eval(variaveis{rato});  
elseif ~isempty(findstr(variaveis{rato},'COH')) %COH ou MPC!!
    COH_tfdata{rato,1}=eval(variaveis{rato}); 
elseif ~isempty(findstr(variaveis{rato},'emg'))
    EMG{rato,1}=eval(variaveis{rato});  
elseif ~isempty(findstr(variaveis{rato},'states_idx'))
    states_idx{rato,1}=eval(variaveis{rato});  
elseif ~isempty(findstr(variaveis{rato},'periods_idx'))
    periods_idx{rato,1}=eval(variaveis{rato});  
end

ratoNames{rato,1}=ratoName;

end

%concatenate and exclude empty
if exist('PFC_tfdata')
PFC_tfdata(find(cellfun(@isempty,PFC_tfdata)))=[]; end
if exist('HPC_tfdata')
HPC_tfdata(find(cellfun(@isempty,HPC_tfdata)))=[]; end
if exist('COH_tfdata')
COH_tfdata(find(cellfun(@isempty,COH_tfdata)))=[]; end
if exist('EMG')
EMG(find(cellfun(@isempty,EMG)))=[]; end
if exist('states_idx')
states_idx(find(cellfun(@isempty,states_idx)))=[]; end
if exist('periods_idx')
periods_idx(find(cellfun(@isempty,periods_idx)))=[]; end

ratoNames = unique(ratoNames);

%group color
% groupcolor = [0 0 0];
groupcolor = [1 0 0];

states = 1:4;
state_labels = {'REM','SWS','WAKE','ACT'};
% state_labels = {'1','2','3','4','5'};
state_colors = [.9 0 0; .1 .1 .9; 0 .5 0; 0 .9 0];  

%frequency vector
f = linspace(.5,55,50); 
f = .5 : .5 : 55; %divided
% f = .5 : .5 : 200; %

% %% Exclude 0s (noise "states")
for rato = 1:numel(ratoNames)
    PFC_tfdata{rato}(:,states_idx{rato}==0) = [];
    HPC_tfdata{rato}(:,states_idx{rato}==0) = [];
    COH_tfdata{rato}(:,states_idx{rato}==0) = [];
    EMG{rato}(states_idx{rato}==0) = [];
    periods_idx{rato}(states_idx{rato}==0) = [];
    states_idx{rato}(states_idx{rato}==0) = [];
end

% %% STATES SIMILARITY
% Define states to measure similarity
stateA = 1; %REM
stateB = 4; %ACT, join WAKE?

% %% Select periods
clear idx_per*
periods_selected = [1:5];

for rato = 1:numel(ratoNames) 
    
%Selected periods
idx_per = 0;
    for per = periods_selected
    idx_per = idx_per+1;
    idx_per_selected(:,idx_per) = periods_idx{rato}==periods_selected(idx_per);
    end
    idx_per_selected = sum(idx_per_selected,2);
    idx_per_excluir = find(idx_per_selected==0);
    
%First periods
qntsPer = 1; %+N periods (1 = 2pers, 4 = 5pers)
    %first period
    idx_per_selected = periods_idx{rato} >= min(periods_idx{rato}) & periods_idx{rato} <= min(periods_idx{rato}+qntsPer); 
%     %last periods
%     idx_per_selected = periods_idx{rato} >= max(periods_idx{rato})-qntsPer & periods_idx{rato} <= max(periods_idx{rato}+1); %penúltimo +1 período
    idx_per_excluir = find(idx_per_selected==0);

%select data
PFC_tfdata{rato}(:,idx_per_excluir) = [];
HPC_tfdata{rato}(:,idx_per_excluir) = [];
COH_tfdata{rato}(:,idx_per_excluir) = [];
EMG{rato}(idx_per_excluir) = [];
states_idx{rato}(idx_per_excluir) = [];
periods_idx{rato}(idx_per_excluir) = [];

clear idx_per*
end

%% New states
% Join or create states (e.g. WAKE+ACT = WAKE)
for rato = 1:numel(ratoNames)  
states_old = states_idx{rato};
periods_rato = periods_idx{rato};

%substitute states
states_new = states_old+10;
states_new(find(states_new==10)) = 0; %noise
states_new(find(states_new==11)) = 1; %REM
states_new(find(states_new==12)) = 2; %SWS
states_new(find(states_new==13)) = 3; %WAKE
states_new(find(states_new==14)) = 3; %ACT
% states_new(find(states_new>=3 & periods_rato==min(periods_rato))) = 4; %WAKE-N (novelty, 1st period)

% states_new(find(states_new>=13 & periods_rato~=max(periods_rato))) = 3; %WAKE-H (habituated, last period)
% states_new(find(states_new>=13 & periods_rato==min(periods_rato))) = 4; %WAKE-N (novelty, 1st period)

states_idx{rato} = states_new;
end

state_labels = {'REM','SWS','WAKE','WAKE-N'};
state_colors = [.9 0 0; .1 .1 .9; 0 .5 0; 0 .9 0;];  

states = 1:3;
stateA = 1;
stateB = 3;

clear states_old states_new;

%% Rats' spectrograms
for rato = 1:numel(ratoNames)
    figure('color','w'),
    subplot(3,1,1)
    imagesc([],f,HPC_tfdata{rato})
    axspct(1)=gca;
        axis xy
%         colorbar
        colormap jet
        caxis([0 100])
        ylabel('Frequency (Hz)')
        title(ratoNames{rato})
    subplot(3,1,2)
    imagesc([],f,PFC_tfdata{rato})
    axspct(2)=gca;
        axis xy
%         colorbar
        colormap jet
        caxis([0 100])
      ylabel('Frequency (Hz)')
      
    subplot(3,1,3)
    plot(EMG{rato},'r')
        axspct(3)=gca;
        ylim([0 mean(EMG{rato})+2*std(EMG{rato})])
        ylabel('EMG')
        
    linkaxes(axspct,'x')   
    xlim([1 length(states_idx{rato})])
    linkaxes(axspct(1:2),'y')   
    
    pause(.1)
end

%% 2D state map
excluirIS = 'y'

fig2d = figure('color','w');
for rato = 1:numel(ratoNames)
    statemapdata(:,:,1) = HPC_tfdata{rato}; 
    statemapdata(:,:,2) = PFC_tfdata{rato}; 
%     statemapdata(:,:,3) = COH_tfdata{rato}; 
    
[ratio1, ratio2,varexp1,varexp2] = statemap2D(statemapdata,f);

% figure,
subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
for state = states
hold on,scatter(ratio2(states_idx{rato}==state),ratio1(states_idx{rato}==state),'filled','markerfacecolor',state_colors(state,:),'sizedata',5)
end    
xlabel('Ratio 2')
ylabel('Ratio 1')
title(ratoNames{rato})

ratios12{rato,1} = [ratio1 ratio2];
varexp12(rato,:) = [varexp1(1) varexp2(1)];

clear statemapdata
end
legend(state_labels{states},'location','southeast');
% close

% %% Exclude transitions (Statemap 2D)

figts=figure('color','w');
fig2dnotjs=figure('color','w');
for rato = 1:numel(ratoNames)

trajspeed = diff(ratios12{rato},1,1);
trajspeed = sqrt(trajspeed(:,1).^2 + trajspeed(:,2).^2);
trajspeed = [0 ; trajspeed];

[tjsy,tjsx] = hist(trajspeed,20); [~,imax] = max(tjsy); 

TrajSpeed{rato} = trajspeed;
TrajSpeedmean(rato,1) = tjsx(imax);
TrajSpeedmean(rato,1) = mean(trajspeed);

set(0,'currentfigure',figts);
% hold on,bar(tjsx,tjsy,1)
hold on,plot(tjsx,tjsy)
    legend(ratoNames)
    title('Trajectory speed distribution')
    ylabel('Counts')
    xlabel('Trajectory speed')

%Threshold for trajectory speed
tjsthr = 0.1;
% tjsthr = max(trajspeed)/10;
% tjsthr = 2*mean(trajspeed);
% tjsthr = mean(trajspeed) + std(trajspeed);
% tjsthr = 2*tjsx(imax);
% tjsthr
TrajSpeedThr(rato) = tjsthr;

idx_notraj = find(trajspeed < tjsthr);

Nnotraj(rato) = length(idx_notraj);
Ntraj(rato) = size(ratios12{rato},1) - Nnotraj(rato);
Ptraj(rato) = 100*Ntraj(rato)/(Ntraj(rato)+Nnotraj(rato));

%2D clean
set(0,'currentfigure',fig2dnotjs)
subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)

scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),'filled','markerfacecolor',[.9 .9 .9],'sizedata',10)

for state = states
hold on,scatter(ratios12{rato}(intersect(idx_notraj,find(states_idx{rato}==state)),2),...
    ratios12{rato}(intersect(idx_notraj,find(states_idx{rato}==state)),1),...
    'filled','markerfacecolor',state_colors(state,:),'sizedata',10)
end
xlabel('Ratio 2')
ylabel('Ratio 1')
% legend(state_labels{states},'location','southeast');
title(ratoNames{rato})

%Exclude Intermediate States (IS; transitions)
    if excluirIS == 'y'
        ratios12{rato} = ratios12{rato}(idx_notraj,:);
        TrajSpeed{rato} = TrajSpeed{rato}(idx_notraj);
        states_idx{rato} = states_idx{rato}(idx_notraj);

        PFC_tfdata{rato} = PFC_tfdata{rato}(:,idx_notraj);
        HPC_tfdata{rato} = HPC_tfdata{rato}(:,idx_notraj);
        COH_tfdata{rato} = COH_tfdata{rato}(:,idx_notraj);
        EMG{rato} = EMG{rato}(idx_notraj);
        periods_idx{rato} = periods_idx{rato}(idx_notraj);
    end
end
legend({'IS',state_labels{states}},'location','southeast');
% close
% close

%% "Standard" sleep/wake estimates (EMG and theta/delta ratio)
figure('color','w')
for rato = 1:numel(ratoNames)
    statemapdata(:,:,1) = HPC_tfdata{rato}; 
%     statemapdata(:,:,2) = PFC_tfdata{rato}; 
%     statemapdata(:,:,3) = COH_tfdata{rato};
[thetadelta,emg] = thetadeltaemg(statemapdata,EMG{rato},f);
ThetaDelta{rato,1} = thetadelta;

subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
for state = states
hold on,scatter(emg(states_idx{rato}==state),thetadelta(states_idx{rato}==state),'filled','markerfacecolor',state_colors(state,:),'sizedata',10)

ThetaDelta_meanstate(rato,state) = mean(thetadelta(states_idx{rato}==state));
EMG_meanstate(rato,state) = mean(emg(states_idx{rato}==state));

end
xlabel('EMG')
ylabel('Theta/Delta')
title(ratoNames{rato})

% figure,
% subplot(2,1,1),hist((thetadelta),1000);
%     title(ratoNames{rato})
%     xlabel('Theta/delta')
% subplot(2,1,2),hist(emg,1000);
%     xlabel('EMG')
clear statemapdata
end
legend(state_labels{states});

results_thetadelta = ThetaDelta_meanstate'; open results_thetadelta
results_emg = EMG_meanstate'; open results_emg
% close all

%% State-maps automatic GMM clustering
clear Pgmm *Pthr tbl *chi2

k = 3; 

figgmdist = figure('color','w');
figgmk = figure('color','w');
figpgmm = figure('color','w');
for rato = 1:numel(ratoNames)

datatoclust = ratios12{rato};
% datatoclust = [PFC_band{rato,1}' ThetaDelta{rato}' EMG{rato}']; datatoclust = zscore(datatoclust,1);
% datatoclust = [ratios12{rato} EMG{rato}']; datatoclust = zscore(datatoclust,1);
    
% idxK = kmeans(datatoclust,k,'distance','sqeuclidean','maxiter',1000,'replicates',10);
% idxK = dbscan(datatoclust,.1,100);
% Z = linkage(datatoclust); idxK = cluster(Z,'maxclust',k);
% figure,gscatter(ratios12{rato}(:,2),ratios12{rato}(:,1),idxK)

GMModel = fitgmdist(datatoclust,k,'replicates',10,...
    'CovarianceType','full','SharedCovariance',true,...
    'ProbabilityTolerance',1e-8,'Regularizationvalue',0,...
    'options',statset('maxiter',1000,'display','final'));

    Pgmm{rato,1} = posterior(GMModel,datatoclust);
    Pthr = 0.80;
    idx_Pthr = find(any(Pgmm{rato}>Pthr,2));
    
    set(0,'currentfigure',figpgmm)    
    subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
    hold on
    scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),'filled','markerfacecolor',[.9 .9 .9],'sizedata',1)
    scatter(ratios12{rato}(idx_Pthr,2),ratios12{rato}(idx_Pthr,1),'filled','k','sizedata',1)
    xlabel('Ratio 2'), ylabel('Ratio 1')
    
idxK = cluster(GMModel,datatoclust);

% figure,scatter(ratios12{rato}(:,1),ratios12{rato}(:,2),'filled','k','sizedata',1)

set(0,'currentfigure',figgmdist)
subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
    hold on,
    scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),'filled','markerfacecolor',[.9 .9 .9],'sizedata',1)
    scatter(GMModel.mu(:,2),GMModel.mu(:,1),'sr')
    % gmPDF = @(x1,x2)reshape(pdf(GMModel,[x1(:) x2(:)]),size(x1));
    gmPDF = @(y,x)arrayfun(@(x,y)pdf(GMModel,[x(:) y(:)]),x,y);
    axgm = gca;
    fcontour(gmPDF,[axgm.XLim axgm.YLim],'meshdensity',200)
    hold off
    title(ratoNames{rato})
    xlabel('Ratio 2'), ylabel('Ratio 1')
    pause(.1)

set(0,'currentfigure',figgmk)
subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
    gscatter(ratios12{rato}(:,2),ratios12{rato}(:,1),idxK); legend off
    title(ratoNames{rato})
    xlabel('Ratio 2'), ylabel('Ratio 1')
    pause(.1)

% clustervsstates 
[tbl(:,:,rato),chi2(rato),pchi2(rato)] = crosstab(states_idx{rato},idxK);    
end

% figure,plot(chi2); title('Chi-squared')
% figure,hist(pchi2,10); title('p-value (chi-squared)')

%% SPECTRAL ESTIMATES (PSDs)

%Calculate relative power?
relativa = 'y'
if relativa == 'y'
flim = [f(1) 55]; 
for rato = 1:numel(ratoNames)
    PFC_tfdata{rato} = relpower(PFC_tfdata{rato},f,flim);
    HPC_tfdata{rato} = relpower(HPC_tfdata{rato},f,flim);
end    
end  

%% Concatenate PSD from all periods and states
for rato = 1:numel(ratoNames) 
    for period = unique(periods_idx{rato})'
        for state = states
        PFC_meanstateper(:,state,period,rato) = mean(PFC_tfdata{rato}(:,find(periods_idx{rato}==period & states_idx{rato}==state)),2);
        HPC_meanstateper(:,state,period,rato) = mean(HPC_tfdata{rato}(:,find(periods_idx{rato}==period & states_idx{rato}==state)),2);
        COH_meanstateper(:,state,period,rato) = mean(COH_tfdata{rato}(:,find(periods_idx{rato}==period & states_idx{rato}==state)),2);
        end
    end
end

% state = 1;
% figure,plot(f,squeeze(mean(HPC_meanstateper(:,state,:,:),4)),'linewidth',2)
% figure,plot(f,squeeze(mean(PFC_meanstateper(:,state,:,:),4)),'linewidth',2)
% figure,plot(f,squeeze(mean(COH_meanstateper(:,state,:,:),4)),'linewidth',2)

% %%
%Calculate average PSD per state
for rato = 1:numel(ratoNames)
    for state = states
    PFC_meanstate(:,state,rato) = mean(PFC_tfdata{rato}(:,states_idx{rato}==state),2);
    HPC_meanstate(:,state,rato) = mean(HPC_tfdata{rato}(:,states_idx{rato}==state),2);
    COH_meanstate(:,state,rato) = mean(COH_tfdata{rato}(:,states_idx{rato}==state),2);
    end
end
 
%Figures of states' averages
figure('color','w'),

% subplot(3,1,1)
figure
for state = states
    hold on,plot(f,nanmean(PFC_meanstate(:,state,:),3),'color',state_colors(state,:),'linewidth',2)
end
    title('PFC')
    ylabel('Power (mV{^2}/Hz)')
    xlabel('Frequency (Hz)')
    xlim([.5 200])
    
% subplot(3,1,2)
figure
for state = states
    hold on,plot(f,nanmean(HPC_meanstate(:,state,:),3),'color',state_colors(state,:),'linewidth',2)
end
    title('HPC')
    ylabel('Power (mV{^2}/Hz)')
    xlabel('Frequency (Hz)')
    xlim([.5 200])
    
% subplot(3,1,3)
figure
for state = states
    hold on,plot(f,nanmean(COH_meanstate(:,state,:),3),'color',state_colors(state,:),'linewidth',2)
end
    title('Coherence')
    ylabel('Coherence')
    xlabel('Frequency (Hz)')
    xlim([.5 200])
    legend(state_labels{states});
    
% close all
 
%% Figures of PSDs' mean and error
states=1:4;
groupcolor = [1 0 0];
    
%PFC  
figure('color','w')
for state = states
    subplot(1,max(states),state)
    boundedline(f,nanmean(PFC_meanstate(:,state,:),3),nanstd(PFC_meanstate(:,state,:),[],3)/sqrt(size(PFC_meanstate,3)),'cmap',groupcolor);
%     plot(f,squeeze(PFC_meanstate(:,state,:))); set(gca,'xscale','linear','yscale','linear')
    xlabel('Frequency (Hz')
    ylabel('Power (mV{^2}/Hz)')
    title(['PFC ' state_labels{state}])
end
    
%HPC
figure('color','w')
for state = states
    subplot(1,max(states),state)
    boundedline(f,nanmean(HPC_meanstate(:,state,:),3),nanstd(HPC_meanstate(:,state,:),[],3)/sqrt(size(HPC_meanstate,3)),'cmap',groupcolor);
%     plot(f,squeeze(HPC_meanstate(:,state,:))); set(gca,'xscale','linear','yscale','linear')
    xlabel('Frequency (Hz')
    ylabel('Power (mV{^2}/Hz)')
    title(['HPC ' state_labels{state}])
end

%HPC-PFC
figure('color','w')
for state = states
    subplot(1,max(states),state)
    boundedline(f,nanmean(COH_meanstate(:,state,:),3),nanstd(COH_meanstate(:,state,:),[],3)/sqrt(size(COH_meanstate,3)),'cmap',groupcolor);
%     plot(f,squeeze(COH_meanstate(:,state,:))); set(gca,'xscale','linear','yscale','linear')
    xlabel('Frequency (Hz')
    ylabel('Coherence')
    title(['HPC-PFC ' state_labels{state}])
end

states = 1:4;
% close all


%% Band power
clear *band*

%Bands' limits
band_labels = {'delta','theta','gamma'};
bandslims(1,:) = [1 4]; %delta
bandslims(2,:) = [4.5 9]; %theta
bandslims(3,:) = [30 50]; %gamma
bands = 1:size(bandslims,1);

%Per segment
for rato = 1:numel(ratoNames)
    
    for band = bands
        PFC_band{rato,band} = mean(PFC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),:),1);
        HPC_band{rato,band} = mean(HPC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),:),1);
        COH_band{rato,band} = mean(COH_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),:),1);
    end
    
figband=figure;
sbpidx=0;
    for state = states
        for band = bands
        PFC_bandstate{rato,band,state} = mean(PFC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1);
        HPC_bandstate{rato,band,state} = mean(HPC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1);
        COH_bandstate{rato,band,state} = mean(COH_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1);

        sbpidx=sbpidx+1;
        set(0,'currentfigure',figband)
        if band==1
        subplot(max(states),max(bands),sbpidx),hist(PFC_bandstate{rato,band,state},100)
        title([ratoNames{rato} ' PFC ' band_labels{band} ' ' state_labels{state}]);
        else
        subplot(max(states),max(bands),sbpidx),hist(HPC_bandstate{rato,band,state},100)
        title([ratoNames{rato} ' HPC ' band_labels{band} ' ' state_labels{state}]);
        end
     end
    end
    %close
end

% %%
%Averages
for rato = 1:numel(ratoNames)
    for band = bands
        for state = states
        PFC_meanbandstate(rato,band,state) = mean( mean(PFC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1) );
        HPC_meanbandstate(rato,band,state) = mean( mean(HPC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1) );
        COH_meanbandstate(rato,band,state) = mean( mean(COH_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state),1) );
        
            for period = unique(periods_idx{rato})'
                PFC_meanbandstateper(rato,band,state,period) = mean( mean(PFC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
                HPC_meanbandstateper(rato,band,state,period) = mean( mean(HPC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
                COH_meanbandstateper(rato,band,state,period) = mean( mean(COH_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
            end
        end
    end
end


% %% Results
data_meanbandstate = COH_meanbandstate;
results_delta = squeeze(data_meanbandstate(:,1,:))'; open results_delta
% results_theta = squeeze(data_meanbandstate(:,2,:))'; open results_theta
% results_gamma = squeeze(data_meanbandstate(:,3,:))'; open results_gamma

state = 1;
data_meanbandstateper = COH_meanbandstateper;
% results_delta = squeeze(data_meanbandstateper(:,1,state,:))'; open results_delta
results_theta = squeeze(data_meanbandstateper(:,2,state,:))'; open results_theta
% results_gamma = squeeze(data_meanbandstateper(:,3,state,:))'; open results_gamma
% close all

%% Band peak estimates
%Bands
clear *bandpk*
band_labels = {'delta','theta'};
bandslims(1,:) = [1 4]; %delta
bandslims(2,:) = [4.5 11]; %theta
bands = 1:size(bandslims,1);

%Per segment
for rato = 1:numel(ratoNames)
    
    for band = bands
        PFC_bandpk{rato,band} = bandpeakpower(PFC_tfdata{rato},f,bandslims(band,:));
        HPC_bandpk{rato,band} = bandpeakpower(HPC_tfdata{rato},f,bandslims(band,:));
        COH_bandpk{rato,band} = bandpeakpower(COH_tfdata{rato},f,bandslims(band,:));
    end
    
    figbandpks=figure;
    sbpidx=0;
    for state = states
        for band = bands
        PFC_bandpkstate{rato,band,state} = PFC_bandpk{rato,band}(states_idx{rato}==state);
        HPC_bandpkstate{rato,band,state} = HPC_bandpk{rato,band}(states_idx{rato}==state);
        COH_bandpkstate{rato,band,state} = COH_bandpk{rato,band}(states_idx{rato}==state);
        
%         PFC_bandpkstate{rato,band,state} = bandpeakpower(PFC_tfdata{rato}(:,states_idx{rato}==state),f,bandslims(band,:));
%         HPC_bandpkstate{rato,band,state} = bandpeakpower(HPC_tfdata{rato}(:,states_idx{rato}==state),f,bandslims(band,:));
%         COH_bandpkstate{rato,band,state} = bandpeakpower(COH_tfdata{rato}(:,states_idx{rato}==state),f,bandslims(band,:));
        
        sbpidx=sbpidx+1;
        set(0,'currentfigure',figbandpks)
        if band==1
        subplot(max(states),max(bands),sbpidx),hist(PFC_bandpkstate{rato,band,state},100)
        title([ratoNames{rato} ' PFC ' band_labels{band} ' ' state_labels{state}]);
        else
        subplot(max(states),max(bands),sbpidx),hist(HPC_bandpkstate{rato,band,state},100)
        title([ratoNames{rato} ' HPC ' band_labels{band} ' ' state_labels{state}]);
        end
        end
        pause(.1)
    end
end

%Averages
for rato = 1:numel(ratoNames)
    for band = bands
        for state = states
        PFC_meanbandpkstate(rato,band,state) = nanmean( PFC_bandpkstate{rato,band,state} );
        HPC_meanbandpkstate(rato,band,state) = nanmean( HPC_bandpkstate{rato,band,state} );
        COH_meanbandpkstate(rato,band,state) = nanmean( COH_bandpkstate{rato,band,state} );

            for period = unique(periods_idx{rato})'
                PFC_meanbandpkstateper(rato,band,state,period) = mean( mean(PFC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
                HPC_meanbandpkstateper(rato,band,state,period) = mean( mean(HPC_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
                COH_meanbandpkstateper(rato,band,state,period) = mean( mean(COH_tfdata{rato}(f>=bandslims(band,1) & f<=bandslims(band,2),states_idx{rato}==state & periods_idx{rato}==period),1)) ;
            end
        end
    end
end

% %% Results
% data_meanbandpkstate = COH_meanbandpkstate;
% results_delta = squeeze(data_meanbandpkstate(:,1,:))'; open results_delta
% results_theta = squeeze(data_meanbandpkstate(:,2,:))'; open results_theta
% results_gamma = squeeze(data_meanbandpkstate(:,3,:))'; open results_gamma

state = 1;
data_meanbandpkstateper = COH_meanbandpkstateper;
% results_delta = squeeze(data_meanbandpkstateper(:,1,state,:))'; open results_delta
results_theta = squeeze(data_meanbandpkstateper(:,2,state,:))'; open results_theta
% results_gamma = squeeze(data_meanbandpkstateper(:,3,state,:))'; open results_gamma
% close all
aleluia

%% Correlation theta x emg
for rato = 1:numel(ratoNames)
    figure,
    for state = states
        HPC_thetaemgcorr(rato,state) = corr(EMG{rato}(states_idx{rato}==state)',HPC_bandstate{rato,2,state}','type','pearson');
    
    subplot(1,max(states),state)
    scatter(EMG{rato}(states_idx{rato}==state)',HPC_bandstate{rato,2,state}')
    title(['state ' num2str(state) ', rato ' ratoNames{rato}])
    end
end

results = HPC_thetaemgcorr(:,4); open results
% close all


plotcorr(EMG_meanstate(:,4),squeeze(HPC_meanbandstate(:,2,4)),1,'k')

%% Compare EMG
figure('color','w'),
hx = 0:1:1000;
EMGall = cell2mat(cellfun(@transpose,EMG,'un',false)); states_idx_all = cell2mat(states_idx);
[hyA]=hist(EMGall(states_idx_all==stateA),hx);
[hyB]=hist(EMGall(states_idx_all==stateB),hx);
hold on,bar(hx,hyA,1,'facecolor',state_colors(stateA,:));
hold on,bar(hx,hyB,1,'facecolor',state_colors(stateB,:));
    ylabel('Counts')
    xlabel('EMG RMS')
xlim([0 200])

figemghist = figure('color','w');
for rato = 1:numel(ratoNames)
subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
[hyA]=hist(EMG{rato}(states_idx{rato}==stateA),hx);
[hyB]=hist(EMG{rato}(states_idx{rato}==stateB),hx);
hold on,bar(hx,hyA,1,'facecolor',state_colors(stateA,:));
hold on,bar(hx,hyB,1,'facecolor',state_colors(stateB,:));
    ylabel('Counts')
    xlabel('EMG RMS')
title(ratoNames{rato})
xlim([0 200])
end

%% Relationships between 2D state-map and other variables
for rato = 1:numel(ratoNames)

figure('color','w'),

subplot(2,3,1)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,PFC_band{rato,1},'filled')
    colormap jet
    c=colorbar; ylabel(c,'Power (mV^{2}/Hz)')
    caxis([0 200])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' PFC delta'])
subplot(2,3,2)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,ThetaDelta{rato,1},'filled')
    colormap jet
    c=colorbar; ylabel(c,'Theta/delta')
    caxis([0 1])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' HPC theta/delta'])
subplot(2,3,3)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,EMG{rato},'filled')
    colormap jet
    c=colorbar; ylabel(c,'RMS')
    caxis([0 100])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' EMG'])
       
subplot(2,3,4)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,COH_band{rato,1},'filled')
    colormap jet
    c=colorbar; ylabel(c,'Coherence')
    caxis([0 .7])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' HPC-PFC delta coherence'])
subplot(2,3,5)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,COH_band{rato,2},'filled')
    colormap jet
    c=colorbar; ylabel(c,'Coherence')
    caxis([0 .7])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' HPC-PFC theta coherence'])
    
subplot(2,3,6)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,TrajSpeed{rato},'filled')
    colormap jet
    colorbar
    c=colorbar; ylabel(c,'Trajectory speed')
    caxis([0 .3])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato} ' Speed'])
    
    pause(.1)
end

%% Relationshops between 2D state-maps with single variable
figure('color','w')
for rato = 1:numel(ratoNames)

C_rato = COH_band{rato,1};   
clabel = 'Coherence';

subplot(round(sqrt(numel(ratoNames))),ceil(sqrt(numel(ratoNames))),rato)
scatter(ratios12{rato}(:,2),ratios12{rato}(:,1),10,C_rato,'filled')
    colormap jet
    c=colorbar;
    ylabel(c,clabel)
%     caxis([0 .7])
%     caxis([0 100])
    xlabel('Ratio 2')
    ylabel('Ratio 1')
    title([ratoNames{rato}])
end

%% Average pointwise Euclidean distance
% ratios
figure('color','w')
for rato = 1:numel(ratoNames)

    for period = unique(periods_idx{rato})'
        Dratiosper = pdist2(ratios12{rato}(states_idx{rato}==stateA & periods_idx{rato}==period,:),...
            ratios12{rato}(states_idx{rato}==stateB & periods_idx{rato}==period ,:),...
        'euclidean');
        Dratiosper = Dratiosper(:);
        
        [hy,hx]=hist(Dratiosper,20);
        hy = hy/sum(hy);
        [~,imax] = max(hy);
        
    MEANDISTPER(rato,period) = hx(imax);
    clear Dratiosper hx hy imax
    end
    % OBS: selecting different periods slightly change ratios12
% end

% state = 1;
% figure,plot(f,squeeze(mean(HPC_meanstateper(:,state,:,:),4)),'linewidth',2)
    Dratios = pdist2(ratios12{rato}(states_idx{rato}==stateA,:),ratios12{rato}(states_idx{rato}==stateB,:),...
        'euclidean');
    Dratios = Dratios(:);
       
    [hy,hx]=hist(Dratios,20);
    hy = hy/sum(hy);
    [~,imax] = max(hy);

    hold on
    plot(hx,hy)
%     plot(hx,hy,'linewidth',2,'color',groupcolor)
%     scatter(hx(imax),hy(imax),'s','markeredgecolor',groupcolor)
        xlabel('Statemap euclidean distance')
        ylabel('Probability (%)')
        title('Statemap distance distribution')
%     xlim([0 5])
    pause(.1)
        
    MEANDIST(rato,1) = hx(imax);
%     MEANDIST(rato,1) = median(Dratios);
end
legend(ratoNames)

results = MEANDIST; open results

%% Results for comparisons

% open PFC_deltapower
% open PFC_thetapower
% open PFC_gammapower
% 
% open HPC_deltapower
% open HPC_thetapower
% open HPC_gammapower
% 
% open COH_deltapower
% open COH_thetapower
% open COH_gammapower

%% Figure nice
cf = gcf;
ca = gca;

set(cf,'color','w')
set(ca,'color','w','xcolor','k','ycolor','k',...
    'linewidth',1,...
    'box','off',...
    'fontname','helvetica')

% legend('off')

% clabel='Power (mV^{2}/Hz)';
% clabel='Theta/delta';
% clabel='RMS';
% clabel='Coherence';
% clabel='Trajectory speed';
% c=colorbar; ylabel(c,clabel)
% set(c,'color','w','xcolor','k','ycolor','k',...
%     'linewidth',1,...
%     'fontname','helvetica')

%% Save figures

%% Nested functions
%RELATIVE POWER
function [RelPower] = relpower(data,f,flim)
f_interesse = f>=flim(1) & f<=flim(end);
RelPower = data./repmat(sum(data(f_interesse,:),1),size(data,1),1);
end

%BAND PEAK POWER
function [BandPeakPower] = bandpeakpower(data,f,bandlim)
f_band = find(f>=bandlim(1) & f<=bandlim(end));

for j = 1:size(data,2)
    [pkpwr,pkidx] = findpeaks(data(f_band,j));
    if ~isempty(pkpwr)
        [~,maxpkidx] = max(pkpwr);
    BandPeakPower(j) = pkpwr(maxpkidx);
    BandPeak(j) = f(f_band(pkidx(maxpkidx)));
    else
    BandPeakPower(j) = NaN;
    BandPeak(j) = NaN;
    end
       
end
end

% THETA/DELTA & EMG
function [thetadelta,emg] = thetadeltaemg(dado,emg,f);
%Define frequencies of ratios
f_theta = find(f>=4.5 & f<=9);
f_delta = find(f>=1 & f<=4);

%Calculate ratios per region
for region = 1:size(dado,3)
thetadelta(region,:) = mean(dado(f_theta,:,region),1) ./ sum(dado(f_delta,:,region),1);
end

%PCA of ratios from all regions
% [~,pcscore] = pca(thetadelta');
% thetadelta = pcscore(:,1);

%Hann window smoothing
% thetadelta = conv(thetadelta,hann(10),'same');
% emg = conv(emg,hann(10),'same')';

% datatoscore = [thetadelta emg];
% Z = linkage(datatoscore);
% scoring = cluster(Z,'maxclust',3);
end

% STATE MAP 2D
function [ratio1, ratio2,varexp1,varexp2] = statemap2D(dado,f);
%Define frequencies of ratios (as in Dzirasa et al., 2006)
f_ratio1n = find(f>=2 & f<=20);
f_ratio1d = find(f>=2 & f<=55);

f_ratio2n = find(f>=2 & f<=4.5);
f_ratio2d = find(f>=2 & f<=9);

%Calculate ratios per region 
for region = 1:size(dado,3);
ratio1(region,:) = sum(dado(f_ratio1n,:,region),1) ./ sum(dado(f_ratio1d,:,region),1);
ratio2(region,:) = sum(dado(f_ratio2n,:,region),1) ./ sum(dado(f_ratio2d,:,region),1);
end

%PCA of ratios from all regions
[~,pcscore,~,~,varexp1] = pca(ratio1');
ratio1 = pcscore(:,1);

[~,pcscore,~,~,varexp2] = pca(ratio2');
ratio2 = pcscore(:,1);

%Hann window smoothing
% originally of 20 secs (adjust according to length)
ratio1 = conv(ratio1,hann(10),'same');
ratio2 = conv(ratio2,hann(10),'same');
end
