% CALCULATES SPECTRAL ESTIMATES PER RAT
% Loads cell_block data
% Calculate power and coherence estimates
% Saves time-frequency data
% 
% Authors: Danilo Benette Marques and Rafael Naime Ruggiero

%Define subject
ratoName = 'R47_HG';
segName = '30s_dividido';
grupoName = 'Grupo A'; 

tic

%Reads variable cell_block
variaveis=who;
for v = 1:numel(variaveis)
    idx_cell(v) = ~isempty(findstr(variaveis{v},'cell_block'));
end
idx_cell = find(idx_cell);

data_block = eval(variaveis{idx_cell});

%Periods' indices
Nper = size(data_block,2);
for per = 1:Nper
    per_cell{per} = ones(size(data_block{4,per}))*per;
end

%Exclude empty periods
idx_excluir = find(cellfun(@isempty,data_block(4,:)));
data_block(:,idx_excluir) = [];
per_cell(:,idx_excluir) = [];
Nper = size(data_block,2);

% %% Divide segments and repeat classified states
dividir = 'y'
if dividir == 'y'
    newsegL = 3000;

    oldsegL = size(data_block{1,1},2);
    repsegnew = oldsegL/newsegL;

    for per = 1:size(data_block,2)
    for signalidx = 1:3
    data_block{signalidx,per} = reshape(data_block{signalidx,per}', newsegL,[])';
    end
    data_block{4,per} = repmat(data_block{4,per},1,repsegnew)';
    data_block{4,per} = data_block{4,per}(:);

    per_cell{1,per} = repmat(per_cell{1,per},1,repsegnew)';
    per_cell{1,per} = per_cell{1,per}(:);
    end
end

% %%
%Concatenate signals
PFC_cell = data_block(1,:);
HPC_cell = data_block(2,:);  
emg_cell = data_block(3,:);
state_cell = data_block(4,:);

state_labels = {'REM','SWS','WAKE','ACT'};

% %% Calculate spectral estimates per period per segment
Fs = 1000;

L = 1*Fs;
noverlap = .5*L;

pf = 0.5; %freq. precision
% f = linspace(.5,55,50);
f = .5:pf:400;

%downsample
downsample = 'n'
if downsample == 'y'
ds = 2; %1000 Hz --> 500 Hz (divide por 2)
    Fs = Fs/ds;
    L = 3*Fs;
    noverlap = .5*Fs;
    
    for per = 1:Nper
        for j = 1:size(PFC_cell{per},1)
           disp(['Downsampling do ' num2str(per) ' da janela ' num2str(j) ]) 
            PFC_ds{per}(j,:) = decimate(PFC_cell{per}(j,:),ds);
            HPC_ds{per}(j,:) = decimate(HPC_cell{per}(j,:),ds);
            emg_ds{per}(j,:) = decimate(HPC_cell{per}(j,:),ds);
            clc
        end
    end
    PFC_cell= PFC_ds;
    HPC_cell = HPC_ds;
    emg_cell = emg_ds;
    clear *_ds
end

%LFP
for per = 1:Nper
    for j = 1:size(PFC_cell{per},1)
       disp(['Calculando espectros do período ' num2str(per) ' da janela ' num2str(j) ]) 
        
    [PFC_pxx{per}(:,j),f] = pwelch(PFC_cell{per}(j,:),hamming(L),noverlap,f,Fs);
    [HPC_pxx{per}(:,j),f] = pwelch(HPC_cell{per}(j,:),hamming(L),noverlap,f,Fs);
    [COH{per}(:,j),f] = mscohere(HPC_cell{per}(j,:),PFC_cell{per}(j,:),hamming(L),noverlap,f,Fs); %arrumar função para input f

    clc
    end
end

%EMG
for per = 1:Nper
    
    emg_filt{per} = eegfilt(emg_cell{per},Fs,30,0);

    for j = 1:size(emg_cell{per},1);
       disp(['Calculando emg do período ' num2str(per) ' da janela ' num2str(j) ]) 
        
       [emg_pxx{per}(:,j),f_emg] = pwelch(emg_filt{per}(j,:),hamming(L),noverlap,[],Fs);
         
    clc
    end
    
    f_emg = find(f_emg>=30 & f_emg<=55 | f_emg>=65 & f_emg<=250);
    emg_rms{per} = sum(emg_pxx{per}(f_emg,:),1);
%     emg_rms{per} = sqrt(mean(emg_pxx{per}(f_emg,:),1));
end

% %% Concatenate all periods
PFC_all = cell2mat(PFC_pxx);
HPC_all = cell2mat(HPC_pxx);
COH_all = cell2mat(COH);
emg_all = cell2mat(emg_rms);
state_all = cell2mat(state_cell');
per_all = cell2mat(per_cell');

% %% Save
eval([ratoName '_PFCtf = PFC_all;']);
eval([ratoName '_HPCtf = HPC_all;']);
eval([ratoName '_COHtf = COH_all;']);
eval([ratoName '_emg = emg_all;']);
eval([ratoName '_states_idx = state_all;']);
eval([ratoName '_periods_idx = per_all;']);

eval(['clearvars -except ' ratoName '* ratoName grupoName segName'])

%Save tf_data in path
save([cd '\tf_data'],'-regexp', '^(?!(grupoName|ratoName|segName)$).')

clc
toc
disp(['Just ran rat ' ratoName ' do ' grupoName])
clear all

