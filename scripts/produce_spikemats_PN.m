function [] = produce_spikemats_PN(session_name)

% [] = produce_spikemats_PN(session_name)
% generate SPIKEMAT and SPIKES matrices for a given session
% -------------------------------------------------------------------------

%% add relevant paths

pars=set_pars_PN();
code_folder=pars.code_folder;
bitcode_folder=pars.bitcode_folder;
input_folder=pars.raw_data_folder;
addpath(code_folder)
addpath(input_folder)
addpath(bitcode_folder)

%% load bitcodes and times

load(['STIM_',session_name,'.mat'])
% STIM contains data (digital bitcode vector) every 66 ms (2 frames at
% 30 Hz) and corresponding bitcode onset times (or every 132 ms! if at 60 Hz).

load(['my_times_',session_name,'.mat'])

BITCODE=data;
% nu_onset=onset-0.033;  % to account for a shift of 1 frame in the anlg signal
nu_onset=onset;

% borders of time window in wich to take trial spike times
PRE_TIME=200/1000;
POST_TIME=1200/1000;
PRE_TIME_bis=800/1000;
POST_TIME_bis=1700/1000;

% borders of time window in wich to take trial spike count
WINDOW_START=50/1000;
WINDOW_STOP=900/1000;

% bitcode diversi tipi di stimolo
bcodes_GR=pars.stimPars.gratings_bitcodes;
bcodes_PL=pars.stimPars.plaids_bitcodes;
bcodes_FL=pars.stimPars.flashes_bitcodes;

%% load Klusta output
[ Spikes_g, Spikes_m, Spikes_n, goodc, muac, noisec, clabel, i_sp, clusters, masksidx ] = extract_spikes_64ch_PN(session_name);
SPIKES.myspikes=cat(2,Spikes_g, Spikes_m, Spikes_n);

clusters_good=clusters(clabel==2);
clusters_mua=clusters(clabel==1);
clusters_noise=clusters(clabel==0);

i_chang=NaN(1,length(clusters_good));
for i=1:length(clusters_good)
    i_chang(i)=mode(masksidx(i_sp==clusters_good(i)));
end
i_chanm=NaN(1,length(clusters_mua));
for i=1:length(clusters_mua)
    i_chanm(i)=mode(masksidx(i_sp==clusters_mua(i)));
end
i_chann=NaN(1,length(clusters_noise));
for i=1:length(clusters_noise)
    i_chann(i)=mode(masksidx(i_sp==clusters_noise(i)));
end

SPIKES.mychannel=cat(2,i_chang, i_chanm, i_chann);

% generate facke channel list for debug purposes
% SPIKES.mychannel=33*ones(1,length(goodc)+length(muac)+length(noisec));

%% bitcode vector artifacts removal

% find and replaces infinites (if any)
[c, ~]=ind2sub(size(nu_onset), find(nu_onset==inf));
istanzainf_on = c;
if isempty(istanzainf_on) == 0
    nu_onset(istanzainf_on) = [];
    BITCODE(istanzainf_on) = [];
end

% find and replaces white screens (if any)
[e, ~]=ind2sub(size(BITCODE), find(BITCODE==1023));
istanzaflash = e;
if isempty(istanzaflash)==0
    nu_onset(istanzaflash) = [];
    BITCODE(istanzaflash) = [];
end

% find and replaces black screens (if any)
[wi, ~]=ind2sub(size(BITCODE), find(BITCODE==0)); 
istanzazero = wi;
if isempty(istanzazero)==0
    nu_onset(istanzazero) = [];
    BITCODE(istanzazero) = [];
end

% find and replace aborted stimuli (isolated bitcodes) (if any)
istanzasingoloflash = [];

for h=[bcodes_GR,bcodes_PL]
    % find indeces of elements of bitcode vector with a given value
    [dio, ~]=ind2sub(size(BITCODE), find(BITCODE==h));
    
    for bibi=1:length(dio)
        sing_flash_idx=dio(bibi);
        % if this is not the last one...
        if sing_flash_idx ~= numel(BITCODE)
            % ... if it is different prom previous and next
            if BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx+1) && BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx-1)
                istanzasingoloflash = [istanzasingoloflash, sing_flash_idx];
                disp(['warning: idx # ', num2str(sing_flash_idx), ' of bcode # ', num2str(h), ' will be removed' ])
            end
          % if this is the last one...
        else
            % ... if it is different prom previous 
            if BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx-1)
                istanzasingoloflash = [istanzasingoloflash, sing_flash_idx];
                disp(['warning: idx # ', num2str(sing_flash_idx), ' of bcode # ', num2str(h), ' will be removed' ])
                
            end
        end
        clear sing_flash_idx
    end
    clear dio poi
    
end
BITCODE(istanzasingoloflash)=[];
nu_onset(istanzasingoloflash) = [];
clear c d e f wi fi

%% get onset and offset time for every bitcode

% create a cronologically ordered list of presented bitcodes
BitCodeS4 = unique_no_sort(BITCODE);
% exclude noise movies
BitCodeS = BitCodeS4(ismember(BitCodeS4,[bcodes_GR,bcodes_PL,bcodes_FL]));

trial_start=cell(zeros());
trial_stop=cell(zeros());
temp_All_Trial_NumberS=cell(zeros());

for bcode=BitCodeS
    
    % find indexes of elements of bitcode vector at the end of a trial 
    % (from the list of all bitcode onset of a particular condition)
    offIdx = find(diff(nu_onset(BITCODE==bcode))>0.5);
    % after an on ther is an off
    we=offIdx+1;
    % the first on is one
    onIdx=[1; we];
    
    tr_count=0;
    mycount=0;
    
    for bibi2=1:length(onIdx)
        
        % index of current trial onset from the list of all bitcode onset of a particaolar condition)
        kaka=onIdx(bibi2);
        
        tr_count=tr_count+1;
        mycount=mycount+1;
        
        % find current bitcode indexes
        [ki, ~]=ind2sub(size(BITCODE), find(BITCODE==bcode));
        
        % current trial onset index from the bitcode vector
        temp_onset_trial_idx=ki(kaka);
        temp_All_Trial_NumberS{bcode, bibi2} = temp_onset_trial_idx;
        % current trial onset
        temp_onset_trial = nu_onset(temp_onset_trial_idx);
        
        % DEVELOPMENT NOTE: METTERE UN IF CHE AGGIUSTA LA STIMULUS LENGTH IN BASE A CHE
        % STIMOLO (BITCODE) E QUESTO E' ANCHE IL POSTO GIUSTO IN CUI
        % SHIFTARE INDIETRO DI 33 o di 66 ms TUTTO IL VETTORE ONSET A
        % SECONDA CHE SIA NOISE O ALTRO ---> Done!
        
% store current trial onsets
if sum(bcode==bcodes_FL)
    stimulusength_FL=0.3;
    % account for bitcode shift at 30 Hz
    temp_onset_trial=temp_onset_trial-0.033; 
    % store current trial onsets
    trial_start{bcode,tr_count} = temp_onset_trial;
    % store current trial offsets
    trial_stop{bcode,tr_count} = temp_onset_trial+stimulusength_FL;
elseif sum(bcode==bcodes_PL) || sum(bcode==bcodes_GR)
    stimulusength_OTHER=0.9;
    % account for bitcode shift at 30 Hz
    temp_onset_trial=temp_onset_trial-0.033;
    % store current trial onsets
    trial_start{bcode,tr_count} = temp_onset_trial;
    % store current trial offsets
    trial_stop{bcode,tr_count} = temp_onset_trial+stimulusength_OTHER;
end
        
    end
end

STIM_START = [];
STIM_STOP = [];

% produce "global" (for every bitcode) onset and offset vector
for zz= BitCodeS
    temp_STIM_START = cell2mat(trial_start(zz,:));
    temp_STIM_STOP = cell2mat(trial_stop(zz,:));
    % concatenate successive bitcodes trial onset time vectors
    STIM_START = [temp_STIM_START, STIM_START];
    STIM_STOP = [temp_STIM_STOP, STIM_STOP];
end

%% recupera le spikes per ogni neurone e trial di ogni bitcode

NeuronS = 1:size(SPIKES.mychannel,2);
SPIKEtimes=cell(zeros());
SPIKEtimes_bis=cell(zeros());
SPIKEcounts=cell(zeros());
basalFR=zeros(1,max(NeuronS));

neurons=0 ;
% for every neuron
for nn=NeuronS
    
    tic
    % get spike timestamps
    T_TIMES=SPIKES.myspikes{nn};
    TIMES=T_TIMES/1000;
    
    neurons = neurons+1;
    grcount=0;
    isicount=0;
    my_bits=unique(BitCodeS);
    
    % and for every bitcode
    for i=my_bits
        grcount=grcount+1;
        
        % recover current bitcode trial onsets
        temp_my_trials=trial_start(i,1:end);
        temp_my_trials=cell2mat(temp_my_trials);
        
        % rename variables
        All_Trial_TimeS = temp_my_trials;
        
        % for every trial of current bitcode
        for tt=1:numel(All_Trial_TimeS)
            
            TIME_START_1=trial_start{i,tt}-PRE_TIME;
            TIME_END_1=trial_start{i,tt}+POST_TIME;
            TIME_WINDOW_ON=trial_start{i,tt}+WINDOW_START;
            TIME_WINDOW_OFF=trial_start{i,tt}+WINDOW_STOP;
            TIME_START_3=trial_start{i,tt}-PRE_TIME_bis;
            TIME_END_3=trial_start{i,tt}+POST_TIME_bis;
            % get spike to put in spike times matrix
            SPIKES_TAKEN_1=TIMES(TIMES<TIME_END_1 & TIMES>TIME_START_1)-TIME_START_1;
            % get spike to put in spike count matrix
            SPIKES_TAKEN_2=TIMES(TIMES<TIME_WINDOW_OFF & TIMES>TIME_WINDOW_ON)-TIME_START_1;
            % get spike to put in spike times matrix bis
            SPIKES_TAKEN_3=TIMES(TIMES<TIME_END_3 & TIMES>TIME_START_3)-TIME_START_3;
            
            SPIKEtimes{i,nn,tt}=SPIKES_TAKEN_1;  %% spike timestamps per trial (RASTER)
            SPIKEcounts{i,nn,tt}=numel(SPIKES_TAKEN_2);  %% spike counts per trial (RASTER)
            SPIKEtimes_bis{i,nn,tt}=SPIKES_TAKEN_3;  %% spike counts per trial (RASTER)
            
            basalFR(nn)=basalFR(nn)+numel(TIMES(TIMES<TIME_START_1 & TIMES>TIME_START_1-0.5));
            isicount=isicount+1;
        end
        
        clear NEURON PSTH_RASTER My_Neurons_PL PsthAndRaster_PL Channel First_Spike Last_Spike TT...
            tt tt0 First_Trial Last_Trial All_Spikes Trial_num T_num TIME_START TIME_END SPIKES_TAKEN...
            My_Spikes PSTH Shape Mean_Firing_Rate q All_Trial_TimeS All_Trial_NumberS
        
    end
    toc
    fprintf(['Neuron ',num2str(nn),' spikes processed \n'])
end
basalFR=basalFR./isicount*0.5;

%% plot response matrix for visual inspection of session quality

SPIKEmean=zeros(size(SPIKEtimes,1),size(SPIKEtimes,2));
SPIKEstd=zeros(size(SPIKEtimes,1),size(SPIKEtimes,2));
% produce mean and std firing rate matrix
for i=my_bits
    for nn=1:size(SPIKEtimes,2)
        SPIKEmean(i,nn)=mean(squeeze(cell2mat(SPIKEcounts(i,nn,:))))-basalFR(nn);
        SPIKEstd(i,nn)=std(squeeze(cell2mat(SPIKEcounts(i,nn,:))));
    end
end
% convert in Hz
SPIKEmean=SPIKEmean./(WINDOW_STOP-WINDOW_START);

tit=[session_name,' response matrix'];
tit=strrep(tit,'_','\_');
f = figure;
set(f,'Position',[10,10,1500,1000]);
imagesc(cat(1,repmat(basalFR,[12,1]),SPIKEmean)); colormap('paruly'); c=colorbar; xlabel('neuron number'); ylabel('bitcode number');
tempmat=SPIKEmean(:,1:numel(goodc));
caxis([min(tempmat(:)),max(tempmat(:))]);
% caxis([quantile(tempmat(:),0.02),quantile(tempmat(:),0.98)]);
set(get(c,'ylabel'),'String', 'mean firing rate (Hz)');
line([0 555],[11.5 11.5],'Color',[1,1,1])
line([0 555],[11.5+100 11.5+100],'Color',[1,1,1])
for kk=1:9
    line([0 555],[11.5+100+12*kk 11.5+100+12*kk],'Color',[1,1,1])
end
line([0 555],[11.5+350 11.5+350],'Color',[1,1,1])
Xticklabels = 0:5:size(SPIKEtimes,2);
Xticks = linspace(0, size(SPIKEmean, 2), numel(Xticklabels)+1);
set(gca,'XTick', Xticks,'XTickLabel', Xticklabels)
Yticklabels = fliplr(cat(2,0,0:12:555));
Yticks = 0:12:555+12;
set(gca, 'YTick', Yticks,'YTickLabel', flipud(Yticklabels(:)));
for kk=1:size(SPIKEtimes,2)
    line([kk kk],[0 12+554],'LineStyle','-.','Color',[0.5,0.5,0.5])
end
for kk=[length(goodc),length(goodc)+length(muac)]
    line([kk+0.5 kk+0.5],[0 12+554],'LineStyle','-','Color',[1,1,1],'LineWidth',2)
end

h = suptitle(tit);
set(gca, 'Visible', 'on');
set(h, 'Visible', 'on', 'FontSize', 15);
set(gcf, 'Color', 'w');
set(gcf, 'PaperPositionMode', 'auto')

imagename=[session_name,'response_matrix.jpg'];
saveas(f,imagename, 'jpg')
save(['SPIKEMAT_',session_name,'.mat'],'SPIKEtimes','SPIKEtimes_bis','SPIKEcounts','SPIKEmean','SPIKEstd','goodc','muac','noisec','basalFR')
save(['SPIKES_',session_name,'.mat'],'SPIKES','goodc','muac','noisec') 

end
