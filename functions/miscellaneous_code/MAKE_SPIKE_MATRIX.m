clear all
close all
clc

%% Set paths and folder

DaysOfRecording_LL = {'18_04_2016','15_04_2016','07_04_2016','25_03_2016','23_03_2016','21_03_2016','17_03_2016'};
Blocks_LL={7,[3,5],[5],4,[4],4,[5]};

DaysOfRecording_V1 = {'02_03_2016','22_02_2016','10_02_2016','29_01_2016','14_01_2016','12_01_2016'};
Blocks_V1={3,3,[6],2,3,[6,3]};

addpath /zocconasphys1/chronic_inv_rec/codes/GIULIO/SPIKE_TRIGGERED_ANALYSIS/NEWCODE/paruly
addpath /zocconasphys1/chronic_inv_rec/codes/GIULIO/SPIKE_TRIGGERED_ANALYSIS/NEWCODE
addpath /zocconasphys1/chronic_inv_rec/codes/NEW_FEDE_GIULIO
addpath /zocconasphys1/chronic_inv_rec/codes/GIULIO/SPIKE_TRIGGERED_ANALYSIS/NEWCODE/functions/
out_folder='/zocconasphys1/chronic_inv_rec/codes/GIULIO/SPIKE_TRIGGERED_ANALYSIS/Global_Results';
mkdir(out_folder);
cd(out_folder)


%% ---------------------------------- V1 ----------------------------------

for ii=1:length(DaysOfRecording_V1)
    for jj=1:length(Blocks_V1{ii})
        
        session_folder=['/zocconasphys1/chronic_inv_rec/Tanks/Giulio_Acute_Recording_', char(DaysOfRecording_V1{ii}),'/Block-', num2str(Blocks_V1{ii}(jj))];
        old=cd(session_folder);
        
        %% Carica dati bitcodes e tempi
        
        load('STIM.mat')
        % STIM contains data (digital bitcode vector) every 66 ms (2 frames at
        % 30 Hz) and corresponding bitcode onset times.
        % % figure; plot(data(1:10000));
        % % figure; plot(onset(1:10000));
        load('SPIKE.mat')
        load('my_times.mat')
        
        BITCODE=data;
        nu_onset=onset-0.033;  % to account for a shift of 1 frame in the anlg signal
        
        % borders of time window in wich to take trial spike times
        PRE_TIME=200/1000;
        POST_TIME=1200/1000;
        
        % borders of time window in wich to take trial spike count
        WINDOW_START=50/1000;
        WINDOW_STOP=900/1000;
        
        % bitcode diversi tipi di stimolo
        bcodes_noise = [1:80];
        bcodes_GR = [100:207];
        bcodes_GR_c = [250:267];
        bcodes_PL = [350:361];
        bcodes_FL = [555,556];
        
        % stimulus duration in seconds
        stimulusength=0.9;
        
        fprintf(['\nSession ',DaysOfRecording_V1{ii},' b',num2str(Blocks_V1{ii}(jj)),' data loaded\n'])
        fprintf(['\n-------------------------------\n'])
        
        %% Elimina artefatti da bitcode vector
        
        % trova ed elimina infiniti
        
        [c, d]=ind2sub(size(nu_onset), find(nu_onset==inf));
        istanzainf_on = c;
        if isempty(istanzainf_on) == 0
            nu_onset(istanzainf_on) = [];
            BITCODE(istanzainf_on) = [];
        end
        
        % trova ed elimina white screens
        [e, f]=ind2sub(size(BITCODE), find(BITCODE==1023));
        istanzaflash = e;
        if isempty(istanzaflash)==0
            nu_onset(istanzaflash) = [];
            BITCODE(istanzaflash) = [];
        end
        
        % trova ed elimina black screens
        [wi, fi]=ind2sub(size(BITCODE), find(BITCODE==0));
        istanzazero = wi;
        if isempty(istanzazero)==0
            nu_onset(istanzazero) = [];
            BITCODE(istanzazero) = [];
        end
        
        % trova ed elimina stimoli abortiti (bitcodes isolati)
        istanzasingoloflash = [];
        
        for h=[bcodes_GR,bcodes_PL]
            % trova indici di tutti gli elementi di bitcode vector con un certo valore
            [dio, poi]=ind2sub(size(BITCODE), find(BITCODE==h));
            
            for bibi=1:length(dio)
                sing_flash_idx=dio(bibi);
                % se non � l'ultimo che...
                if sing_flash_idx ~= numel(BITCODE)
                    % ... sia diverso dal precedente e dal successivo
                    if BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx+1) && BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx-1)
                        istanzasingoloflash = [istanzasingoloflash, sing_flash_idx];
                        disp(['warning: idx # ', num2str(sing_flash_idx), ' of bcode # ', num2str(h), ' will be removed' ])
                    end
                    % se non � l'ultimo che...
                else
                    % ... sia diverso dal precedente
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
        
        %% Recupera tempi di onset e offset per ogni trial di ogni bitcode
        
        % create a cronologically ordered list of presented bitcodes
        BitCodeS4 = unique_no_sort(BITCODE);
        % exclude nois movies
        BitCodeS = BitCodeS4(ismember(BitCodeS4,[bcodes_GR,bcodes_PL]));
        
        trial_start=cell(zeros());
        trial_stop=cell(zeros());
        temp_All_Trial_NumberS=cell(zeros());
        
        for bcode=BitCodeS;
            
            % find indexes of elements of bitcode vector at the end of a trial (from the list of all bitcode onset of a particaolar condition)
            offIdx = find(diff(nu_onset(BITCODE==bcode))>0.5);
            % after an on ther is an off
            we=offIdx+1;
            % the first on is one
            onIdx=[1; we];
            
            tr_count=0;
            mycount=0;
            temp_onset_trial_idx=0;
            temp_onset_trial=0;
            
            for bibi2=1:length(onIdx)
                
                % index of current trial onset from the list of all bitcode onset of a particaolar condition)
                kaka=onIdx(bibi2);
                
                tr_count=tr_count+1;
                mycount=mycount+1;
                
                % find current bitcode indexes
                [ki, li]=ind2sub(size(BITCODE), find(BITCODE==bcode));
                
                % current trial onset index from the bitcode vector
                temp_onset_trial_idx=ki(kaka);
                temp_All_Trial_NumberS{bcode, bibi2} = temp_onset_trial_idx;
                % current trial onset
                temp_onset_trial = nu_onset(temp_onset_trial_idx);
                
                % store current trial onsets
                trial_start{bcode,tr_count} = temp_onset_trial;
                % store current trial offsets
                trial_stop{bcode,tr_count} = temp_onset_trial+stimulusength;
                
            end
        end
        
        STIM_START = [];
        STIM_STOP = [];
        
        % produce "global" (for every bitcode)onset and offset vector
        for zz= BitCodeS;
            temp_STIM_START = cell2mat(trial_start(zz,:));
            temp_STIM_STOP = cell2mat(trial_stop(zz,:));
            % concatenate successive bitcodes trial onset time vectors
            STIM_START = [temp_STIM_START, STIM_START];
            STIM_STOP = [temp_STIM_STOP, STIM_STOP];
        end
        
        % create vectors containing every trial onset and offset in cronological order
        aha = sort(STIM_START);
        baha = sort(STIM_STOP);
        
        %% Recupera le spikes per ogni neurone e trial di ogni bitcode
        
        NeuronS = 1:size(SPIKES.mychannel,2);
        BitCodeS3 = BITCODE(ismember(BITCODE,[bcodes_PL,bcodes_GR]));
        BIT_START=BitCodeS3;
        NumTrials = length(bcodes_PL);
        ChannelS = SPIKES.mychannel;
        
        nu_onsets_gr=cell(zeros());
        nu_offsets_gr=cell(zeros());
        First_Spike = cell(zeros());
        Last_Spike = cell(zeros());
        
        temp_nu_onset_all = [];
        temp_nu_offset_all = [];
        
        
        Trial_num=cell(zeros());
        SPIKEtimes=cell(zeros());
        SPIKEcounts=cell(zeros());
        basalFR=zeros(1,max(NeuronS));
        
        neurons=0 ;
        % for every neuron
        for nn=NeuronS;
            
            tic
            % get spike timestamps
            T_TIMES=SPIKES.myspikes{nn};
            TIMES=T_TIMES/1000;
            
            neurons = neurons+1;
            grcount=0;
            isicount=0;
            PSTH=cell(zeros());
            PSTH_25=cell(zeros());
            my_bits=unique(BitCodeS);
            
            % and for every bitcode
            for i=my_bits
                grcount=grcount+1;
                trialcount = 0;
                tr_count=0;
                
                % recover current bitcode trial onsets
                temp_my_trials=trial_start(i,1:end);
                temp_my_trials=cell2mat(temp_my_trials);
                my_trials = numel(temp_my_trials);
                
                % recover current bitcode trial onset index from the bitcode vector
                temp_trial_numbers = temp_All_Trial_NumberS(i,1:end);
                temp_trial_numbers=cell2mat(temp_trial_numbers);
                trial_numbers = temp_trial_numbers;
                
                % rename variables
                All_Trial_TimeS = temp_my_trials;
                All_Trial_NumberS = temp_trial_numbers;
                STIM_TI=All_Trial_TimeS;
                TRIAL_TI=All_Trial_NumberS;
                
                % for every trial of current bitcode
                for tt=1:numel(All_Trial_TimeS)
                    
                    TIME_START_1=trial_start{i,tt}-PRE_TIME;
                    TIME_END_1=trial_start{i,tt}+POST_TIME;
                    TIME_WINDOW_ON=trial_start{i,tt}+WINDOW_START;
                    TIME_WINDOW_OFF=trial_start{i,tt}+WINDOW_STOP;
                    % get spike to put in spike times matrix
                    SPIKES_TAKEN_1=TIMES(TIMES<TIME_END_1 & TIMES>TIME_START_1)-TIME_START_1;
                    % get spike to put in spike count matrix
                    SPIKES_TAKEN_2=TIMES(TIMES<TIME_WINDOW_OFF & TIMES>TIME_WINDOW_ON)-TIME_START_1;
                    
                    SPIKEtimes{i,nn,tt}=SPIKES_TAKEN_1;  %% spike timestamps per trial (RASTER)
                    SPIKEcounts{i,nn,tt}=numel(SPIKES_TAKEN_2);  %% spike counts per trial (RASTER)
                    
                    basalFR(nn)=basalFR(nn)+numel(TIMES(TIMES<TIME_START_1 & TIMES>TIME_START_1-0.5));
                    isicount=isicount+1;
                end
                
                clear q All_Trial_TimeS All_Trial_NumberS
                clear NEURON PSTH_RASTER My_Neurons_PL PsthAndRaster_PL Channel First_Spike Last_Spike TT tt tt0 First_Trial Last_Trial All_Spikes Trial_num T_num TIME_START TIME_END SPIKES_TAKEN My_Spikes PSTH Shape Mean_Firing_Rate
                
            end
            toc
            fprintf(['Neuron ',num2str(nn),' spikes processed \n'])
        end
        basalFR=basalFR./isicount*0.5;
        
        %% Plot response matrix
        
        cd(old)
        SPIKEmean=zeros(size(SPIKEtimes,1),size(SPIKEtimes,2));
        SPIKEstd=zeros(size(SPIKEtimes,1),size(SPIKEtimes,2));
        % produce mean and std firing rate matrix
        for i=my_bits
            for nn=1:size(SPIKEtimes,2)
                SPIKEmean(i,nn)=mean(squeeze(cell2mat(SPIKEcounts(i,nn,:))))-basalFR(nn);;
                SPIKEstd(i,nn)=std(squeeze(cell2mat(SPIKEcounts(i,nn,:))));
            end
        end
        % convert in Hz
        SPIKEmean=SPIKEmean./(WINDOW_STOP-WINDOW_START);
        
        % plot figure
        tit=['V1 response matrix ',DaysOfRecording_V1{ii},' b',num2str(Blocks_V1{ii}(jj))];
        tit=strrep(tit,'_','\_');
        f = figure;
        set(f,'Position',[10,10,1500,1000]);
        imagesc(cat(1,repmat(basalFR,[12,1]),SPIKEmean)); colormap('paruly'); c=colorbar; caxis([0,30]); xlabel('neuron number'); ylabel('bitcode number');
        set(get(c,'ylabel'),'String', 'mean firing rate (Hz)');        
        line([0 361],[11.5 11.5],'Color',[1,1,1])
        line([0 361],[11.5+100 11.5+100],'Color',[1,1,1])
        for kk=1:9
        line([0 361],[11.5+100+12*kk 11.5+100+12*kk],'Color',[1,1,1])
        end
        line([0 361],[11.5+350 11.5+350],'Color',[1,1,1])         
        Xticklabels = 0:5:size(SPIKEtimes,2);
        Xticks = linspace(0, size(SPIKEmean, 2), numel(Xticklabels)+1);
        set(gca,'XTick', Xticks,'XTickLabel', Xticklabels)
        Yticklabels = fliplr(cat(2,0,0:12:361));
        Yticks = 0:12:361+12;
        set(gca, 'YTick', Yticks,'YTickLabel', flipud(Yticklabels(:)));
        for kk=1:size(SPIKEtimes,2)
        line([kk kk],[0 12+360],'LineStyle','-.','Color',[0.5,0.5,0.5])
        end
        h = suptitle(tit);
        set(gca, 'Visible', 'on');
        set(h, 'Visible', 'on', 'FontSize', 15);
        set(gcf, 'Color', 'w');
        set(gcf, 'PaperPositionMode', 'auto')
        
        % save results
        imagename=['response_matrix_',DaysOfRecording_V1{ii},'_b',num2str(Blocks_V1{ii}(jj)),'.jpg'];
        saveas(f,imagename, 'jpg')
        save(['SPIKEMAT_',DaysOfRecording_V1{ii},'_b',num2str(Blocks_V1{ii}(jj)),'.mat'],'SPIKEtimes','SPIKEcounts','SPIKEmean','SPIKEstd')
        
        fprintf(['\nSession ',DaysOfRecording_V1{ii},' b',num2str(Blocks_V1{ii}(jj)),' spike matrix produced!\n'])
        fprintf(['\n-------------------------------\n'])
        
    end
end

%
close all
%

%% ---------------------------------- LL ----------------------------------

for ii=1:length(DaysOfRecording_LL)
    for jj=1:length(Blocks_LL{ii})
        
        session_folder=['/zocconasphys1/chronic_inv_rec/Tanks/Giulio_Acute_Recording_', char(DaysOfRecording_LL{ii}),'/Block-', num2str(Blocks_LL{ii}(jj))];
        old=cd(session_folder);
        
        %% Carica dati bitcodes e tempi
        
        load('STIM.mat')
        % STIM contains data (digital bitcode vector) every 66 ms (2 frames at
        % 30 Hz) and corresponding bitcode onset times.
        % % figure; plot(data(1:10000));
        % % figure; plot(onset(1:10000));
        load('SPIKE.mat')
        load('my_times.mat')
        
        BITCODE=data;
        nu_onset=onset-0.033;  % to account for a shift of 1 frame in the anlg signal
        
        % borders of time window in wich to take trial spike times
        PRE_TIME=200/1000;
        POST_TIME=1200/1000;
        
        % borders of time window in wich to take trial spike count
        WINDOW_START=50/1000;
        WINDOW_STOP=900/1000;
        
        % bitcode diversi tipi di stimolo
        bcodes_noise = [1:80];
        bcodes_GR = [100:207];
        bcodes_GR_c = [250:267];
        bcodes_PL = [350:361];
        bcodes_FL = [555,556];
        
        % stimulus duration in seconds
        stimulusength=0.9;
        
        fprintf(['\nSession ',DaysOfRecording_LL{ii},' b',num2str(Blocks_LL{ii}(jj)),' data loaded\n'])
        fprintf(['\n-------------------------------\n'])
        
        %% Elimina artefatti da bitcode vector
        
        % trova ed elimina infiniti
        
        [c, d]=ind2sub(size(nu_onset), find(nu_onset==inf));
        istanzainf_on = c;
        if isempty(istanzainf_on) == 0
            nu_onset(istanzainf_on) = [];
            BITCODE(istanzainf_on) = [];
        end
        
        % trova ed elimina white screens
        [e, f]=ind2sub(size(BITCODE), find(BITCODE==1023));
        istanzaflash = e;
        if isempty(istanzaflash)==0
            nu_onset(istanzaflash) = [];
            BITCODE(istanzaflash) = [];
        end
        
        % trova ed elimina black screens
        [wi, fi]=ind2sub(size(BITCODE), find(BITCODE==0));
        istanzazero = wi;
        if isempty(istanzazero)==0
            nu_onset(istanzazero) = [];
            BITCODE(istanzazero) = [];
        end
        
        % trova ed elimina stimoli abortiti (bitcodes isolati)
        istanzasingoloflash = [];
        
        for h=[bcodes_GR,bcodes_PL]
            % trova indici di tutti gli elementi di bitcode vector con un certo valore
            [dio, poi]=ind2sub(size(BITCODE), find(BITCODE==h));
            
            for bibi=1:length(dio)
                sing_flash_idx=dio(bibi);
                % se non � l'ultimo che...
                if sing_flash_idx ~= numel(BITCODE)
                    % ... sia diverso dal precedente e dal successivo
                    if BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx+1) && BITCODE(sing_flash_idx) ~= BITCODE(sing_flash_idx-1)
                        istanzasingoloflash = [istanzasingoloflash, sing_flash_idx];
                        disp(['warning: idx # ', num2str(sing_flash_idx), ' of bcode # ', num2str(h), ' will be removed' ])
                    end
                    % se non � l'ultimo che...
                else
                    % ... sia diverso dal precedente
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
        
        %% Recupera tempi di onset e offset per ogni trial di ogni bitcode
        
        % create a cronologically ordered list of presented bitcodes
        BitCodeS4 = unique_no_sort(BITCODE);
        % exclude nois movies
        BitCodeS = BitCodeS4(ismember(BitCodeS4,[bcodes_GR,bcodes_PL]));
        
        trial_start=cell(zeros());
        trial_stop=cell(zeros());
        temp_All_Trial_NumberS=cell(zeros());
        
        for bcode=BitCodeS;
            
            % find indexes of elements of bitcode vector at the end of a trial (from the list of all bitcode onset of a particaolar condition)
            offIdx = find(diff(nu_onset(BITCODE==bcode))>0.5);
            % after an on ther is an off
            we=offIdx+1;
            % the first on is one
            onIdx=[1; we];
            
            tr_count=0;
            mycount=0;
            temp_onset_trial_idx=0;
            temp_onset_trial=0;
            
            for bibi2=1:length(onIdx)
                
                % index of current trial onset from the list of all bitcode onset of a particaolar condition)
                kaka=onIdx(bibi2);
                
                tr_count=tr_count+1;
                mycount=mycount+1;
                
                % find current bitcode indexes
                [ki, li]=ind2sub(size(BITCODE), find(BITCODE==bcode));
                
                % current trial onset index from the bitcode vector
                temp_onset_trial_idx=ki(kaka);
                temp_All_Trial_NumberS{bcode, bibi2} = temp_onset_trial_idx;
                % current trial onset
                temp_onset_trial = nu_onset(temp_onset_trial_idx);
                
                % store current trial onsets
                trial_start{bcode,tr_count} = temp_onset_trial;
                % store current trial offsets
                trial_stop{bcode,tr_count} = temp_onset_trial+stimulusength;
                
            end
        end
        
        STIM_START = [];
        STIM_STOP = [];
        
        % produce "global" (for every bitcode)onset and offset vector
        for zz= BitCodeS;
            temp_STIM_START = cell2mat(trial_start(zz,:));
            temp_STIM_STOP = cell2mat(trial_stop(zz,:));
            % concatenate successive bitcodes trial onset time vectors
            STIM_START = [temp_STIM_START, STIM_START];
            STIM_STOP = [temp_STIM_STOP, STIM_STOP];
        end
        
        % create vectors containing every trial onset and offset in cronological order
        aha = sort(STIM_START);
        baha = sort(STIM_STOP);
        
        %% Recupera le spikes per ogni neurone e trial di ogni bitcode
        
        NeuronS = 1:size(SPIKES.mychannel,2);
        BitCodeS3 = BITCODE(ismember(BITCODE,[bcodes_PL,bcodes_GR]));
        BIT_START=BitCodeS3;
        NumTrials = length(bcodes_PL);
        ChannelS = SPIKES.mychannel;
        
        nu_onsets_gr=cell(zeros());
        nu_offsets_gr=cell(zeros());
        First_Spike = cell(zeros());
        Last_Spike = cell(zeros());
        
        temp_nu_onset_all = [];
        temp_nu_offset_all = [];
        
        
        Trial_num=cell(zeros());
        SPIKEtimes=cell(zeros());
        SPIKEcounts=cell(zeros());
        basalFR=zeros(1,max(NeuronS));
        
        neurons=0 ;
        % for every neuron
        for nn=NeuronS;
            
            tic
            % get spike timestamps
            T_TIMES=SPIKES.myspikes{nn};
            TIMES=T_TIMES/1000;
            
            neurons = neurons+1;
            grcount=0;
            isicount=0;
            PSTH=cell(zeros());
            PSTH_25=cell(zeros());
            my_bits=unique(BitCodeS);
            
            % and for every bitcode
            for i=my_bits
                grcount=grcount+1;
                trialcount = 0;
                tr_count=0;
                
                % recover current bitcode trial onsets
                temp_my_trials=trial_start(i,1:end);
                temp_my_trials=cell2mat(temp_my_trials);
                my_trials = numel(temp_my_trials);
                
                % recover current bitcode trial onset index from the bitcode vector
                temp_trial_numbers = temp_All_Trial_NumberS(i,1:end);
                temp_trial_numbers=cell2mat(temp_trial_numbers);
                trial_numbers = temp_trial_numbers;
                
                % rename variables
                All_Trial_TimeS = temp_my_trials;
                All_Trial_NumberS = temp_trial_numbers;
                STIM_TI=All_Trial_TimeS;
                TRIAL_TI=All_Trial_NumberS;
                
                % for every trial of current bitcode
                for tt=1:numel(All_Trial_TimeS)
                    
                    TIME_START_1=trial_start{i,tt}-PRE_TIME;
                    TIME_END_1=trial_start{i,tt}+POST_TIME;
                    TIME_WINDOW_ON=trial_start{i,tt}+WINDOW_START;
                    TIME_WINDOW_OFF=trial_start{i,tt}+WINDOW_STOP;
                    % get spike to put in spike times matrix
                    SPIKES_TAKEN_1=TIMES(TIMES<TIME_END_1 & TIMES>TIME_START_1)-TIME_START_1;
                    % get spike to put in spike count matrix
                    SPIKES_TAKEN_2=TIMES(TIMES<TIME_WINDOW_OFF & TIMES>TIME_WINDOW_ON)-TIME_START_1;
                    
                    SPIKEtimes{i,nn,tt}=SPIKES_TAKEN_1;  %% spike timestamps per trial (RASTER)
                    SPIKEcounts{i,nn,tt}=numel(SPIKES_TAKEN_2);  %% spike counts per trial (RASTER)
                    
                    basalFR(nn)=basalFR(nn)+numel(TIMES(TIMES<TIME_START_1 & TIMES>TIME_START_1-0.5));
                    isicount=isicount+1;
                end
                
                clear q All_Trial_TimeS All_Trial_NumberS
                clear NEURON PSTH_RASTER My_Neurons_PL PsthAndRaster_PL Channel First_Spike Last_Spike TT tt tt0 First_Trial Last_Trial All_Spikes Trial_num T_num TIME_START TIME_END SPIKES_TAKEN My_Spikes PSTH Shape Mean_Firing_Rate
                
            end
            toc
            fprintf(['Neuron ',num2str(nn),' spikes processed \n'])
        end
        basalFR=basalFR./isicount*0.5;
        
        %% Plot response matrix
        
        cd(old)
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
        
        % plot figure
        tit=['LL response matrix ',DaysOfRecording_LL{ii},' b',num2str(Blocks_LL{ii}(jj))];
        tit=strrep(tit,'_','\_');
        f = figure;
        set(f,'Position',[10,10,1500,1000]);
        imagesc(cat(1,repmat(basalFR,[12,1]),SPIKEmean)); colormap('paruly'); c=colorbar; caxis([0,30]); xlabel('neuron number'); ylabel('bitcode number');
        set(get(c,'ylabel'),'String', 'mean firing rate (Hz)');        
        line([0 361],[11.5 11.5],'Color',[1,1,1])
        line([0 361],[11.5+100 11.5+100],'Color',[1,1,1])
        for kk=1:9
        line([0 361],[11.5+100+12*kk 11.5+100+12*kk],'Color',[1,1,1])
        end
        line([0 361],[11.5+350 11.5+350],'Color',[1,1,1])         
        Xticklabels = 0:5:size(SPIKEtimes,2);
        Xticks = linspace(0, size(SPIKEmean, 2), numel(Xticklabels)+1);
        set(gca,'XTick', Xticks,'XTickLabel', Xticklabels)
        Yticklabels = fliplr(cat(2,0,0:12:361));
        Yticks = 0:12:361+12;
        set(gca, 'YTick', Yticks,'YTickLabel', flipud(Yticklabels(:)));
        for kk=1:size(SPIKEtimes,2)
        line([kk kk],[0 12+360],'LineStyle','-.','Color',[0.5,0.5,0.5])
        end
        h = suptitle(tit);
        set(gca, 'Visible', 'on');
        set(h, 'Visible', 'on', 'FontSize', 15);
        set(gcf, 'Color', 'w');
        set(gcf, 'PaperPositionMode', 'auto')
        
        % save results
        imagename=['response_matrix_',DaysOfRecording_LL{ii},'_b',num2str(Blocks_LL{ii}(jj)),'.jpg'];
        saveas(f,imagename, 'jpg')
        save(['SPIKEMAT_',DaysOfRecording_LL{ii},'_b',num2str(Blocks_LL{ii}(jj)),'.mat'],'SPIKEtimes','SPIKEcounts','SPIKEmean','SPIKEstd')
        
        fprintf(['\nSession ',DaysOfRecording_LL{ii},' b',num2str(Blocks_LL{ii}(jj)),' spike matrix produced!\n'])
        fprintf(['\n-------------------------------\n'])
        
    end
end