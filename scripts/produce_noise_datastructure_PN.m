function [] = produce_noise_datastructure_PN( sind, bind )

% [] = produce_gaussian_noise_movies_datastructure_PN( sind, bind )
% generate movies spiking response datastructure for a given session
% and save results in the appropriate folders
%--------------------------------------------------------------------------

% set pars
pars=set_pars_PN();

% add useful folders
processed_data_folder=pars.processed_data_folder;
raw_data_folder=pars.raw_data_folder;
bcode_folder=pars.bitcode_folder;
code_folder=pars.code_folder;
cd(processed_data_folder)
addpath(bcode_folder);
addpath(code_folder);

% initialize logging
logging=[];
log_idx=1;

% retrive session name
load('indexing.mat')
session_list=pars.listSessions(1,:);
current_session = session_list{sind};
ssidx = (M(:,1)==sind);
blocks = unique(M(ssidx,2));
current_block=blocks(bind);
session_name=[current_session,'_b',num2str(current_block)];

%% load digital_bcode_list_noiseonly and spikes

load([bcode_folder,filesep,'STIM_',session_name,'.mat'])
load([processed_data_folder,filesep,'SPIKES_',session_name,'.mat'])
load([bcode_folder,filesep,'my_times_',session_name,'.mat'])

bcodes_noise = pars.stimPars.noises_bitcodes;    % load noises_bcodes
digital_bcode=data;                                    % load digital_bcode data
nu_onset=onset-0.033;                            % to account for a shift of 1 frame in the anlg signal
clear data name onset offset

% cd to data folder
data_folder=raw_data_folder;
cd(data_folder)

%% clean digital_bcode

% elimina digital_bcode_list_noiseonly whitescreen
[whitebit_id , ~]=ind2sub(size(digital_bcode), find(digital_bcode==1023));
if isempty(whitebit_id)==0
    nu_onset(whitebit_id) = [];
    digital_bcode(whitebit_id) = [];
end

% elimina digital_bcode_list_noiseonly blackscreen
[blackbit_id , ~]=ind2sub(size(digital_bcode), find(digital_bcode==0));
if isempty(blackbit_id)==0
    nu_onset(blackbit_id) = [];
    digital_bcode(blackbit_id) = [];
end

% reduce indexing matrix to this session only
M_red=M(M(:,1)==sind,:);
M_red=M_red(M_red(:,2)==current_block,:);

% compute list of presented noise bitcodes
digital_bcode_list = unique_no_sort(digital_bcode); % list of all presentend bitcodes in presentation order
digital_bcode_list_noiseonly = digital_bcode_list(ismember(digital_bcode_list,bcodes_noise));  % list of all presentend noise bitcodes in presentation order

% check and print whether all noises were played
n_noise_check=(numel(bcodes_noise)-numel(digital_bcode_list_noiseonly));
bool_all_noises_were_played=(n_noise_check==0);
if bool_all_noises_were_played
    fprintf(['\nAll noise chunks were played (',session_name,')\n'])
else
    warning_msg=['\nWARNING: ',num2str(n_noise_check),' chunks were skipped (',session_name,')\n'];
    fprintf(warning_msg)
    logging{log_idx}=warning_msg;
    log_idx=log_idx+1;
end

%% get real onset and offsets

nu_onsets_movies = cell([numel(digital_bcode_list_noiseonly),1]);
nu_offsets_movies = cell([numel(digital_bcode_list_noiseonly),1]);
bool_warning_frames_skipped=zeros(size(digital_bcode_list_noiseonly));
for current_noise_bitcode=digital_bcode_list_noiseonly
    
    % grab current bitcode onsets
    [ki , ~]=ind2sub(size(digital_bcode), find(digital_bcode==current_noise_bitcode));
    current_nu_onset = nu_onset(ki);
    
    % check and print whether frames were skipped
    num_onsets=numel(current_nu_onset);
    bool_all_frames_were_played=(num_onsets==909);
    if bool_all_frames_were_played
        fprintf(['\nAll frames of noise ',num2str(current_noise_bitcode),' were played (',session_name,')\n'])
        bool_warning_frames_skipped(current_noise_bitcode)=0;
    else
        warning_msg=['\nWARNING: during noise ',num2str(current_noise_bitcode),' some frames were skipped (',session_name,')\n'];
        fprintf(warning_msg)
        logging{log_idx}=warning_msg;
        log_idx=log_idx+1;
        bool_warning_frames_skipped(current_noise_bitcode)=1;
    end
    
    % create vector containing real frame duration (i.e. "period")
    nu_period=zeros(size(1,num_onsets));
    for wi = 1:num_onsets
        if wi ~= num_onsets
            nu_period(wi) = (current_nu_onset(wi+1)-(current_nu_onset(wi)))/2; % compute "real period"
        else
            nu_period(wi) = 0.033; % assign arbitrary period to the last one
        end
    end
    nu_offsets_movies{current_noise_bitcode} = current_nu_onset+nu_period';
    nu_onsets_movies{current_noise_bitcode} = current_nu_onset;
    
end

%% get spikes frame by frame

% get within block neuron list
neuron_list = 1:size(SPIKES.mychannel,2);

% initialize otput cell arrays
ts_by_frames=cell([numel(digital_bcode_list_noiseonly),numel(neuron_list)]);
ts_by_frames_startime=zeros(numel(digital_bcode_list_noiseonly),numel(neuron_list));
ts_by_frames_endtime=zeros(numel(digital_bcode_list_noiseonly),numel(neuron_list));
psth_by_frames=zeros([numel(digital_bcode_list_noiseonly),numel(neuron_list),2*size(nu_onsets_movies{1},1)]);
psth_by_frames_startime=zeros(numel(digital_bcode_list_noiseonly),numel(neuron_list));
psth_by_frames_endtime=zeros(numel(digital_bcode_list_noiseonly),numel(neuron_list));

for nn=neuron_list % for every neuron
    
    % load spiketimes and convert into seconds
    T_TIMES=SPIKES.myspikes{nn};
    TIMES=T_TIMES/1000;
    clear T_TIMES
    
    for i=digital_bcode_list_noiseonly % still not sorted
        
        if M_red(nn,8)==0 % process good units only
            
            curent_on = nu_onsets_movies{i}; % get current on times
            curent_off = nu_offsets_movies{i}; % get current off times
            frames_refresh=sort([curent_on;curent_off]); % merge on and offset to recover frame refresh time stamps (bitcode on and off every 2 refreshes!
            
            % get and cout spikes falling into each frame
            numberofframes = numel(frames_refresh);
            SPIKES_TAKEN_frames=cell([numberofframes,1]);
            num_SPIKES_TAKEN_frames=zeros([numberofframes,1]);
            for frnum = 1:numberofframes
                if frnum<numberofframes
                    SPIKES_TAKEN_frames{frnum} = TIMES(TIMES>frames_refresh(frnum) & TIMES<frames_refresh(frnum+1));
                elseif frnum==numberofframes
                    SPIKES_TAKEN_frames{frnum} = TIMES(TIMES>frames_refresh(frnum) & TIMES<frames_refresh(frnum)+0.033);
                end
                num_SPIKES_TAKEN_frames(frnum)=numel(SPIKES_TAKEN_frames{frnum});
            end
            % handle truncated chunks with less than 1818 frames
            if bool_warning_frames_skipped(i) && max(diff(curent_on))<0.1 && numel(curent_on)<909 % frames skipped at the end
                num_SPIKES_TAKEN_frames=[num_SPIKES_TAKEN_frames;zeros(1818-numberofframes,1)]; %#ok<AGROW>
            elseif bool_warning_frames_skipped(i) && max(diff(curent_on(1:909)))<0.1 && numel(curent_on)>909 % extra false frames
                num_SPIKES_TAKEN_frames=num_SPIKES_TAKEN_frames(1:1818);
                num_SPIKES_TAKEN_frames(end-1:end)=[0,0]; %last two frames can be wrong so put them to 0
            end
            
            % fill otput cell arrays
            ts_by_frames{i,nn}=SPIKES_TAKEN_frames;
            ts_by_frames_startime(i,nn)=frames_refresh(1);
            ts_by_frames_endtime(i,nn)=frames_refresh(end);
            psth_by_frames(i,nn,:)=num_SPIKES_TAKEN_frames;
            psth_by_frames_startime(i,nn)=frames_refresh(1);
            psth_by_frames_endtime(i,nn)=frames_refresh(end);
            
            fprintf(['\nNeuron ',num2str(nn),' out of ', num2str(max(neuron_list)),' , noise chunck number ',num2str(i),' out of ',num2str(max(digital_bcode_list_noiseonly)), ' processed\n'])      % '/',num2str(max(FirstTrialdigital_bcode_list_noiseonly))])
            
        else
            
            % fill otput cell arrays
            ts_by_frames{i,nn}=[];
            ts_by_frames_startime(i,nn)=NaN;
            ts_by_frames_endtime(i,nn)=NaN;
            psth_by_frames(i,nn,:)=NaN(2*size(nu_onsets_movies{1},1),1);
            psth_by_frames_startime(i,nn)=NaN;
            psth_by_frames_endtime(i,nn)=NaN;
            
            fprintf(['\nNeuron ',num2str(nn),' out of ', num2str(max(neuron_list)),' , noise chunck number ',num2str(i),' out of ',num2str(max(digital_bcode_list_noiseonly)), ' skipped (',num2str(M_red(nn,8)),')\n'])      % '/',num2str(max(FirstTrialdigital_bcode_list_noiseonly))])
            
        end
        
    end
    
end

% save results
header = {'noise_bitcode','within_block_neuron_number','frame_number'}; %#ok<NASGU>
save([processed_data_folder,filesep,['NOISEMAT_',session_name],'.mat'],'ts_by_frames','ts_by_frames_startime','ts_by_frames_endtime','psth_by_frames','psth_by_frames_startime','psth_by_frames_endtime','logging','header','-v7.3')
cd(processed_data_folder)

end

