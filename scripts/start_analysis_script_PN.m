% ------------------------- START ANALYSIS SCRIPT ------------------------- 

clear 
close all
clc

%% set analysis parameters

startingsession=1; % modify starting session to easily add new one
overwrite=0;
overwrite_indexing=0;
analysis_path='C:\Users\labuser\Desktop\code_development_scratchfolder\PN_analysis';
cd(analysis_path);

%% set pars

addpath(genpath(analysis_path));
pars=set_pars_PN();
session_list=pars.listSessions(1,:);
home=pars.processed_data_folder;
cd(home)

%% core loop

fprintf('         ...STARTING DATA EXTRACTION...\n')
fprintf('\n------------------------------------------------------------\n')

%%  do spike triggered analysis

numSS=length(session_list);
for sind=startingsession:numSS
    
    session = session_list{sind};
    blocks=pars.listSessions{2,sind};
    numBL=length(blocks);
    
    for bind=1:numBL
        
        current_block=blocks(bind);
        sname=[session,'_b',num2str(current_block)];
        if exist(['SPIKEMAT_',sname,'.mat'],'file')==2 && exist(['SPIKES_',sname,'.mat'],'file')==2  && overwrite==0
            fprintf(['\n',sname,' data already extracted, skipping... \n'])
        else
            produce_spikemats_PN(sname);
            fprintf(['\n',sname,' data extraction completed... \n'])
        end
        fprintf('------------------------------------------------------------\n\n')
        
    end
    
end

%% do indexing

fprintf('------------------------------------------------------------\n')
if exist('Indexing.mat','file')==2 && overwrite_indexing==0
    fprintf('data already indexed, skipping... \n\n')
else
    fprintf('starting indexing... \n')
    tic
    do_indexing_PN;
    toc
    fprintf('indexing done... \n')
end
fprintf('------------------------------------------------------------\n\n')

fprintf(' ...STARTING NOISE RESPONSES PREPROCESSING ...\n\n')
fprintf('------------------------------------------------------------\n')

%%  do spike triggered analysis

numSS=length(session_list);

for sind=startingsession:numSS
    
    session = session_list{sind};
    blocks=pars.listSessions{2,sind};
    numBL=length(blocks);
    
    for bind=1:numBL
        
        current_block=blocks(bind);
        sname=[session,'_b',num2str(current_block)];
        
        if exist(['NOISEMAT_',sname,'.mat'],'file')==2 && overwrite==0
            fprintf([sname,' noise response datastructure already exist, skipping... \n'])
            fprintf('------------------------------------------------------------\n')
        else
            produce_noise_datastructure_PN( sind, bind );
            fprintf('\n------------------------------------------------------------\n')
            fprintf([sname,' noise response datastructure created... \n'])
        end
        fprintf('------------------------------------------------------------\n')
        spike_triggered_average_gaussian_PN_MEDIUM( sind, bind )
        fprintf([sname,' spike triggered anlysis performed... \n'])
        fprintf('------------------------------------------------------------\n')
        
    end
    
end

%%  analyze tuning

% do_tuning_analysis_PN

%%  select neurons to be included into the analysis

% do_selectivity_analysis_PN