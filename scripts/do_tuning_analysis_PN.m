% ------------------------- TUNING ANALYSIS SCRIPT -------------------------
close all
clear
clc

%% set paths

pars = set_pars_PN();
data_folder=pars.processed_data_folder;
addpath(data_folder);
code_folder=pars.code_folder;
addpath(code_folder);
oldfold=cd(data_folder);

%% set prams

SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
numSF=pars.stimPars.numSF;
numTF=pars.stimPars.numTF;
numDIR=pars.stimPars.numDIR;
stimulustype=pars.stimPars.stimulustype;
numstimulustype=pars.stimPars.numstimulustype;
gratings_bitcodes=pars.stimPars.gratings_bitcodes;
plaids_bitcodes=pars.stimPars.plaids_bitcodes;
nreduxbit=numel(gratings_bitcodes)+numel(plaids_bitcodes);
ntrials=pars.stimPars.ntrials;
winwidth=pars.stimPars.winwidth;


%% load indexing and set session list

load('Indexing.mat')
session_list=pars.listSessions(1,:);

%% core loop

numSS=numel(unique(M(:,1)));
% loop over sessions
for ss=1:numSS
    
    % retrive current session date and block number
    session = session_list{ss};
    % current session indexes
    ssidx=find(M(:,1)==ss);
    % retrive their blocktype
    blocks = unique(M(ssidx,2));
    
    numBL=length(blocks);
    % loop over blocks
    for bb=1:numBL
        % find first and last index of current block
        blidx=intersect(find(M(:,2)==blocks(bb)),ssidx);
        
        % compile current block and session SPIKEMAT filename
        folder = fullfile(data_folder);
        filename=fullfile(folder, ['SPIKEMAT_',session,'_b',num2str(blocks(bb) ),'.mat']);
        
        % load SPIKEMAT
        tic
        fprintf( [ 'Loading session ', session, ' block ', num2str(blocks(bb) ), '\n' ])
        S=load(filename);
        fprintf('Done.\n')
        toc
        
        % create ouput folder
        dirnam=['tuning_results_',session,'_b',num2str(blocks(bb))];
        mkdir(dirnam);
        oldd=cd(dirnam);
        dirnam_bis=['tuning_curves'];
        mkdir(dirnam_bis);
        cd(dirnam_bis);
        
        % perform tuning analysis on current session
        [ tuning_curve(:,blidx,:,:,:), c_tuning_curve(:,blidx,:,:,:), d_tuning_curve(:,blidx,:,:,:),...
            tuning_curve_matrix(:,:,blidx,:,:,:), pref_DIR(blidx,:,:,:),...
            OSI(blidx,:,:,:), DSI(blidx,:,:,:), DI(blidx,:,:,:), c_pref_DIR(blidx,:,:,:),...
            c_OSI(blidx,:,:,:), c_DSI(blidx,:,:,:), c_DI(blidx,:,:,:), ...
            d_pref_DIR(blidx,:,:,:), d_OSI(blidx,:,:,:), d_DSI(blidx,:,:,:), d_DI(blidx,:,:,:), ...
            sigperm_counts(:,:,blidx,:,:,:), muperm_counts(:,:,blidx,:,:,:),...
            response_counts(:,:,blidx,:,:,:), spontaneous_counts(:,:,blidx,:,:,:),...
            sigperm_fr(:,blidx,:,:,:), muperm_fr(:,blidx,:,:,:), response_fr(:,blidx,:,:,:),...
            spontaneous_fr(:,blidx,:,:,:), tcorr(blidx,:,:,:), sigtrials(:,blidx,:,:,:), ...
            tuning_matrix(:,:,blidx,:), tuning_matrix_error(:,:,blidx,:), pref_SF(blidx,:),...
            pref_TF(blidx,:), c_tuning_matrix(:,:,blidx,:), c_tuning_matrix_error(:,:,blidx,:),...
            c_pref_SF(blidx,:), c_pref_TF(blidx,:), ...
            d_tuning_matrix(:,:,blidx,:), d_tuning_matrix_error(:,:,blidx,:),...
            d_pref_SF(blidx,:), d_pref_TF(blidx,:), ...
            Zp(blidx,:,:), Zc(blidx,:,:), Rp(blidx,:,:), Rc(blidx,:,:), PI(blidx,:,:), ...
            c_Zp(blidx,:,:), c_Zc(blidx,:,:), c_Rp(blidx,:,:), c_Rc(blidx,:,:), c_PI(blidx,:,:), ...
            d_Zp(blidx,:,:), d_Zc(blidx,:,:), d_Rp(blidx,:,:), d_Rc(blidx,:,:), d_PI(blidx,:,:), ...
            modulation_index(:,blidx,:,:,:), modulation_index_bis(:,blidx,:,:,:), ...
            tuning_curve_z(:,blidx,:,:,:), basal_fr(blidx), basal_fr_std(blidx)] ...
            = analyze_tuning_PN( S ); %#ok<*SAGROW>
        
        cd(oldd)
        
        fprintf( [ '\nTuning session ', session, ' block ', num2str(blocks(bb) ), ' analyzed\n\n' ])
        
    end
end


%% save results

save('Tuning.mat',...
    'tuning_curve', 'c_tuning_curve', 'd_tuning_curve',...
    'tuning_curve_z','basal_fr','basal_fr_std',...
    'tuning_curve_matrix', 'pref_DIR',...
    'OSI', 'DSI', 'DI', 'c_pref_DIR',...
    'c_OSI', 'c_DSI', 'c_DI', 'd_pref_DIR',...
    'd_OSI', 'd_DSI', 'd_DI', ...
    'sigperm_counts', 'muperm_counts',...
    'response_counts', 'spontaneous_counts',...
    'sigperm_fr', 'muperm_fr', 'response_fr',...
    'spontaneous_fr', 'tcorr', 'sigtrials', 'tuning_matrix',...
    'tuning_matrix_error', 'pref_SF', 'pref_TF',...
    'c_tuning_matrix', 'c_tuning_matrix_error', 'c_pref_SF', 'c_pref_TF',...
    'd_tuning_matrix', 'd_tuning_matrix_error', 'd_pref_SF', 'd_pref_TF',...
    'Zp', 'Zc', 'Rp', 'Rc', 'PI', ...
    'c_Zp', 'c_Zc', 'c_Rp', 'c_Rc', 'c_PI', ...
    'd_Zp', 'd_Zc', 'd_Rp', 'd_Rc', 'd_PI', ...
    'modulation_index', 'modulation_index_bis')

cd(oldfold)

%% plot raster mosaics

% for h=1:2
%     for hh=1:2
%         for hhh=2
%             plot_raster_mosaics( h,hh,stimulustype{hhh} )
%         end
%     end
% end