function pars = set_pars_PN()

% pars = set_pars_PN()
% set analysis paths and parameters along with stimuli and sessions metadata
% -------------------------------------------------------------------------

%% initialize pars structure

pars = struct();

%% set paths

current_dir=pwd;
if not(isunix)
current_dir=split(current_dir,'\PN_analysis');
current_dir=current_dir{1};
pars.raw_data_folder = [current_dir,'\PN_analysis\raw_data'];
pars.processed_data_folder = [current_dir,'\PN_analysis\processed_data'];
pars.code_folder = genpath([current_dir,'\PN_analysis\functions']);
pars.stim_folder = [current_dir,'\PN_analysis\stims'];
pars.bitcode_folder = [current_dir,'\PN_analysis\bitcodes'];
else
current_dir=strrep(current_dir,'/PN_analysis/scripts','');
pars.raw_data_folder = [current_dir,'/PN_analysis/raw_data'];
pars.processed_data_folder = [current_dir,'/PN_analysis/processed_data'];
pars.code_folder = genpath([current_dir,'/PN_analysis/functions']);
pars.stim_folder = [current_dir,'/PN_analysis/stims'];
pars.bitcode_folder = [current_dir,'/PN_analysis/bitcodes'];
end

%% list of sessions/blocks and their type

templist={'02_03_2017','04_09_2017','14_09_2017','18_09_2017'...
};

listSessions = cell(9,size(templist,2));
listSessions(1,:) = templist;

listSessions(2,:) = {[2],[10],[7],[8] ...
};                              % block numbers
listSessions(3,:) = {[1],[1],[1],[1] ...
};                              % blocktype (1=main protocol, 0=RF mapping)
listSessions(4,:) = {0,1,1,1 ...
};                                % targeted area (0=V1 fentanest, 1=LM fentanest, 2=RL fentanest, 3=AL fentanest)
listSessions(5,:) = {0,0,0,0 ...
};                                % recording type (1=head-fixed,0=acute)
listSessions(6,:) = {1,1,1,1 ...
};                                % targeted depth (1=horizontal penetration,0=vertical penetration)

listSessions(7,:) = {[zeros(1,200),zeros(1,200)],...
    [ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)]...
};      % area manual annotation - before first main block - (0=V1 fentanest, 1=LM fentanest, 2=RL fentanest, 3=AL fentanest)
listSessions(8,:) = {[zeros(1,200),zeros(1,200)],...
    [ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)]...
};      % area manual annotation - before second main block - (0=V1 fentanest, 1=LM fentanest, 2=RL fentanest, 3=AL fentanest)
listSessions(9,:) = {[zeros(1,200),zeros(1,200)],...
    [ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)],[ones(1,200),ones(1,200)]...
};      % area manual annotation - before third main block - (0=V1 fentanest, 1=LM fentanest, 2=RL fentanest, 3=AL fentanest)
pars.listSessions = listSessions;

%% stimulus parameters

pars.stimPars = [];

pars.stimPars.SF=[0.02,0.04];                        % (Spatial Frequencies)
pars.stimPars.TF=[2,6];                              % (Temporal Frequencies)
pars.stimPars.DIR=0:30:330;                          % (DIRections)
pars.stimPars.numSF=length(pars.stimPars.SF);
pars.stimPars.numTF=length(pars.stimPars.TF);
pars.stimPars.numDIR=length(pars.stimPars.DIR);

pars.stimPars.stimulustype={'grating','plaid'};
pars.stimPars.numstimulustype=length(pars.stimPars.stimulustype);

pars.stimPars.gratings_bitcodes=[100:111,112:123,124:135,136:147];   % (first increasing SF, second increasing TF, DIR increasing within chunk)
pars.stimPars.plaids_bitcodes=[350:361,362:373,374:385,386:397];     % (first increasing SF, second increasing TF, DIR increasing within chunk)
pars.stimPars.flashes_bitcodes=[554,555];                            % (first black, second white)
pars.stimPars.noises_bitcodes=[1:40];                                % (30 s = 909 frames chunks of gaussianly correlated 0.04 cpd-matched 5s contrast modulated noise)

pars.stimPars.ntrials=20;
pars.stimPars.winwidth=0.85;
pars.stimPars.pre_time=0.8;
pars.stimPars.post_time=0.8;
pars.stimPars.stim_time=0.9;
pars.stimPars.pre_delay=0.05;
pars.stimPars.post_delay=0.1;

pars.neuPars.ctypes={'good','mua','noise'};
pars.neuPars.clabels={2,1,0};

end