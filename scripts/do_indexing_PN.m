function M = do_indexing_PN()

% M = do_indexing_PN()
% assign univoque integer index to all neurons in a session/block;
% this should be done after each modification of the dataset such
% as addition or removal of sessions. This require a careful check
% in set_pars when defining listSessions
%--------------------------------------------------------------------------

% set parameters
pars = set_pars_PN();
% reassign important variables from the pars structure for ease of use
listSessions = pars.listSessions;
dataFolder = pars.processed_data_folder;

% initialize neuron index
index = 0;
% initialize indexing matrix
M = [];

for i = 1:size(listSessions,2)
    
    % get current session parameters
    session                = listSessions{1,i};
    blocks                 = listSessions{2,i};
    blocktypes             = listSessions{3,i};
    area                   = listSessions{4,i};
    recordingtype          = listSessions{5,i};
    depth                  = listSessions{6,i};
    
    maincounter=0;
    
    for j = 1: numel(blocks)
        
        if blocktypes(j)==1
            maincounter=maincounter+1;
            manual_area=listSessions{7+(maincounter-1),i};
        else
            manual_area=99*ones(1,100);
        end
        
        % compile current block data folder
        folder = fullfile(dataFolder);
        filename=fullfile(folder, ['SPIKEMAT_',session,'_b',num2str(blocks(j) ),'.mat']);
        
        if exist(filename,'file')==2
            
            % load data
            tic
            fprintf( [ 'Loading session ', session, ' block ', num2str(blocks(j) ), '\n' ])
            S = load(filename);
            fprintf('Done.\n')
            toc
            
            % get different cluster types relative indeces
            ctypelabels=cat(2,zeros(1,length(S.goodc)),1*ones(1,length(S.muac)));
            ctypelabels=cat(2,ctypelabels,2*ones(1,length(S.noisec)));
            clabels=[S.goodc;S.muac;S.noisec];
            
            for k = 1: size(S.SPIKEmean,2)
                
                % increment neuron index
                index = index + 1;
                
                M(index, 1) = i; % session (index)
                M(index, 2) = blocks(j); % block
                M(index, 3) = blocktypes(j); % blocktype
                M(index, 4) = area(j); % targeted area
                M(index, 5) = recordingtype; % spherical correction (0 corrected, 1 not corrected)
                M(index, 6) = depth; % depth (0 deep, 1 superficial)
                M(index, 7) = manual_area(k); % manually annotated area
                M(index, 8) = ctypelabels(k); % cluster label (within session as assigned by Klusta)
                M(index, 9) = clabels(k); % clustertype label (0=good 1=mua 2=noise) (WARNING: 8 e 9 scambiati? 27/11/2017)
                M(index, 10) = k; % within session index

            end
            
        else
            
        end
        
    end
    
end

% save indexing matrix
save('Indexing.mat','M')
end

