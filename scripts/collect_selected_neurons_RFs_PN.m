function [] = collect_selected_neurons_RFs_PN()

% [] = collect_selected_neurons_RFs_PN()
% collect RFs of selected neurons into a single datastructure per area 
% and evaluate their quality computing power and local contrast
%--------------------------------------------------------------------------

%% set pars an paths

pars = set_pars_PN();

% set useful paths
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(genpath(data_folder));
cd(data_folder);

% set useful ariables
session_list=pars.listSessions(1,:);

% set stim pars
maxlag=pars.STA_depth;
stimwidth=pars.STA_width;
stimheight=pars.STA_height;
interpolation_factor=pars.interp_factor;
crop_pixel_size=pars.crop_pixel_size;

% load indexing
load('Indexing.mat')

% loop over areas
target_areas={'V1','LM','RL'};
for target_area_idx=1:length(target_areas)
    
    % load neuron selection results
    area=target_areas{target_area_idx};
    load([area,'_selected_neurons',filesep,'selected_population_',area,'.mat'])
       
    % initialize storage variables
    rsta=zeros(stimheight,stimwidth,maxlag,length(selectedsi));
    wsta=zeros(stimheight,stimwidth,maxlag,length(selectedsi));
    Zrsta=zeros(stimheight,stimwidth,maxlag,length(selectedsi));
    Zwsta=zeros(stimheight,stimwidth,maxlag,length(selectedsi));
    neuron_number=zeros(length(selectedsi),1);
    totalspikes=zeros(length(selectedsi),1);
    power=zeros(maxlag,length(selectedsi));
    bestfr_power=zeros(length(selectedsi),1);
    contrast=zeros(maxlag,length(selectedsi));
    bestfr_contrast=zeros(length(selectedsi),1);
    gostap=zeros(length(selectedsi),1);
    gostac=zeros(length(selectedsi),1);
    
    % loop over selected neurons
    for nidx=1:length(selectedsi)
        
        % get absolute neuron index
        n=selectedsi{nidx}(1);
        
        % retrive session name
        sind=M(n,1);
        ssidx = (M(:,1)==sind);
        current_block=M(n,2);
        current_session=session_list{sind};
        session_name=[current_session,'_b',num2str(current_block)];
        nind=M(n,10);

        % name of matfile to load
        matf=['STA_results_',session_name,filesep,'neuron_',num2str(nind),'_medium_results',filesep,'medium_ST_results_neuron_',num2str(nind),'.mat'];
        % load RF data
        load(matf)
        
        % add current neuren RF data to the corresponding field of the data structure
        rsta(:,:,:,nidx)=Dstafr;
        wsta(:,:,:,nidx)=Dwstafr;
        Zrsta(:,:,:,nidx)=DZstafr;
        Zwsta(:,:,:,nidx)=DZwstafr;
        neuron_number(nidx)=n;
        totalspikes(nidx)=Dtotspikes;
        
        % get RF shape params for every frame
        for ff=1:maxlag
            
            % select STA frame and interpolate it
            inputfr=Zwsta(:,:,ff,nidx);
            % set interpolation factr
            interp_factor=interpolation_factor;
            % interpolate
            outputfr = interpolate_RF_frame( inputfr, interp_factor );
            
            % compute STA power
            currfr_pow=smoothn(outputfr.^2,10000);
            power(ff,nidx)=mean(currfr_pow(:));
            
            % crop current STA frame
            [ ~, crop_ridx, crop_cidx  ] = apply_crop( outputfr , crop_pixel_size ,0 );
            cropped_frame=outputfr(crop_ridx, crop_cidx);
            % compute STA contrast
            [ ~, ~, ~, contrast(ff,nidx), ~] = get_shape_params( cropped_frame, 0 );
         
        end
        
        % select power best frames
        [maxpo,bestfr_power(nidx)]=randmax(power(:,nidx));
        gostap(nidx)=maxpo;
        % select contrast best frames
        [maxco,bestfr_contrast(nidx)]=randmax(contrast(:,nidx));
        gostac(nidx)=maxco;

        message=['\nNeuron ',num2str(nind),' (',session_name,') added to RF data structure\n'];
        fprintf(message)
        close all
        
    end
    
    save(['RFs_datastructure_',area],'rsta','wsta','Zrsta','Zwsta','neuron_number',...
        'totalspikes','power','bestfr_power','contrast','bestfr_contrast','gostap','gostac')
    
    % cd back to original folder
    message=['\n-----',area,' RFs collected -----\n'];
    fprintf(message)
    
end