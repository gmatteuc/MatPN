function [ tuning_curve, c_tuning_curve, tuning_curve_matrix, pref_DIR,...
    OSI, DSI, DI, c_pref_DIR, c_OSI, c_DSI, c_DI, ...
    sigperm_counts, muperm_counts, response_counts,...
    spontaneous_counts, sigperm_fr, muperm_fr, response_fr, spontaneous_fr,...
    tcorr, sigtrials ] = get_direction_tuning_PN_NEW( nn, SF, TF, stimulustype, S)

% [ tuning_curve, c_tuning_curve, tuning_curve_matrix, pref_DIR,...
%     OSI, DSI, DI, c_pref_DIR, c_OSI, c_DSI, c_DI, ...
%     sigperm_counts, muperm_counts, response_counts,...
%     spontaneous_counts, sigperm_fr, muperm_fr, response_fr, spontaneous_fr,...
%     tcorr, sigtrials ] = get_direction_tuning_PN_NEW( nn, SF, TF, stimulustype, S)
%
%Comupute tuning curve and parameters (selectivity indices) for neuron nn
%and stimulustype of spatial frequency SF and temporal frequency TF

%--------------------------------------------------------------------------

% set some pars
pars = set_pars_PN();
DIR=pars.stimPars.DIR;
ntrials=pars.stimPars.ntrials;
ndirs=pars.stimPars.numDIR;
pre_time=pars.stimPars.pre_time;
post_time=pars.stimPars.post_time;
stim_time=pars.stimPars.stim_time;
pre_delay=pars.stimPars.pre_delay;
post_delay=pars.stimPars.post_delay;

% select bitcodes for the tuning curve
selected_idx = get_indexes_PN( SF, TF, 'all', stimulustype  )';

% skip (output NaNs) parameters combinations for which ter is no bitcode
if isempty(selected_idx)
    
    tuning_curve=NaN(12,1);
    c_tuning_curve=NaN(12,1);
    tuning_curve_matrix=NaN(ndirs,ntrials);
    pref_DIR=NaN;
    OSI=NaN;
    DSI=NaN;
    DI=NaN;
    c_pref_DIR=NaN;
    c_OSI=NaN;
    c_DSI=NaN;
    c_DI=NaN;
    sigperm_counts=NaN(ndirs,ntrials);
    muperm_counts=NaN(ndirs,ntrials);
    response_counts=NaN(ndirs,ntrials);
    spontaneous_counts=NaN(ndirs,ntrials);
    sigperm_fr=NaN(ndirs,1);
    muperm_fr=NaN(ndirs,1);
    response_fr=NaN(ndirs,1);
    spontaneous_fr=NaN(ndirs,1);
    tcorr=NaN;
    
    message='WARNING: skipped bitcode!!!';
    fprintf(message);
    
else
    
    
    %% get desired tuning curve - "raw"
        
    % preallocate response significance variables
    sigperm_counts=zeros(ndirs,ntrials);
    muperm_counts=zeros(ndirs,ntrials);
    response_counts=zeros(ndirs,ntrials);
    spontaneous_counts=zeros(ndirs,ntrials);
    mask=zeros(ndirs,ntrials);
    sigperm_fr=zeros(ndirs,1);
    muperm_fr=zeros(ndirs,1);
    response_fr=zeros(ndirs,1);
    spontaneous_fr=zeros(ndirs,1);
    
    % set permutation n
    nperm=500;
    masktreshold=3; % in sigma units
    
    % loop over directions
    for  dind=1:length(DIR)
     
    % initialize spike counters  
    spc=0;
    spc_spont=0;
    
    outside_duration=(pre_time)+(post_time-post_delay); % pre- and post- counting window duration
    inside_duration=stim_time-pre_delay+post_delay; % counting window duration
    
    % loop over trials
    for trial=1:ntrials
        
        if trial<=size(S.SPIKEtimes_bis( selected_idx(dind), nn, :),3) % if trial exist
            sp=S.SPIKEtimes_bis{ selected_idx(dind), nn, trial}-pre_time;
        else % if trial does not exist (protocol aborted)
            sp=[];
        end
        
        % count spikes outside stimulus
        sp_tris_pre=sp(sp<0);
        sp_tris_post=sp(sp>(stim_time+post_delay));
        count_spont=numel(sp_tris_post)+numel(sp_tris_pre);
        spc_spont=spc_spont+count_spont;
        
        % count spikes during stimulus
        sp_bis=sp(sp<(stim_time+post_delay));
        sp_bis=sp_bis(sp_bis>pre_delay);
        count=numel(sp_bis);
        spc=spc+count;
                
        % -------------------- trialwise permutation ----------------------
        
        if count_spont==0
            
            % store permutation test results - null case
            sigperm_counts(dind,trial)=eps;
            muperm_counts(dind,trial)=eps;
            response_counts(dind,trial)=count;
            spontaneous_counts(dind,trial)=0;
            mask(dind,trial)=0;
            
        else
            
            % set permutation sample num = number of expected spontaneous spikes
            smpls=round(count_spont*((outside_duration+inside_duration)/(outside_duration)));
            % generate permutation samples
            permsamp=(stim_time+post_time)*rand(nperm,smpls)-pre_time;
            % trialwise permutation test for evoked count
            permsamp_bis1=zeros(size(permsamp)); permsamp_bis1(permsamp<(stim_time+post_delay))=1;
            permsamp_bis2=zeros(size(permsamp)); permsamp_bis2(permsamp>pre_delay)=1;
            permsamp_bis=permsamp_bis1.*permsamp_bis2;
            perm_count=sum(permsamp_bis,2);
            
            % store permutation test results - ok case 
            sigperm_counts(dind,trial)=std(perm_count); % permutation estimated stds
            muperm_counts(dind,trial)=mean(perm_count); % permutation estimated means
            response_counts(dind,trial)=count;          % observed responses
            spontaneous_counts(dind,trial)=round(count_spont*...
                ((inside_duration)/(outside_duration))); % observed background
            if (count-mean(perm_count))>masktreshold*std(perm_count)
                mask(dind,trial)=1;
            end
            
            % -------------------------------------------------------------
            
        end
    
    end

    % -------------------- across trials permutation ----------------------
    
    % global (i.e. across trials averaged) result
    smpls_all=round(spc_spont*((outside_duration+inside_duration)/(outside_duration)));
    % generate permutation samples
    permsamp=(stim_time+post_time)*rand(nperm,smpls_all)-pre_time;
    % global permutation test for evoked fr
    permsamp_bis1=zeros(size(permsamp)); permsamp_bis1(permsamp<(stim_time+post_delay))=1;
    permsamp_bis2=zeros(size(permsamp)); permsamp_bis2(permsamp>pre_delay)=1;
    permsamp_bis=permsamp_bis1.*permsamp_bis2;
    perm_count=sum(permsamp_bis,2);
    % store global permutation test results
    sigperm_fr(dind)=std((perm_count/inside_duration)/ntrials);              % permutation estimated stds
    muperm_fr(dind)=mean((perm_count/inside_duration)/ntrials);              % permutation estimated means
    response_fr(dind)=((spc)/inside_duration)/ntrials;                       % observed responses
    spontaneous_fr(dind)=((spc_spont)/outside_duration)/ntrials;             % observed background
    
    % -------------------------------------------------------------
     
    end
    
    % compute background corrected tc
    tuning_curve=max(0,response_fr-spontaneous_fr);
    % compute background corrected tc matrix
    tuning_curve_matrix=max(0,response_counts-spontaneous_counts)./inside_duration;
  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % ------------------------- to inspect results ------------------------- 
% %     if sum(nn==[3:6])
% %         f=figure; set(f,'Position',[10,10,1500,1000]);
% %         zscoretot=(response_fr-muperm_fr)./sigperm_fr;
% %         zscmap=(response_counts'-muperm_counts')./sigperm_counts';
% %         subplot(1,3,1); imagesc(zscmap); axis square; colorbar; subplot(1,3,2); imagesc(mask');
% %         colorbar; axis square;
% %         subplot(1,3,3); plot(0:30:330,zscoretot,'-b','Linewidth',2); hold on;
% %         plot(0:30:330,nanmean(zscmap,1),'Color',[0,0.7,0],'Linewidth',2); axis square;
% %         plot([0,330],[3,3],'--k','Linewidth',1); plot([0,330],[-3,-3],'--k','Linewidth',1);
% %         title(['n',num2str(nn),' - SF',num2str(SF),' - TF',num2str(TF)])
% %         pause(2);
% %     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% get desired tuning curve - "corrected"
    
    % select responsive trials (at least one signif. condition)
    sign_trial_mask=NaN(ndirs,ntrials);
    % significant tral vector (boolean)
    sigtrials=(sum(mask,1)~=0);
    sign_trial_mask(:,sigtrials)=ones(ndirs,sum(sigtrials));
    % recompute tuning curve omitting  nonresponsive trials
    c_tuning_curve=(nansum( tuning_curve_matrix.*sign_trial_mask,2 )...
        ./ntrials)./inside_duration;  

    %% compute reliability metric
    
    % calculate mean across trial correlation coefficient as a reliability measure
    tc_corr_mat=corr(tuning_curve_matrix);
    % do not take into account the diagonal
    tc_corr_mat(eye(size(tc_corr_mat,1))==1)=NaN;
    % compute mean of the correlation matrix ignoring NaNs
    tcorr=nanmean(tc_corr_mat(:));
    
    
    %% compute direction and orientation selectivity indeces
    
    % on "raw" tuning curves
    [ OSI,DSI,DI,pref_DIR,~  ] = compute_SIs( tuning_curve );
    
    % on "raw" tuning curves
    [ c_OSI,c_DSI,c_DI,c_pref_DIR,~  ] = compute_SIs( c_tuning_curve );
       
    
end
end