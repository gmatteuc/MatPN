function [ tuning_curve, tuning_curve_error, pref_DIR, DSI, OSI, DCIr, DCIa, OCIr, OCIa, c_tuning_curve, c_tuning_curve_error, c_pref_DIR, c_DSI, c_OSI, c_DCIr, c_DCIa, c_OCIr, c_OCIa, m_tcorr, r_tuning_curve, t_tuning_matrix ] = get_direction_tuning_PN( nn, SF, TF, stimulustype, S)

%Comupute tuning curve and parameters (selectivity indices) for neuron nn
%and stimulustype of spatial frequency SF and temporal frequency TF

%--------------------------------------------------------------------------

% set plotting
boolplot=0;

% set some pars
pars = set_pars_PN();
DIR=pars.stimPars.DIR;
ntrials=pars.stimPars.ntrials;
ndirs=pars.stimPars.numDIR;

% select bitcodes for the tuning curve
selected_idx = get_indexes_PN( SF, TF, 'all', stimulustype  )';

% skip (output NaNs) parameters combinations for which ter is no bitcode
if isempty(selected_idx)
    tuning_curve=NaN(12,1);
    tuning_curve_error=NaN(12,1);
    pref_DIR=NaN;
    DSI=NaN;
    OSI=NaN;
    DCIr=NaN;
    DCIa=NaN;
    OCIr=NaN;
    OCIa=NaN;
    c_tuning_curve=NaN(12,1);
    c_tuning_curve_error=NaN(12,1);
    c_pref_DIR=NaN;
    c_DSI=NaN;
    c_OSI=NaN;
    c_DCIr=NaN;
    c_DCIa=NaN;
    c_OCIr=NaN;
    c_OCIa=NaN;
    m_tcorr=NaN;
    r_tuning_curve=NaN(12,1);
    t_tuning_matrix=NaN(12,20);
else
    
    
    %% get desired tuning curve
    
    
    % preallocate response significance variables
    sigperm_counts=zeros(ndirs,ntrials);
    muperm_counts=zeros(ndirs,ntrials);
    response_counts=zeros(ndirs,ntrials);
    spontaneous_counts=zeros(ndirs,ntrials);
    mask=zeros(ndirs,ntrials);
    sigperm_fr=zeros(ndirs);
    muperm_fr=zeros(ndirs);
    response_fr=zeros(ndirs);
    spontaneous_fr=zeros(ndirs);
    
    % --- DAVIDE CORRECTED ---
    
%     % get spike counts of current neuron and conditon
%     spike_count=S.SPIKEcounts(selected_idx,nn,:);
%     % fill empty trials with 0s
%     idxtorep=cellfun('isempty',spike_count);
%     spike_count(idxtorep)={0};
%     % handle ntrial<20 exception: if protocol is incomplete add zeros coresponding to missing trials
%     if size(spike_count,3)<20
%         % store incomplete spike count matrix
%         temp=spike_count;
%         % create complete but empty spike count matrix
%         spike_count=cell(size(spike_count,1),size(spike_count,2),20);
%         idxtofill=cellfun('isempty',spike_count);
%         spike_count(idxtofill)={0};
%         % overwrite existing trials
%         spike_count(:,:,1:size(temp,3))=temp;
%     else
%     end
%     % transform cell into a matrix
%     temp_spike_counts=squeeze(cell2mat(spike_count));
%     % create  response mak
%     response_mask=zeros(size(selected_idx,1),20);
%     for trial=1:20
%         % trials with no response above basal fr are eliminated
%         if sum(squeeze(temp_spike_counts(:,trial))>=S.basalFR(nn))
%             response_mask(:,trial)=ones(size(selected_idx))';
%         else
%             response_mask(:,trial)=zeros(size(selected_idx))';
%         end
%     end
%     temp_spike_counts(response_mask==0)=NaN;
%     % mean ignoring NaNs to compute "real" FR
%     tuning_curve=nanmean(temp_spike_counts,2);
%     tuning_curve(isnan(tuning_curve))=0;
%     % standard deviation NaNs to compute "real" FR standard deviation
%     tuning_curve_error=nanstd(temp_spike_counts,0,2);
%     tuning_curve_error(isnan(tuning_curve_error))=0;
    
    %-------------------------------------------------------
    
    % --- GIULIO CORRECTED ----
    
    
    % set permutation n
    nperm=500;
    % initialize storage variables
    tuning_curve=zeros(length(DIR),1);
    tuning_curve_error=zeros(length(DIR),1);
    % loop over directions
    for  dind=1:length(DIR)
    spp=cell(1,20);
    spc=0;
    spc_spont=0;
    outside_duration=(0+0.8)+(1.7-1); % pre- and post- counting window duration
    inside_duration=1-0.05; % counting window duration
    % loop over trials
    for trial=1:20
        if trial<=size(S.SPIKEtimes_bis( selected_idx(dind), nn, :),3) % if trial exist
            sp=S.SPIKEtimes_bis{ selected_idx(dind), nn, trial}-0.8;
        else % if trial does not exist (protocol aborted)
            sp=[];
        end
        spp{trial}=sp;
        % count spikes outside stimulus
        sp_tris_pre=sp(sp<0);
        sp_tris_post=sp(sp>1);
        count_spont=numel(sp_tris_post)+numel(sp_tris_pre);
        spc_spont=spc_spont+count_spont;
        % count spikes during stimulus
        sp_bis=sp(sp<1);
        sp_bis=sp_bis(sp_bis>0.05);
        count=numel(sp_bis);
        spc=spc+count;
        
        if count_spont==0
            % store permutation test results - null case
            sigperm_counts(dind,trial)=0;
            muperm_counts(dind,trial)=0;
            response_counts(dind,trial)=0;
            spontaneous_counts(dind,trial)=0;
            mask(dind,trial)=0;
        else
            smpls=round(count_spont*((1.7+0.8)/(0.7+0.8)));
            % trialwise permutation test for spontaneous count
              permsamp=1.7*rand(nperm,smpls)-0.8;
%             permsamp_pre=zeros(size(permsamp)); permsamp_pre(permsamp<0)=1;
%             permsamp_post=zeros(size(permsamp)); permsamp_post(permsamp>1)=1;
%             perm_count_spont=sum(permsamp_pre,2)+sum(permsamp_post,2);
            % trialwise permutation test for evoked count
            permsamp_bis1=zeros(size(permsamp)); permsamp_bis1(permsamp<1)=1;
            permsamp_bis2=zeros(size(permsamp)); permsamp_bis2(permsamp>0.05)=1;
            permsamp_bis=permsamp_bis1.*permsamp_bis2;
            perm_count=sum(permsamp_bis,2);
            % store permutation test results - ok case 
            sigperm_counts(dind,trial)=std(perm_count);
            muperm_counts(dind,trial)=mean(perm_count);
            response_counts(dind,trial)=count;
            spontaneous_counts(dind,trial)=count_spont;
            % create masks
            if abs(count-mean(perm_count))>4*std(perm_count)
                mask(dind,trial)=1;
            else
                mask(dind,trial)=0;
            end
        end
    
    end
    % in there are no spikes in the counting window
    if  isempty(spc)
        spc=0;
    end
    % in there are no spikes outside of the counting window
    if isempty(spc_spont) 
        spc_spont=0;
    end 
    
    % global (i.e. across trials averaged) result
    smpls_all=round(spc_spont*((outside_duration+inside_duration)/(outside_duration)));
    % global permutation test for spontaneous fr
    permsamp=1.7*rand(nperm,smpls_all)-0.8;
%     permsamp_pre=zeros(size(permsamp)); permsamp_pre(permsamp<0)=1;
%     permsamp_post=zeros(size(permsamp)); permsamp_post(permsamp>1)=1;
%     perm_count_spont=sum(permsamp_pre,2)+sum(permsamp_post,2);
    % global permutation test for evoked fr
    permsamp_bis1=zeros(size(permsamp)); permsamp_bis1(permsamp<1)=1;
    permsamp_bis2=zeros(size(permsamp)); permsamp_bis2(permsamp>0.05)=1;
    permsamp_bis=permsamp_bis1.*permsamp_bis2;
    perm_count=sum(permsamp_bis,2);
    % store global permutation test results
    sigperm_fr(dind)=std((perm_count/inside_duration)/20);
    muperm_fr(dind)=mean((perm_count/inside_duration)/20);
    response_fr(dind)=((spc)/inside_duration)/20;
    spontaneous_fr(dind)=((spc_spont)/outside_duration)/20;
    tuning_curve(dind)=abs(response_fr(dind)-spontaneous_fr(dind));
    
    end
   
    % PREALLOCA STORAGE VARIABLES!
    % ricalcolati zscoretot come media suitrial di zscmap e vedi chi viene meglio
    % decidi quali variabili fare uscire da function
    
    % inspect results
% %     if sum(nn==[10:11])
% %         f=figure; set(f,'Position',[10,10,1500,1000]);
% %         zscoretot=(response_fr-muperm_fr)./sigperm_fr;
% %         zscmap=(response_counts'-muperm_counts')./sigperm_counts';
% %         subplot(1,3,1); imagesc(zscmap); colorbar; subplot(1,3,2); imagesc(mask');
% %         subplot(1,3,3); plot(0:30:330,zscoretot,'-b','Linewidth',2); 
% %         title(['n',num2str(nn),' SF',num2str(SF),' TF',num2str(TF)])
% %         pause(2); axis square; 
% %     end
    
    %-------------------------------------------------------
    
    % --- CORRECTED ---
    
    % fill not played trials with zeros
    spike_count=S.SPIKEcounts(:,:,:);
    idxtorep=cellfun('isempty',spike_count);
    spike_count(idxtorep)={0};
    
    % handle ntrial<20 exception
    if size(spike_count,3)<20
        temp=spike_count;
        spike_count=cell(size(spike_count,1),size(spike_count,2),20);
        idxtofill=cellfun('isempty',spike_count);
        spike_count(idxtofill)={0};
        spike_count(:,:,1:size(temp,3))=temp;
    else
    end
    
    % couting responsive trials building "response masks" (1= responsive trial 0 otherwise)
    response_mask=zeros(size(selected_idx,1),20);
    for trial=1:20
        for i=1:length(selected_idx)
            idx=selected_idx(i);
            % if a given trial is not above basal FR count (or above a given FR) it as a non responsive trial
            if squeeze(cell2mat(spike_count(idx,nn,trial)))>=S.basalFR(nn)
                response_mask(i,trial)=1;
            else
                response_mask(i,trial)=0;
            end
        end
    end
    % substitute with NaNs spike counts of non responsive trials
    temp_spike_counts=squeeze(cell2mat(spike_count(selected_idx,nn,:)));
    t_tuning_matrix=temp_spike_counts;
    temp_spike_counts(response_mask==0)=NaN;
    % store "binarized" (counting responsive trial as 1 0 otherwise) tuning curve i.e. "r tuning curve"
    r_tuning_curve=sum(response_mask,2);
    % mean ignoring NaNs to compute "real" FR
    c_tuning_curve=nanmean(temp_spike_counts,2);
    c_tuning_curve(isnan(c_tuning_curve))=0;
    % standard deviation NaNs to compute "real" FR standard deviation
    c_tuning_curve_error=nanstd(temp_spike_counts,0,2);
    c_tuning_curve_error(isnan(c_tuning_curve_error))=0;
    
    % calculate mean across trial correlation coefficient as a reliability measure
    tc_corr_mat=corr(t_tuning_matrix);
    % do not take into account the diagonal
    tc_corr_mat(eye(size(tc_corr_mat,1))==1)=NaN;
    % compute mean of the correlation matrix ignoring NaNs
    m_tcorr=nanmean(tc_corr_mat(:));
    
    %------------------------------------------------------
    
    %% compute preferred directions
    
    % find preferred direction (normal)
    tuning_curve_n=tuning_curve+0.00001*rand(size(tuning_curve));
    pref_idx=find((tuning_curve_n)==max(tuning_curve_n));
    pref_DIR=DIR(pref_idx);
    
    %-------------------------------------------------------
    
    % find preferred direction (corrected)
    c_tuning_curve_n=c_tuning_curve+0.00001*rand(size(c_tuning_curve));
    c_pref_idx=find((c_tuning_curve_n)==max(c_tuning_curve_n));
    c_pref_DIR=DIR(c_pref_idx);
    
    %% compute direction and orientation selectivity indeces
    
    % handle indexes wrap-around (normal)
    null_idx=pref_idx+6;
    if null_idx>12
        null_idx=pref_idx-6;
    else
    end
    orth_idx1=pref_idx+3;
    if orth_idx1>12
        orth_idx1=orth_idx1-12;
    else
    end
    orth_idx2=pref_idx-3;
    if orth_idx2<=0
        orth_idx2=12-abs(orth_idx2);
    else
    end
    
    % calculate DSI (normal)
    DSI = (tuning_curve(pref_idx)-tuning_curve(null_idx))/(tuning_curve(pref_idx));
    % calculate OSI (normal)
    OSI = ((tuning_curve(pref_idx)+tuning_curve(null_idx))-(tuning_curve(orth_idx1)+tuning_curve(orth_idx2)))/((tuning_curve(pref_idx)+tuning_curve(null_idx)));
    
    %-------------------------------------------------------
    
    % handle indexes wrap-around (corrected)
    c_null_idx=c_pref_idx+6;
    if c_null_idx>12
        c_null_idx=c_pref_idx-6;
    else
    end
    c_orth_idx1=c_pref_idx+3;
    if c_orth_idx1>12
        c_orth_idx1=c_orth_idx1-12;
    else
    end
    c_orth_idx2=c_pref_idx-3;
    if c_orth_idx2<=0
        c_orth_idx2=12-abs(c_orth_idx2);
    else
    end
    
    % calculate DSI  (corrected)
    c_DSI = (c_tuning_curve(c_pref_idx)-c_tuning_curve(c_null_idx))/(c_tuning_curve(c_pref_idx));
    % calculate OSI  (corrected)
    c_OSI = ((c_tuning_curve(c_pref_idx)+c_tuning_curve(c_null_idx))-(c_tuning_curve(c_orth_idx1)+c_tuning_curve(c_orth_idx2)))/((c_tuning_curve(c_pref_idx)+c_tuning_curve(c_null_idx)));
    
    %% compute circular variance - DIR
    
    % find direction circular variance (normal)
    Img = 1i; ang = 1; anum = 0; denum = 0;
    for kkk = 1 : length(tuning_curve)
        
        anum = anum + (tuning_curve(kkk) * exp( Img * degtorad( DIR(kkk) )));
        denum = denum +  tuning_curve(kkk);
        
    end
    atemp = anum / denum;                                                       % Get temporary CirVarIdx
    DCIr = abs( atemp );                                                           % Get radius
    DCIa = radtodeg(angle( atemp ) / ang);
    
    %-------------------------------------------------------
    
    % find direction circular variance (corrected)
    Img = 1i; ang = 1; c_anum = 0; c_denum = 0;
    for kkk = 1 : length(c_tuning_curve)
        
        c_anum = c_anum + (c_tuning_curve(kkk) * exp( Img * degtorad( DIR(kkk) )));
        c_denum = c_denum +  c_tuning_curve(kkk);
        
    end
    c_atemp = c_anum / c_denum;                                                       % Get temporary CirVarIdx
    c_DCIr = abs( c_atemp );                                                           % Get radius
    c_DCIa = radtodeg(angle( c_atemp ) / ang);
    
    %% compute circular variance - ORI
    
    % find orientation circular variance (normal)
    Img = 2i;   ang = 2;  anum = 0; denum = 0;
    tuning_curve_ori=tuning_curve;
    for kkk = 1 : length(tuning_curve)/2
        
        anum = anum + (tuning_curve_ori(kkk) * exp( Img * degtorad( DIR(kkk) )));
        denum = denum +  tuning_curve_ori(kkk);
        
    end
    atemp = anum / denum;                                                       % Get temporary CirVarIdx
    OCIr = abs( atemp );                                                           % Get radius
    OCIa = radtodeg(angle( atemp ) / ang);                                       % Get angle in radians
    
    %-------------------------------------------------------
    
    % find orientation circular variance (corrected)
    Img = 2i;   ang = 2;  anum = 0; denum = 0;
    c_tuning_curve_ori=c_tuning_curve;
    for kkk = 1 : length(c_tuning_curve)/2
        
        c_anum = c_anum + (c_tuning_curve_ori(kkk) * exp( Img * degtorad( DIR(kkk) )));
        c_denum = c_denum +  c_tuning_curve_ori(kkk);
        
    end
    c_atemp = c_anum / c_denum;                                                       % Get temporary CirVarIdx
    c_OCIr = abs( c_atemp );                                                           % Get radius
    c_OCIa = radtodeg(angle( c_atemp ) / ang);                                       % Get angle in radians
    
    %% visualize and plot some exemples
    
    if boolplot==1
        
        switch stimulustype
            case 'grating'
                
                % plot only for good neurons
                if nn<=numel(S.goodc)
                    
                    fname=strrep(['Neuron_',num2str(nn),'_SF',num2str(SF),'_TF',num2str(TF),'_grating'],'.','');
                    if exist([fname,'.jpg'])~=2
                        
                        % always re-use the same figure to avoid popping up
                        if not(isempty(get(0,'children')))
                            f1=get(0,'children');
                        else
                            f1 = figure;
                        end
                        set(f1,'Position',[10,10,1500,1000]);
                        
                        % plot tuning curves
                        sb1=subplot(666,666,666);
                        par=parula;
                        coffset=10;
                        plot(tuning_curve,'Color',par(end-coffset,:),'LineWidth',3.5); hold on; plot(c_tuning_curve,'Color',par(coffset,:),'LineWidth',3.5); plot(r_tuning_curve./20,'k','LineWidth',2.5); plot(ones(size(r_tuning_curve)),'--k');
                        xlim([1,12]); ylim([0,20]); legend('original tc','corrected tc','% responsive tc'); ylabel('firing rate (Hz)'); xlabel('direction (deg)');
                        set(gca,'XTickLabel',num2str(round(linspace(0,330,12)')));
                        title(['Neuron ',num2str(nn),' SF=',num2str(SF),' TF=',num2str(TF),' - grating tuning curves']);
                        text(10.4,15,['DSI=',num2str(DSI,'%.2f')],'FontSize',10);
                        text(10.4,16,['OSI=',num2str(OSI,'%.2f')],'FontSize',10);
                        text(9,15,['cDSI=',num2str(c_DSI,'%.2f')],'FontSize',10);
                        text(9,16,['cOSI=',num2str(c_OSI,'%.2f')],'FontSize',10);
                        text(9,14,['mtcorr=',num2str(m_tcorr,'%.2f')],'FontSize',10);
                        set(sb1,'Position',[.05,.13,.41,.75]);
                        
                        % plot trial/direction spike count matrix ("t_tuning_matrix")
                        sb2=subplot(666,666,666);
                        imagesc(flipud(t_tuning_matrix')); colormap('parula');
                        ylabel('trial number'); xlabel('direction (deg)');
                        set(gca,'XTickLabel',num2str([30;90;150;210;270;330]));
                        set(gca,'YTickLabel',num2str(flipud((2:2:20)')-1));
                        hc=colorbar; ylabel(hc,'spike count');
                        title(['spike count across trials and directions']);
                        set(sb2,'Position',[.51,.13,.41,.75]);
                        
                        % save figure
                        set(gcf, 'PaperPositionMode', 'auto')
                        saveas(gcf,fname,'jpeg')
                        clf(f1)
                        
                    else
                    end
                    
                else
                end
            case 'plaid'
                
                % plot only for good neurons
                if nn<=numel(S.goodc)
                    
                    fname=strrep(['Neuron_',num2str(nn),'_SF',num2str(SF),'_TF',num2str(TF),'_plaid'],'.','');
                    if exist([fname,'.jpg'])~=2
                        
                        % always re-use the same figure to avoid popping up
                        if not(isempty(get(0,'children')))
                            f1=get(0,'children');
                        else
                            f1 = figure;
                        end
                        set(f1,'Position',[10,10,1500,1000]);
                        
                        % plot tuning curves
                        sb1=subplot(666,666,666);
                        par=winter;
                        coffset=10;
                        plot(tuning_curve,'Color',par(end-coffset,:),'LineWidth',3.5); hold on; plot(c_tuning_curve,'Color',par(coffset,:),'LineWidth',3.5); plot(r_tuning_curve./20,'k','LineWidth',2.5); plot(ones(size(r_tuning_curve)),'--k');
                        xlim([1,12]); ylim([0,20]); legend('original tc','corrected tc','% responsive tc'); ylabel('firing rate (Hz)'); xlabel('direction (deg)');
                        set(gca,'XTickLabel',num2str(round(linspace(0,330,12)')));
                        title(['Neuron ',num2str(nn),' SF=',num2str(SF),' TF=',num2str(TF),' - plaid tuning curves']);
                        text(10.4,15,['DSI=',num2str(DSI,'%.2f')],'FontSize',10);
                        text(10.4,16,['OSI=',num2str(OSI,'%.2f')],'FontSize',10);
                        text(9,15,['cDSI=',num2str(c_DSI,'%.2f')],'FontSize',10);
                        text(9,16,['cOSI=',num2str(c_OSI,'%.2f')],'FontSize',10);
                        text(9,14,['mtcorr=',num2str(m_tcorr,'%.2f')],'FontSize',10);
                        set(sb1,'Position',[.05,.13,.41,.75]);
                        
                        % plot trial/direction spike count matrix ("t_tuning_matrix")
                        sb2=subplot(666,666,666);
                        imagesc(flipud(t_tuning_matrix')); colormap('winter');
                        ylabel('trial number'); xlabel('direction (deg)');
                        set(gca,'XTickLabel',num2str([30;90;150;210;270;330]));
                        set(gca,'YTickLabel',num2str(flipud((2:2:20)')-1));
                        hc=colorbar; ylabel(hc,'spike count');
                        title(['spike count across trials and directions']);
                        set(sb2,'Position',[.51,.13,.41,.75]);
                        
                        % save figure
                        set(gcf, 'PaperPositionMode', 'auto')
                        saveas(gcf,fname,'jpeg')
                        clf(f1)
                        
                    else
                    end
                    
                else
                end
                
        end
    end
    
end
end