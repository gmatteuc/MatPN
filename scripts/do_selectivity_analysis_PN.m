% ------------------------- DO NEURON SELECTION -------------------------

clear
close all
clc

% set pars
pars = set_pars_PN();
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
types=pars.neuPars.ctypes;
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);
cd(data_folder);

% loop over areas
target_areas={'V1','LM','RL'};

for target_area_idx=1:2%:length(target_areas)
    
    % set target area
    target_group=target_areas{target_area_idx};
    outfold=[target_group,'_selected_neurons'];
    outfold=strrep(outfold,'','_');
    mkdir(outfold);
    oldfold=cd(outfold);
    
    % load tuning analysis results
    load('Tuning.mat')
    load('Indexing.mat')
%     trial_matrix=response_counts-muperm_counts./sigperm_counts;
    
    %%  -------------------------- neuron selection --------------------------
    
    % get good single units
    goodidx_su=find(M(:,8)==0);
    
    % get target group index
    switch target_group
        case 'V1'
            goodidx_ar=find(M(:,4)==0);
        case 'LM'
            goodidx_ar=find(M(:,4)==1);
%         case 'AL'
%             goodidx_ar=find(M(:,4)==2);
        case 'RL'
            goodidx_ar=find(M(:,4)==2);
    end
    % intersect to get good neuron index
    goodidxtot=intersect(goodidx_su,goodidx_ar,'stable');
    
    % selecting neurons ...
    
    % initialize selection storage variables
    neuroncounter=0;
    selectedsi=cell(0);
    gooosi=[];
    gooocia=[];
    goodsi=[];
    goodsf=[];
    goodtf=[];
    gooddir=[];
    gooddir_p=[];
    goodmidx=[];
    goodtcorr=[];
    gooc_osi=[];
    gooc_dsi=[];
    goozp=[];
    goozc=[];
    goopi=[];
    gooc_zp=[];
    gooc_zc=[];
    gooc_pi=[];
    goodprds=[];
    goodsigtrials=[];
    manual_rejection_list=[];
    
    % create selected sessions vector containing all sessions of selected area
    switch target_group
        case 'V1'
            targetlabel=0;
            selectedsessions=find(cell2mat(listSessions(4,:))==0);
        case 'LM'
            targetlabel=1;
            selectedsessions=find(cell2mat(listSessions(4,:))==1);
%         case 'AL'
%             targetlabel=0;
%             selectedsessions=find(cell2mat(listSessions(4,:))==99);
        case 'RL'
            targetlabel=2;
            selectedsessions=find(cell2mat(listSessions(4,:))==2);
    end
    
    % loop over neurons and conditions
    alreadytaken=0;
    for k=1:size(OSI,1)
        for i=1:size(OSI,2)
            for j=1:size(OSI,3)
                
                if max( sum(k==goodidxtot) ) ...
                        && tuning_matrix(i,j,k)== max(max(tuning_matrix(:,:,k))) ...
                        && not(sum(k==manual_rejection_list)) ...
                        && not(sum(k==alreadytaken(:,1)))%... 
%                         && sum(M(k,1)==selectedsessions) ...                   
%                         && M(k,7)>=targetlabel ...
%                         && sum(sigtrials(:,k,i,j,1))>2 ...
%                         && max(nanmean(trial_matrix(:,:,k,i,j,1),2))>2

                    [~,idxmaxdir]=randmax(c_tuning_curve(:,k,i,j,1));
                    [~,idxmaxdir_p]=randmax(c_tuning_curve(:,k,i,j,2));
                    
                    gooosi=[gooosi;OSI(k,i,j,1)];
                    goodsi=[goodsi;DSI(k,i,j,1)];
                    selectedsi=[selectedsi;[k,i,j]];
                    alreadytaken=cell2mat(selectedsi);
                    goodsf=[goodsf;SF(i)];
                    goodtf=[goodtf;TF(j)];
                    gooddir=[gooddir;idxmaxdir];
                    gooddir_p=[gooddir_p;idxmaxdir_p];
                    goodmidx=[goodmidx;modulation_index(idxmaxdir,k,i,j,1)];
                    goodtcorr=[goodtcorr;tcorr(k,i,j)];
                    gooc_osi=[gooc_osi;c_OSI(k,i,j,1)];
                    gooc_dsi=[gooc_dsi;c_DSI(k,i,j,1)];
                    goozp=[goozp;Zp(k,i,j)];
                    goozc=[goozc;Zc(k,i,j)];
                    goopi=[goopi;PI(k,i,j)];
                    gooc_zp=[goozp;c_Zp(k,i,j)];
                    gooc_zc=[goozc;c_Zc(k,i,j)];
                    gooc_pi=[gooc_pi;c_PI(k,i,j)];                    
%                     goodprds=[goodprds;(c_DI(k,i,j,1)>=0.5 && max(c_tuning_curve_z(:,k,i,j,2))>=3)];
                    goodprds=[goodprds;(c_DSI(k,i,j,1)>=0.33 && max(c_tuning_curve_z(:,k,i,j,2))>=3)];
                    goodsigtrials=[goodsigtrials,sum(sigtrials(:,k,i,j,1))];
                    
                    % count selected neurons
                    neuroncounter=neuroncounter+1;
                    
                else
                end
            end
        end
    end
    
    fprintf(['\n--- ',target_group,' selected neurons = ',num2str(neuroncounter),' ---\n'])
    
    % save population analysis results
    save(['selected_population_',target_group],'gooosi','goodsi','selectedsi'...
        ,'goodsf','goodtf','gooddir','goodmidx'...
        ,'goozp','goozc','goopi','goodprds');
    
    
    %%  -------------------------- results plotting --------------------------

    %% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    % sostituire con function "mosaic...'
%     hhh=1;
%     n=goodn2n('LM',hhh);
%     n_sn=[listSessions{1,M(n,1)},'_b',num2str(M(n,2))];
%     finame=['SPIKEMAT_',n_sn];
%     S=load(finame);
%     SFind=selectedsi{hhh}(2);
%     sf=SF(SFind);
%     TFind=selectedsi{hhh}(3);
%     tf=TF(TFind);
%     dir=DIR(gooddir(hhh));
%     selected_idx = get_indexes_PN( sf, tf, dir, 'grating' );
%     % find whidthin block index
%     ssidx=find(M(:,1)==M(n,1)); % indexes of this session
%     blidx=intersect(find(M(:,2)==M(n,2)),ssidx);  % indexes of this block
%     nnum=find(blidx==n); % find position of selected n
%     % collect spike times
%     S_ts=[];
%     for ii=1:size(S.SPIKEtimes_bis,3)
%         S_ts=[S_ts;S.SPIKEtimes_bis{ selected_idx, nnum, ii}-0.8];
%     end
%     T_s=0.010;
%     hedges=-0.8:T_s:1.7;
%     % produce psth
%     [psth]=hist(S_ts,hedges);
%     figure; plot(hedges,psth)
%     % collect spike times
%     figure;
%     hold on
%     spp=cell(1,20);
%     for trial=1:20
%         if trial<=size(S.SPIKEtimes_bis( selected_idx, nnum, :),3)
%             sp=S.SPIKEtimes_bis{ selected_idx, nnum, trial}-0.8;
%         else
%             sp=[];
%         end
%         if not(isempty(sp))
%             plot(sp,1.05+0.05*trial,'.k', 'MarkerSize',15)
%         else
%         end
%         box on
%         spp{trial}=sp;
%     end
%     raster=spp;
%% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    
    % plotting osi and dsi distribution
    plot_OSI_DSI_distribution( gooc_osi, gooc_dsi, target_group );
    
    % plotting preferred sf and tf distribution
    plot_SF_TF_distribution( goodsf, goodtf, SF, TF, target_group );
    
    % plotting mi distribution
    plot_MI_distribution( goodmidx, target_group );
     
    % plotting pattern and component analysis
    plot_pattern_component( gooc_zc, gooc_zp, gooc_pi, goodprds, selectedsi, target_group );
    

    % plotting tuning curves - single neurons
    plot_tuning_curves( selectedsi,c_tuning_curve,c_tuning_curve_z,response_counts,c_OSI,c_DSI,Zc,Zp,PI,modulation_index, sigtrials, tcorr, basal_fr )
     
% %     % plotting tuning curves - roseplot mosaics
% %     plot_roseplots( selectedsi, gooddir, target_group );
    
    % plotting PSTH and power spectrum
    plot_rasters( selectedsi, gooddir, gooddir_p )
        
    cd(oldfold);
    
end
% %%
% nnn=198;
% aa=cell2mat(selectedsi);
% gnnn=(find(aa(:,1)==nnn));
% 
% 
% figure; subplot(1,2,1); imagesc(tuning_curve_matrix(:,:,nnn,aa(gnnn,2),aa(gnnn,3),1)'); colorbar; axis equal; title('grating');
% subplot(1,2,2); imagesc(tuning_curve_matrix(:,:,nnn,aa(gnnn,2),aa(gnnn,3),2)'); colorbar; axis equal; title('plaid');
% 
% figure; plot(tuning_curve(:,nnn,aa(gnnn,2),aa(gnnn,3),1),'k'); hold on; plot(tuning_curve(:,nnn,aa(gnnn,2),aa(gnnn,3),2),'--k');
% plot(c_tuning_curve(:,nnn,aa(gnnn,2),aa(gnnn,3),1),'b'); 
%     plot(c_tuning_curve(:,nnn,aa(gnnn,2),aa(gnnn,3),2),'--b');