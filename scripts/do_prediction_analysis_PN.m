% ------------------------- DO PREDICTION ANALYSIS -------------------------

clear
close all
clc

pars = set_pars_PN();

% set stimulus parameters
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
types=pars.neuPars.ctypes;
frame_dur=pars.stimPars.frame_duration;
stim_dur=pars.stimPars.stim_time;
prediction_isi=pars.prediction_isi;

% set analysis parameters
maxlag=pars.STA_depth;
stimwidth=pars.STA_width;
stimheight=pars.STA_height;
interp_factor=pars.interp_factor;
stimlength=pars.stimPars.noiselength;
crop_pixel_size=pars.crop_pixel_size;
contrast_th=pars.contrast_threshold; % 25/07/2017 version used 5.5, now (01/12/2017) 0

% set useful paths
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);
cd(data_folder);
stimulustypes=pars.stimPars.stimulustype;

% load tuning analysis results and indexing
load('Tuning.mat')
load('Indexing.mat')

% loop over areas
target_areas={'V1','LM','RL'};
for target_area_idx=2:length(target_areas)
    
    % set area
    are=target_areas{target_area_idx};
    
    % load neuron selection results
    sel_res_file=[are,'_selected_neurons',filesep,'selected_population_',are,'.mat'];
    load(sel_res_file)
    
    % sta analysis results
    sta_res_file=[data_folder,filesep,'RFs_datastructure_',are];
    load(sta_res_file)
    
    % create output folder
    outfold=[are,'_prediction_results'];
    outfold=strrep(outfold,'','_');
    mkdir(outfold);
    oldfold=cd(outfold);
    
    % select good sta neurons onto which perform the analysis
    bgosta=zeros(size(gostac));
    bgosta(gostac>contrast_th)=1;
    
    % initialize prediction storage variables
    stimulus_binnumber=floor(stim_dur/frame_dur);
    isi_binnumber=prediction_isi;
    prediction_binnumber=isi_binnumber + stimulus_binnumber + isi_binnumber;
    pred_FR=zeros(prediction_binnumber,length(DIR),length(selectedsi),numel(stimulustypes));
    pred_COUNT=zeros(length(DIR),length(selectedsi),numel(stimulustypes));
    filter=zeros(stimheight*interp_factor,stimwidth*interp_factor,maxlag);
    explvar_sD=zeros(length(selectedsi),numel(stimulustypes));
    explvar_D=zeros(length(selectedsi),numel(stimulustypes));
    explvar_T=zeros(length(selectedsi),numel(stimulustypes));
    pDIRdiff=zeros(length(selectedsi),numel(stimulustypes));
    pORIdiff=zeros(length(selectedsi),numel(stimulustypes));
    PI_O=zeros(length(selectedsi),1);
    Zp_O=zeros(length(selectedsi),1);
    Zc_O=zeros(length(selectedsi),1);
    Rp_O=zeros(length(selectedsi),1);
    Rc_O=zeros(length(selectedsi),1);
    PI_P=zeros(length(selectedsi),1);
    Zp_P=zeros(length(selectedsi),1);
    Zc_P=zeros(length(selectedsi),1);
    Rp_P=zeros(length(selectedsi),1);
    Rc_P=zeros(length(selectedsi),1);
    good_contrast=zeros(length(selectedsi),1);
    good_power=zeros(length(selectedsi),1);
    
    % loop over selected neurons 
    for i=1:length(selectedsi)
        
        % get current neuron STA ------------------------------------------------------------------
        selectedSTA=Zwsta(:,:,:,i);
        n=goodn2n(are,i);
        
        % loop over stimulus types
        for k=1:numel(stimulustypes)
            
            % get the preferred stimulus ------------------------------------------------------------------
            
            % get preferred condition
            pSF=SF(1);
            pTF=TF(selectedsi{i}(3));
            if k==2 % if plaid
                [~,idxmax]=randmax(tuning_curve(:,n,pSF==SF,pTF==TF,k));
                pDIR=DIR(idxmax);
            else % if grating
                [~,idxmax]=randmax(tuning_curve(:,n,pSF==SF,pTF==TF,k));
                pDIR=DIR(idxmax);
            end
            
            % load stimulus frames
            stimulustype=stimulustypes{k};
            pstimul = get_stimulus_PN( pSF, pTF, pDIR, stimulustype );
            
            % get the neural filter ------------------------------------------------------------------
            
            % set cropping params
            sigma=min([size(selectedSTA,1)*interp_factor,size(selectedSTA,2)*interp_factor])/6;
            filter = zeros(crop_pixel_size,crop_pixel_size,size(selectedSTA,3));
            
            % decide cropping window on best frame
            inputfr = selectedSTA(:,:,bestfr_contrast(1));
            binp = interpolate_RF_frame( inputfr, interp_factor );
            [ ~, crop_ridx, crop_cidx  ] = apply_crop( binp,crop_pixel_size,0 );
            
            for ff=1:size(selectedSTA,3)
                
                inputfr = selectedSTA(:,:,ff);
                binp = interpolate_RF_frame( inputfr, interp_factor );
                % stack into filter matrix
                filter(:,:,ff) = binp(crop_ridx,crop_cidx);
                
            end
            
            % compute predictions for every direction to produce a tuning curve ------------------------------------------------------------------
            
            for j=1:length(DIR)
                
                % get stimulus
                stimul = get_stimulus_PN( pSF, pTF, DIR(j),stimulustype );
                % crop stimulus
                inp = stimul(crop_ridx,crop_cidx,:);
                stimul = inp;
                % compute prediction
                [pred_FR(:,j,i,k),pred_COUNT(j,i,k),pred_FR_time]=neural_filtering(filter,stimul);
                
                fprintf(['filter response at ',stimulustype,' ( DIR=',num2str(DIR(j)),' SF=',num2str(pSF),' TF=',num2str(pTF),' ) neuron ',num2str(n),' produced ...\n']);
            end
            
            % evaluate goodness of prediction   ------------------------------------------------------------------
            
            % get observed psth and raster
            [ psth, raster, hedges ] = get_psth_PN(  n, pSF, pTF, pDIR, stimulustype  );
            % extract firing rate prediction
            tinf=find(pred_FR_time<1);
            tsup=find(pred_FR_time>0);
            tind=intersect(tinf,tsup);
            spred_FR=pred_FR(:,DIR==pDIR,i,k);
            spred_FRth=max((spred_FR(tind)/max(spred_FR(:))),zeros(length(tind),1));
            % resample psth at the same sr of the prediction
            new_samplenum=length(spred_FRth);
            old_samplenum=length(psth);
            psth_resampled = resample(psth, new_samplenum,old_samplenum);
            
            % better to resample psth at a lower sr than doing the opposite
            % resample psth at the same sr of the prediction after smoothing
            % (in order to make dynamics timescale comparable and reduce fluctuations)
            
            % smooth to match timescales
            gaussFilter = gausswin(10,1/1.1);
            gaussFilter = gaussFilter / sum(gaussFilter);
            spsth=conv(psth, gaussFilter,'same');
            new_samplenum=length(spred_FRth);
            old_samplenum=length(psth);
            spsth_resampled = resample(spsth, new_samplenum,old_samplenum);
            
            % compare smoothed observed and predicted dynamics
            explvar_sdynamics=corr(spred_FRth,spsth_resampled').^2;
            if isnan(explvar_sdynamics)
                explvar_sdynamics=0;
            else
            end
            explvar_sD(i,k)=explvar_sdynamics;
            
            % compare observed and predicted dynamics
            explvar_dynamics=corr(spred_FRth,psth_resampled').^2;
            if isnan(explvar_dynamics)
                explvar_dynamics=0;
            else
            end
            explvar_D(i,k)=explvar_dynamics;
            
            % compare observed and predicted tuning curve
            explvar_tuning=(corr(tuning_curve( :,n,pSF==SF,pTF==TF,k),pred_COUNT(:,i,k))).^2;
            if isnan(explvar_tuning)
                explvar_tuning=0;
            else
            end
            explvar_T(i,k)=explvar_tuning;
            
            % compare observed and predicted preferred direction
            pred_pDIR=DIR(pred_COUNT(:,i,k)==max(pred_COUNT(:,i,k)));
            pred_pDIR=pred_pDIR(1);
            pDIRdiff(i,k)=rad2deg(angleDiff(deg2rad(pDIR),deg2rad(pred_pDIR)));
            
            % compare observed and predicted preferred orientation
            if rad2deg(angleDiff(deg2rad(pDIR),deg2rad(pred_pDIR)))>=90
                pORIdiff(i,k)=rad2deg(angleDiff(deg2rad(pDIR+180),deg2rad(pred_pDIR)));
            elseif rad2deg(angleDiff(deg2rad(pDIR),deg2rad(pred_pDIR)))<-90
                pORIdiff(i,k)=rad2deg(angleDiff(deg2rad(pDIR-180),deg2rad(pred_pDIR)));
            else
                pORIdiff(i,k)=rad2deg(angleDiff(deg2rad(pDIR),deg2rad(pred_pDIR)));
            end
            
            % compute predicted pattern index  ------------------------------------------------------------------
            
            if k==2 %do this only when you are analizing plaids (before grating)
                tuning_curve_grating_P=pred_COUNT(:,i,1);
                tuning_curve_plaid_P=pred_COUNT(:,i,2);
                PI_O(i)=PI(n,pSF==SF,pTF==TF);
                Zp_O(i)=Zp(n,pSF==SF,pTF==TF);
                Zc_O(i)=Zc(n,pSF==SF,pTF==TF);
                Rp_O(i)=Rp(n,pSF==SF,pTF==TF);
                Rc_O(i)=Rc(n,pSF==SF,pTF==TF);
                [ PI_P(i), Zp_P(i), Zc_P(i), Rp_P(i), Rc_P(i) ] = get_pattern_index( tuning_curve_grating_P,tuning_curve_plaid_P );
            else
            end
            
            % plot predictions ------------------------------------------------------------------
            
            f1 = figure;
            set(f1,'Position',[10,10,1500,1000]);
            sb1=subplot(666,666,1);
            % plot PSTH
            T_s=1/length(psth);
            ds=psth-mean(psth);
            time=0:T_s:(length(ds)-1)*T_s;
            area(hedges,psth/(max(psth)),'EdgeColor','k','LineWidth',2.5,'FaceColor','k');
            hold on
            if k==2
                plot(pred_FR_time(tind),spred_FRth,'-b','LineWidth',3.5)
            else
                plot(pred_FR_time(tind),spred_FRth,'-g','LineWidth',3.5)
            end
            plot(pred_FR_time(tind),max(spred_FRth(:)).*(spsth_resampled/(max(spsth_resampled))),'--k','LineWidth',3.5)
            % draw raster
            spp=raster;
            spk=0;
            for kkk=1:numel(raster)
                if not(isempty(raster{kkk}))
                    plot(raster{kkk},1.05+0.05*kkk,'.k', 'MarkerSize',15)
                    spk=spk+numel(raster{kkk});
                else
                end
            end
            xlim([-0.2,1.2]);
            ylim([-0,4]);
            plot([0,0],[0,5],'--k', 'LineWidth',2)
            plot([1,1],[0,5],'--k', 'LineWidth',2)
            tt=text(0.05,3.5,['DIR = ',num2str(pDIR),' ',' spike count = ',num2str(spk)],'FontSize',14);
            set(gca,'FontSize',10);
            str = sprintf(['expl. var. = ',num2str(explvar_sD(i,k),'%.2f')]);
            tx=text(0.05,3.0,str);
            hlabelx=get(gca,'Xlabel');
            set(hlabelx,'String','time [s]','FontWeight','bold','FontSize',10,'color','k')
            hlabely=get(gca,'Ylabel');
            set(hlabely,'String','normalized firing rate','FontWeight','bold','FontSize',10,'color','k')
            axis square
            % select bitcodes for the tuning curve
            selected_idx = get_indexes_PN( pSF, pTF, pDIR, stimulustype );
            sessionname=[listSessions{1,M(n,1)},'_b',num2str(M(n,2))];
            sname=strrep(sessionname,'_',' ');
            title([stimulustype,' - PSTH neuron ',num2str(n),' session ',sname])
            set(sb1,'Position',[.52,0.4,.5,.5]);
            
            sb2=subplot(666,666,2);
            predv=[squeeze(pred_COUNT(:,i,k))'./max(squeeze(pred_COUNT(:,i,k))),squeeze(pred_COUNT(1,i,k))'./max(squeeze(pred_COUNT(:,i,k)))];
            obsv=[tuning_curve( :,n,pSF==SF,pTF==TF,k)'./max(squeeze(tuning_curve( :,n,pSF==SF,pTF==TF,k))),tuning_curve( 1,n,pSF==SF,pTF==TF,k )./max(squeeze(tuning_curve( :,n,pSF==SF,pTF==TF,k )))];
            
            if k==2
                p1=polar([degtorad(DIR),2*pi],predv,'-b');
            else
                p1=polar([degtorad(DIR),2*pi],predv,'-g');
            end
            hold on
            p2=polar([degtorad(DIR),2*pi],obsv,'-k');
            set(p1, 'linewidth', 3.5);
            set(p2, 'linewidth', 3.5);
            str = sprintf(['expl. var. = ',num2str(explvar_T(i,k),'%.2f')]);
            tx=text(0.5,1.1,str);
            str2 = sprintf(['pref. dir. diff. = ',num2str(pDIRdiff(i,k),'%.0f')]);
            tx2=text(0.8,0.9,str2);
            str3 = sprintf(['pref. ori. diff. = ',num2str(pORIdiff(i,k),'%.0f')]);
            tx3=text(1.1,0.7,str3);
            % select best frame
            imaxfp=bestfr_power(i);
            if k==2
                str4 = sprintf(['Zp_O = ',num2str(Zp_O(i),'%.01f'),' Zc_O = ',num2str(Zc_O(i),'%.01f')]);
                tx4=text(1.4,0.5,str4);
                str5 = sprintf(['Zp_P = ',num2str(Zp_P(i),'%.01f'),' Zc_P = ',num2str(Zc_P(i),'%.01f')]);
                tx5=text(1.4,0.3,str5);
                str6 = sprintf(['obs. gDSI = ',num2str(DSI(selectedsi{i}(1),selectedsi{i}(2),selectedsi{i}(3),1),'%.02f')]);
                tx6=text(1.4,0.1,str6);
%                 str7 = sprintf(['ori = ',num2str(orientation(imaxfp,i),'%.02f'),' deg at f',num2str(imaxfp,'%.0f')]);
%                 tx7=text(1.4,-0.1,str7);
%                 str8 = sprintf(['apect ratio = ',num2str(aspect_ratio(imaxfp,i),'%.02f')]);
%                 tx8=text(1.4,-0.3,str8);
%                 str9 = sprintf(['eccentricity = ',num2str(eccentricity(imaxfp,i),'%.02f')]);
%                 tx9=text(1.4,-0.5,str9);
            else
            end
            
            axis off
            set(sb2,'Position',[.02,0.4,.5,.5]);
            
            hold on
            
            % loop over frames
            for jj=1:size(selectedSTA,3)
                sb3=subplot(666,666,1);
                fram=interpolate_RF_frame( selectedSTA(:,:,jj), interp_factor );
                fram=fram(crop_ridx,crop_cidx);
                l1=imagesc(fram); colormap('gray'); caxis([-6, 6]); set(gca,'dataAspectRatio',[1 1 1]); axis off
                set(sb3,'Position',[.02+0.095*(jj-1),0.1,.09,.09]);
            end
            hold on
            % loop over frames
            for jj=1:size(selectedSTA,3)
                sb4=subplot(666,666,1);
                fram=pstimul(crop_ridx,crop_cidx,25+jj);
                imagesc(fram); colormap(gray); set(gca,'dataAspectRatio',[1 1 1]); axis off
                set(sb4,'Position',[.02+0.095*(jj-1),0.2,.09,.09]);
                ylimit=get(gca,'ylim');
                xlimit=get(gca,'xlim');
                strrsq = sprintf(['STA cont = ',num2str(contrast(jj,i),'%.2f')]);
                tx3=text(0.03*xlimit(2),2.30*ylimit(2),strrsq);
                ylimit=get(gca,'ylim');
                xlimit=get(gca,'xlim');
                strpow = sprintf(['STA pow = ',num2str(power(jj,i),'%.1f')]);
                tx4=text(0.03*xlimit(2),-.21*ylimit(2),strpow);
            end
            
            % save
            fname=['prediction ','(n_',num2str(n),' goodn_',num2str(i),') ',are,' ',stimulustypes{k}];
            fname=strrep(fname,'.','');
            % saveas(gcf,fname, 'epsc')
            set(gcf, 'PaperPositionMode', 'auto')
            saveas(gcf,fname, 'jpg')
            
        end
        
        close all
        
    end
    % save results
    save(['prediction_explained_variance_',are],'explvar_T','explvar_D','explvar_sD','pDIRdiff','pORIdiff',...
        'PI_O','Zp_O','Zc_O','Rp_O','Rc_O','PI_P','Zp_P','Zc_P','Rp_P','Rc_P')
    
    % plot disributions
    plot_prediction_results_PN(  are  )
    
    cd(oldfold)
    
end
