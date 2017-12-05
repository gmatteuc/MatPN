function [] = plot_rasters( selectedsi, gooddir, gooddir_p )

% plot_rasters( selectedsi, gooddir )
%
% plot grating and plaid rasters
%-----------------------------------------------------------------------

% set pars
pars=set_pars_PN;
data_folder=pars.processed_data_folder;
load(fullfile(data_folder,'Tuning.mat'));
load(fullfile(data_folder,'Indexing.mat'));
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
stimulustypes={'grating','plaid'};
pre_time=pars.stimPars.pre_time;
post_time=pars.stimPars.post_time;
stim_time=pars.stimPars.stim_time;
post_delay=pars.stimPars.post_delay;
ntrials=pars.stimPars.ntrials;

for goondn=1:size(selectedsi,1)
    for stimidx=1:2
        
        % switch colors
        if stimidx==1
            colo=[0,0,0.9]; % grating
        else
            colo=[0,0.7,0]; % plaids
        end
        
        % get current neuron parameters
        stimulustype=stimulustypes{stimidx};
        nn=selectedsi{goondn}(1);
        SFind=selectedsi{goondn}(2);
        sf=SF(SFind);
        TFind=selectedsi{goondn}(3);
        tf=TF(TFind);
        % switch colors
        if stimidx==1
            dir=DIR(gooddir(goondn)); % grating
        else
            dir=DIR(gooddir_p(goondn)); % plaids
        end
        
        % select bitcodes for chosen condition
        selected_idx = get_indexes_PN( sf, tf, dir, stimulustype  );
        % retrive within session index
        [intn,sessionname]=n2intn(nn);
        % get raster data
        finame=['SPIKEMAT_',sessionname];
        S=load(finame);
        
        % rebuild psth and raster
        S_ts=[];
        for i=1:size(S.SPIKEtimes_bis,3)
            S_ts=[S_ts;S.SPIKEtimes_bis{ selected_idx, intn, i}-pre_time];
        end
        % take spikes falling into analysis window only
        S_ts=S_ts(S_ts<=(stim_time+1*post_delay));
        S_ts=S_ts(S_ts>=0);
        % set sampling frequency
        T_s=0.010;
        hedges=0:T_s:(stim_time+1*post_delay);
        % bin psth
        [psth]=hist(S_ts,hedges);
        % compute F1z
        [ ~, pow, f, N_f ] = get_F1z( psth, tf, T_s );
        % get index of target frequency
        fidx=find(f<=tf);
        fidx=fidx(end);
        sigspect=std(pow);
        meanspect=mean(pow);
        
        % plot result
        f1 = figure;
        set(f1,'Position',[10,10,1500,1000]);
        
        % draw raster
        subplot(1,2,1);
        bbbb=bar(hedges,5*psth/(max(psth)));
        set(bbbb,'FaceColor',colo)
        hold on
        spp=[];
        for trial=1:ntrials
            if trial<=size(S.SPIKEtimes_bis( selected_idx, intn, :),3)
                sp=S.SPIKEtimes_bis{ selected_idx, intn, trial}-pre_time;
            else
                sp=[];
            end
            if not(isempty(sp))
                plot(sp,7.5+0.5*trial,'.k', 'MarkerSize',15)
            else
            end
            box on
            spp=cat(1,spp,sp);
        end
        xlim([-pre_time,stim_time+post_time]);
        ylim([-5,ntrials+5]);
        plot([0,0],[-5,ntrials+5],'--', 'LineWidth',3,'Color',colo)
        plot([stim_time,stim_time],[-5,ntrials+5],'--', 'LineWidth',3,'Color',colo)
        text(-0.5,22,['DIR = ',num2str(dir),' ',' number of spikes = ',num2str(round(ntrials*response_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),...
            ' (',num2str(round(ntrials*spontaneous_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),')'],'FontSize',12);
        set(gca,'FontSize',10);
        hlabelx=get(gca,'Xlabel');
        set(hlabelx,'String','time from onset [s]','FontWeight','bold','FontSize',12,'color','k')
        hlabely=get(gca,'Ylabel');
        set(hlabely,'String','trial number','FontWeight','bold','FontSize',12,'color','k')
        axis square
        sname=strrep(sessionname,'_',' ');
        title(['PSTH n ',num2str(nn),' goodn ',num2str(goondn),' - session ',sname]);
        
        % draw spectrum
        subplot(1,2,2);
        plot(f(1:N_f/2+1),pow(1:N_f/2+1),'k','LineWidth',3)
        hold on
        plot(f(fidx),pow(fidx),'o','LineWidth',6,'Color',colo);
        plot(f(1:N_f/2+1),(meanspect+sigspect)*ones(size(pow(1:N_f/2+1))),'--','LineWidth',1,'Color',colo);
        plot(f(1:N_f/2+1),(meanspect-sigspect)*ones(size(pow(1:N_f/2+1))),'--','LineWidth',1,'Color',colo);
        plot(f(1:N_f/2+1),(meanspect)*ones(size(pow(1:N_f/2+1))),'.-','LineWidth',0.5,'Color',colo);
        ylimit=get(gca,'ylim');
        xlimit=get(gca,'xlim');
        text(0.6*xlimit(2),0.85*ylimit(2),['target TF = ',num2str(tf),' Hz'],'FontSize',14);
        set(gca,'FontSize',10);
        hlabelx=get(gca,'Xlabel');
        set(hlabelx,'String','f [Hz]','FontWeight','bold','FontSize',10,'color','k')
        hlabely=get(gca,'Ylabel');
        set(hlabely,'String','PSD','FontWeight','bold','FontSize',10,'color','k')
        hold off
        axis square
        F1z_orig=modulation_index(gooddir(goondn), nn, find(SF==sf), find(TF==tf), stimidx);
        F1F0_orig=modulation_index_bis(gooddir(goondn), nn, find(SF==sf), find(TF==tf), stimidx);
        % save plot
        title(['Power spectrum (F1z = ',num2str(F1z_orig),', F1F0 = ',num2str(F1F0_orig),') ',stimulustype]);
        fname=[stimulustype,' modulation analysis ',' n_',num2str(nn)];
        fname=strrep(fname,'.','');
        fname=strrep(fname,' ','_');
        set(gcf, 'PaperPositionMode', 'auto')
        saveas(gcf,fname, 'jpg')
 
    % ---------------------------------------------------------------------
    
    % set plotname
    fname=['raster mosaic ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' n',num2str(nn),''];
    fname=strrep(fname,'.','');
    fname=strrep(fname,' ','_');
    
    % initialize figure
    f1 = figure;
    set(f1,'Position',[10,10,1500,1000]);
    poscount1=0;
    poscount2=0;
    
    for dind=1:length(DIR)
        
        % select direction
        dir=DIR(dind);
        % select bitcode
        selected_idx = get_indexes_PN( sf, tf, dir, stimulustype );
        
        % initialize subplot
        sb1=subplot(666,666,666);
        hold on
        
        spp=cell(1,ntrials);
        for trial=1:ntrials
            if trial<=size(S.SPIKEtimes_bis( selected_idx, intn, :),3)
                sp=S.SPIKEtimes_bis{ selected_idx, intn, trial}-pre_time;
            else
                sp=[];
            end
            if not(isempty(sp))
                plot(sp,trial,'.k', 'MarkerSize',20)
            else
            end
            box on
            spp{trial}=sp;
        end
        xlim([-pre_time,post_time+stim_time]);
        ylim([-5,ntrials+5]);
        plot([0,0],[-5,ntrials+5],'--', 'LineWidth',3,'Color',colo)
        plot([0.9,0.9],[-5,25],'--', 'LineWidth',3,'Color',colo)
        text(-0.5,22,['DIR = ',num2str(dir),' ',' number of spikes = ',num2str(round(ntrials*response_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),...
            ' (',num2str(round(ntrials*spontaneous_fr(dir==DIR,nn,sf==SF,tf==TF,stimidx))),')'],'FontSize',12);
        set(gca,'FontSize',10);
        hlabelx=get(gca,'Xlabel');
        set(hlabelx,'String','time from onset [s]','FontWeight','bold','FontSize',12,'color','k')
        hlabely=get(gca,'Ylabel');
        set(hlabely,'String','trial number','FontWeight','bold','FontSize',12,'color','k')
            
        % displace subplot to correct position
        hold off
        set(sb1,'Position',[.00+0.23*(poscount1)+0.07,.01+0.32*(poscount2)+0.1,0.20,0.20]);
        poscount1=poscount1+1;
        if poscount1==4 % change row after 4 subplots
            poscount2=poscount2+1;
            poscount1=0;
        end
        
    end
    
     % ---------------------------------------------------------------------
     
     % save plot
     saveas(f1,fname,'jpg')
     close all
    end
    
    message=['\nNeuron ',num2str(intn),' raster mosaic produced\n'];
    fprintf(message)
    clear S
    
end
end

