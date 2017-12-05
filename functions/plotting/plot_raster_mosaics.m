function [] = plot_raster_mosaics( SFind,TFind,stimulustype )

%Plot tuning curves for every SF and TF of every neuron in every session
%
%--------------------------------------------------------------------------

% load indexing
load('Indexing.mat')

% set pars and paths
pars = set_pars_PN();
SF=pars.stimPars.SF;
listSessions=pars.listSessions;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;

% switch color
colo{1}=[0,0,0.7];
colo{2}=[0,0.7,0];
switch stimulustype
    
    case 'grating'
        coloidx=1;
    case 'plaid'
        coloidx=2;
        
end

numSS=size(pars.listSessions(1,:),2);
for sind=1:numSS
    
    % select session
    session = listSessions{1,sind};
    ssidx=find(M(:,1)==sind);
    blocks = unique(M(ssidx,2));
    numBL=length(blocks);
    
    for bind=1:numBL
        
        % select block
        bbidx=find(M(:,2)==blocks(bind));
        sidx=intersect(ssidx,bbidx);
        
        % get datafile name
        n_sn=[session,'_b',num2str(blocks(bind))];
        finame=['SPIKEMAT_',n_sn];
        S=load(finame);
        
        % get sf and tf values
        sf=SF(SFind);
        tf=TF(TFind);
        
        % create ouput folder
        dirnam=['tuning_results_',session,'_b',num2str(blocks(bind))];
        mkdir(dirnam);
        oldd=cd(dirnam);
        
        for nind=1:length(sidx)
            
            tic
            % set sessionwise neuron index
            intind=sidx(nind);
            
            % set plotname
            fname=['raster mosaic ',stimulustype,' SF=',num2str(sf),' TF=',num2str(tf),' n_',num2str(intind),''];
            fname=strrep(fname,'.','');
            fname=strrep(fname,' ','_');
            if exist([fname,'.jpg'],'file')~=1
                
                % plot mosaic for current neuron --------------------------
                
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
                    
                    spp=cell(1,20);
                    spc=0;
                    spc_spont=0;
                    for trial=1:20
                        if trial<=size(S.SPIKEtimes_bis( selected_idx, intind, :),3)
                            sp=S.SPIKEtimes_bis{ selected_idx, intind, trial}-0.8;
                        else
                            sp=[];
                        end
                        if not(isempty(sp))
                            plot(sp,trial,'.k', 'MarkerSize',20)
                        else
                        end
                        box on
                        spp{trial}=sp;
                        % count spikes outside stimulus
                        sp_tris_pre=sp(sp<0);
                        sp_tris_post=sp(sp>1);
                        spc_spont=spc_spont+numel(sp_tris_post)+numel(sp_tris_pre);
                        % count spikes during stimulus
                        sp_bis=sp(sp<1);
                        sp_bis=sp_bis(sp_bis>0.05);
                        spc=spc+numel(sp_bis);
                    end
                    xlim([-0.8,1.7]);
                    ylim([-5,25]);
                    plot([0,0],[-5,25],'--', 'LineWidth',3,'Color',colo{coloidx})
                    plot([0.9,0.9],[-5,25],'--', 'LineWidth',3,'Color',colo{coloidx})
                    tt=text(-0.5,22,['DIR = ',num2str(dir),' ',' number of spikes = ',num2str(spc),' (',num2str(round((spc_spont*0.9)/1.5)),')'],'FontSize',12);
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
                
                %----------------------------------------------------------

            else
            end
            
            saveas(f1,fname,'jpg')
            close all
            message=['\nNeuron ',num2str(intind),' raster mosaic produced\n'];
            fprintf(message)
            toc
            
        end
        cd(oldd)
    end
end
end
