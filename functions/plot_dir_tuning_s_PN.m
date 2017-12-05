function [] = plot_dir_tuning_s_PN( SFind,TFind,stimulustype )

%Plot tuning curves for every SF and TF of every neuron in every session
%
%--------------------------------------------------------------------------

%% load data

load('Tuning.mat')
load('Indexing.mat')

%% set pars and paths

pars = set_pars_PN();
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
types=pars.neuPars.ctypes;
code_folder=pars.codeFolder;
addpath(code_folder);
paruly_folder=pars.parulyFolder;
addpath(paruly_folder);

%% switch stimulus type boolean

switch stimulustype
    case 'grating'
        boolst=1;
    case 'plaid'
        boolst=2;
end

%% prepare colors

lev=3;
lev2=6;
cmap = paruly;
cmapV1=cmap(lev,:);
cmapLL=cmap(end-lev2,:);

%% core loop

numSS=size(pars.listSessions(1,:),2);
for sind=1:numSS;
    
    session = listSessions{1,sind};
    ssidx=find(M(:,1)==sind);
    blocks = unique(M(ssidx,2));
    numBL=length(blocks);
    
    for bind=1:numBL;
        
        bbidx=find(M(:,2)==blocks(bind));
        sidx=intersect(ssidx,bbidx);
        
        % create ouput folder
        dirnam=['tuning_results_',session,'_b',num2str(blocks(bind))];
        mkdir(dirnam);
        oldd=cd(dirnam);
        
        %% ||||||||||||||| to loop over cluster types  |||||||||||||||
        for z=1:length(types)
            
            type=types{z};
            
            switch z
                case 1 % good neurons
                    goodidxtot=find(M(:,8)==0);
                    goodidx=intersect(sidx,goodidxtot);
                    clindeces=M(goodidx,9);
                    rtuning_curve=tuning_curve( :,goodidx,:,:,:);
                    rtuning_curve_error=tuning_curve_error( :,goodidx,:,:,:);
                    rDSI=DSI(goodidx,:,:,:);
                    rOSI=OSI(goodidx,:,:,:);
                case 2 % muas
                    muaidxtot=find(M(:,8)==1);
                    muaidx=intersect(sidx,muaidxtot);
                    clindeces=M(muaidx,9);
                    rtuning_curve=tuning_curve( :,muaidx,:,:,:);
                    rtuning_curve_error=tuning_curve_error( :,muaidx,:,:,:);
                    rDSI=DSI(muaidx,:,:,:);
                    rOSI=OSI(muaidx,:,:,:);
                case 3 % noise
                    noiseidxtot=find(M(:,8)==2);
                    noiseidx=intersect(sidx,noiseidxtot);
                    clindeces=M(noiseidx,9);
                    rtuning_curve=tuning_curve( :,noiseidx,:,:,:);
                    rtuning_curve_error=tuning_curve_error( :,noiseidx,:,:,:);
                    rDSI=DSI(noiseidx,:,:,:);
                    rOSI=OSI(noiseidx,:,:,:);
            end
            
            % calculate number of plots needed
            plotnum=floor(numel(clindeces)/16);
            rplotnum=mod(numel(clindeces),16);
            
            switch stimulustype
                case 'grating'
                    
                    % |||||||||||||||||| if plotting GRATINGS ||||||||||||||||||
                    
                    % do complete plots
                    if plotnum~=0
                        for hh=1:plotnum
                            f1 = figure;
                            set(f1,'Position',[10,10,1500,1000]);
                            poscount1=0;
                            poscount2=0;
                            for nn=(1:16)+16*(hh-1)
                                tc=rtuning_curve( :,nn,SFind,TFind,boolst);
                                tce=rtuning_curve_error( :,nn,SFind,TFind,boolst);
                                dsi=rDSI(nn,SFind,TFind,boolst);
                                osi=rOSI(nn,SFind,TFind,boolst);

                                tempmaxdir=find(tuning_curve(:,nn,SFind,TFind,boolst)==tuning_matrix(SFind,TFind,nn,boolst));
                                if size(tempmaxdir)>1
                                    id=randperm(size(tempmaxdir));
                                else
                                    id=1;
                                end
                                mi=modulation_index(nn,SFind,TFind,tempmaxdir(id));

                                sb1=subplot(666,666,666);
                                pp=polar([degtorad(DIR),2*pi],[tc',tc(1)]); set(pp,'Color',cmapV1); set(pp, 'linewidth', 3);
                                hold on;
                                [X1,Y1] = pol2cart([degtorad(DIR),2*pi],[tc',tc(1)]+[tce',tce(1)]); % Convert polar to cartesian coordinates
                                [X2,Y2] = pol2cart([degtorad(DIR),2*pi],[tc',tc(1)]-[tce',tce(1)]); % Convert polar to cartesian coordinates
                                h2=fill([X1,fliplr(X2)],[Y1,fliplr(Y2)],cmapV1); % This will fill in your shape.
                                set(h2,'facealpha',.45); set(h2,'EdgeColor','None');
                                hold off;
                                set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
                                sb2=subplot(666,666,666);
                                tc=tc+0.00001*rand(size(tc));
                                if z==1
                                str = sprintf(['n ',num2str(goodidx(nn)),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nF1z=',num2str(mi,'%0.1f')]);
                                else
                                str = sprintf(['norder ',num2str(nn),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nF1z=',num2str(mi,'%0.1f')]);
                                end
                                tx=text(0.5,0.1,str); axis off
                                set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
                                poscount1=poscount1+1;
                                if poscount1==4;
                                    poscount2=poscount2+1;
                                    poscount1=0;
                                end
                            end
                            h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(hh),')']);
                            set(gca, 'Visible', 'off');
                            set(h, 'Visible', 'on', 'FontSize', 15);
                            set(gcf, 'Color', 'w');
                            set(gcf, 'PaperPositionMode', 'auto')
                            fname=['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(hh),')'];
                            fname=strrep(fname,'.','');
                            saveas(gcf,fname, 'jpg')
                           
                        end
                    else
                    end
                    close all
                    
                    % do incomplete plots
                    f1 = figure;
                    set(f1,'Position',[10,10,1500,1000]);
                    poscount1=0;
                    poscount2=0;
                    for nn=16*(plotnum)+1:16*(plotnum)+rplotnum
                        tc=rtuning_curve( :,nn,SFind,TFind,boolst);
                        tce=rtuning_curve_error( :,nn,SFind,TFind,boolst);
                        dsi=rDSI(nn,SFind,TFind,boolst);
                        osi=rOSI(nn,SFind,TFind,boolst);
                        
                        tempmaxdir=find(tuning_curve(:,nn,SFind,TFind,boolst)==tuning_matrix(SFind,TFind,nn,boolst));
                        if size(tempmaxdir)>1
                            id=randperm(size(tempmaxdir));
                        else
                            id=1;
                        end
                        mi=modulation_index(nn,SFind,TFind,tempmaxdir(id));
                        
                        sb1=subplot(666,666,666);
                        pp=polar([degtorad(DIR),2*pi],[tc',tc(1)]); set(pp,'Color',cmapV1); set(pp, 'linewidth', 3);
                        hold on;
                        [X1,Y1] = pol2cart([degtorad(DIR),2*pi],[tc',tc(1)]+[tce',tce(1)]); % Convert polar to cartesian coordinates
                        [X2,Y2] = pol2cart([degtorad(DIR),2*pi],[tc',tc(1)]-[tce',tce(1)]); % Convert polar to cartesian coordinates
                        h2=fill([X1,fliplr(X2)],[Y1,fliplr(Y2)],cmapV1); % This will fill in your shape.
                        set(h2,'facealpha',.45); set(h2,'EdgeColor','None');
                        hold off;
                        set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
                        sb2=subplot(666,666,666);
                        tc=tc+0.00001*rand(size(tc));
                        if z==1
                        str = sprintf(['n ',num2str(goodidx(nn)),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nF1z=',num2str(mi,'%0.1f')]);
                        else
                        str = sprintf(['norder ',num2str(nn),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nF1z=',num2str(mi,'%0.1f')]);
                        end
                        tx=text(0.5,0.1,str); axis off
                        set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
                        poscount1=poscount1+1;
                        if poscount1==4;
                            poscount2=poscount2+1;
                            poscount1=0;
                        end
                    end
                    h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(plotnum+1),')']);
                    set(gca, 'Visible', 'off');
                    set(h, 'Visible', 'on', 'FontSize', 15);
                    set(gcf, 'Color', 'w');
                    set(gcf, 'PaperPositionMode', 'auto')
                    fname=['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(plotnum+1),')'];
                    fname=strrep(fname,'.','');
                    saveas(gcf,fname, 'jpg')
                    saveas(gcf,fname, 'epsc')
                    
                case 'plaid'
                    
                    % |||||||||||||||||| if plotting PLAIDS ||||||||||||||||||
                    
                    % do complete plots
                    if plotnum~=0
                        for hh=1:plotnum
                            f1 = figure;
                            set(f1,'Position',[10,10,1500,1000]);
                            poscount1=0;
                            poscount2=0;
                            for nn=(1:16)+16*(hh-1)
                                tc=rtuning_curve( :,nn,SFind,TFind,boolst);
                                tce=rtuning_curve_error( :,nn,SFind,TFind,boolst);
                                tcg=rtuning_curve( :,nn,SFind,TFind,1);
                                tceg=rtuning_curve_error( :,nn,SFind,TFind,1);
                                dsi=rDSI(nn,SFind,TFind,boolst);
                                osi=rOSI(nn,SFind,TFind,boolst);
                                
                                sb1=subplot(666,666,666);
                                
                                data=tc;
                                data=[data',data(1)];
                                databis=tcg;
                                databis=[databis',databis(1)];
                                centers=deg2rad(DIR);
                                centers=[centers,centers(1)];
                                colmap=winter;
                                alphaval=0.5;
                                
                                handle=cart2rose(centers,databis);
                                xh = get(handle,'Xdata');
                                yh = get(handle,'Ydata');
                                g=patch(xh,yh,'y');
                                par=colmap;
                                set(g,'FaceColor',par(1,:));
                                set(g,'FaceAlpha',alphaval);
                                set(g,'EdgeColor',par(1,:)/2.1);
                                set(g,'LineWidth',2);
                                
                                hold on
                                
                                handle=cart2rose(centers,data);
                                xh = get(handle,'Xdata');
                                yh = get(handle,'Ydata');
                                g=patch(xh,yh,'y');
                                par=colmap;
                                set(g,'FaceColor',par(end,:));
                                set(g,'FaceAlpha',alphaval);
                                set(g,'EdgeColor',par(end,:)/2.1);
                                set(g,'LineWidth',2);
                                
                                set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
                                sb2=subplot(666,666,666);
                                tc=tc+0.00001*rand(size(tc));
                                if z==1
                                str = sprintf(['n ',num2str(goodidx(nn)),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind))]);
                                else
                                str = sprintf(['norder ',num2str(nn),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind))]);
                                end
                                tx=text(0.5,0.1,str); axis off
                                set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
                                poscount1=poscount1+1;
                                if poscount1==4;
                                    poscount2=poscount2+1;
                                    poscount1=0;
                                end
                            end
                            h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(hh),')']);
                            set(gca, 'Visible', 'off');
                            set(h, 'Visible', 'on', 'FontSize', 15);
                            set(gcf, 'Color', 'w');
                            set(gcf, 'PaperPositionMode', 'auto')
                            fname=['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(hh),')'];
                            fname=strrep(fname,'.','');
                            saveas(gcf,fname, 'jpg')
                            saveas(gcf,fname, 'epsc')
                        end
                    else
                    end
                    close all
                    
                    % do incomplete plots
                    f1 = figure;
                    set(f1,'Position',[10,10,1500,1000]);
                    poscount1=0;
                    poscount2=0;
                    for nn=16*(plotnum)+1:16*(plotnum)+rplotnum
                        tc=rtuning_curve( :,nn,SFind,TFind,boolst);
                        tce=rtuning_curve_error( :,nn,SFind,TFind,boolst);
                        tcg=rtuning_curve( :,nn,SFind,TFind,1);
                        tceg=rtuning_curve_error( :,nn,SFind,TFind,1);
                        dsi=rDSI(nn,SFind,TFind,boolst);
                        osi=rOSI(nn,SFind,TFind,boolst);
                        sb1=subplot(666,666,666);
                        
                        data=tc;
                        data=[data',data(1)];
                        databis=tcg;
                        databis=[databis',databis(1)];
                        centers=deg2rad(DIR);
                        centers=[centers,centers(1)];
                        colmap=winter;
                        alphaval=0.5;
                        
                        handle=cart2rose(centers,databis);
                        xh = get(handle,'Xdata');
                        yh = get(handle,'Ydata');
                        g=patch(xh,yh,'y');
                        par=colmap;
                        set(g,'FaceColor',par(1,:));
                        set(g,'FaceAlpha',alphaval);
                        set(g,'EdgeColor',par(1,:)/2.1);
                        set(g,'LineWidth',2);
                        
                        hold on
                        
                        handle=cart2rose(centers,data);
                        xh = get(handle,'Xdata');
                        yh = get(handle,'Ydata');
                        g=patch(xh,yh,'y');
                        par=colmap;
                        set(g,'FaceColor',par(end,:));
                        set(g,'FaceAlpha',alphaval);
                        set(g,'EdgeColor',par(end,:)/2.1);
                        set(g,'LineWidth',2);
                        
                        set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
                        sb2=subplot(666,666,666);
                        tc=tc+0.00001*rand(size(tc));
                        if z==1
                        str = sprintf(['n_',num2str(goodidx(nn)),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind))]);
                        else
                        str = sprintf(['norder ',num2str(nn),'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind))]);
                        end
                        tx=text(0.5,0.1,str); axis off
                        set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
                        poscount1=poscount1+1;
                        if poscount1==4;
                            poscount2=poscount2+1;
                            poscount1=0;
                        end
                    end
                    h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(plotnum+1),')']);
                    set(gca, 'Visible', 'off');
                    set(h, 'Visible', 'on', 'FontSize', 15);
                    set(gcf, 'Color', 'w');
                    set(gcf, 'PaperPositionMode', 'auto')
                    fname=['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(plotnum+1),')'];
                    fname=strrep(fname,'.','');
                    saveas(gcf,fname, 'jpg')
                    saveas(gcf,fname, 'epsc')
            end
            
        end
        cd(oldd)
    end
end
