function [] = plot_roseplots( selectedsi, gooddir, target_group  )

% plot_tuning_curves()
%
% plot grating and plaid tuning curve of each selected neuron
%-----------------------------------------------------------------------

pars=set_pars_PN;
data_folder=pars.processed_data_folder;
load(fullfile(data_folder,'Tuning.mat'));
load(fullfile(data_folder,'Indexing.mat'));
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;

    plotnum=floor(length(selectedsi)/16);
    rplotnum=mod(length(selectedsi),16);
    
    if plotnum~=0
        for hh=1:plotnum
            f1 = figure;
            set(f1,'Position',[10,10,1500,1000]);
            poscount1=0;
            poscount2=0;
            for nnn=(1:16)+16*(hh-1)
                nn=selectedsi{nnn}(1);
                SFind=selectedsi{nnn}(2);
                TFind=selectedsi{nnn}(3);
                tc=tuning_curve( :,nn,SFind,TFind,2);
                tce=tuning_curve_error( :,nn,SFind,TFind,2);
                dsi=DSI(nn,SFind,TFind,2);
                osi=OSI(nn,SFind,TFind,2);
                mi=modulation_index(find(DIR==gooddir(nnn)),nn,SFind,TFind,1);
                tcg=tuning_curve( :,nn,SFind,TFind,1);
                tceg=tuning_curve_error( :,nn,SFind,TFind,1);
                
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
                
                hold off;
                set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
                
                sb2=subplot(666,666,666);
                
                tc=tc+0.00001*rand(size(tc));
                sessionname=[listSessions{1,M(nn,1)},'_b',num2str(M(nn,2))];
                sname=strrep(sessionname,'_',' ');
                str = sprintf(['n ',num2str(nn),' o ',num2str(nnn), ' - ', sname,'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nZc=',num2str(Zc(nnn),'%0.1f'),' Zp=',num2str(Zp(nnn),'%0.1f')]);
                tx=text(0.5,0.1,str); axis off
                set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
                poscount1=poscount1+1;
                
                if poscount1==4;
                    poscount2=poscount2+1;
                    poscount1=0;
                end
            end
            %                             h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(hh),')']);
            %                             set(gca, 'Visible', 'off');
            %                             set(h, 'Visible', 'on', 'FontSize', 15);
            %                             set(gcf, 'Color', 'w');
            set(gcf, 'PaperPositionMode', 'auto')
            fname=['selected plaid tuning curves ','(',num2str(hh),') ',target_group];
            fname=strrep(fname,'.','');
            saveas(gcf,fname, 'jpg')
            
        end
    else
    end
    close all
    
    % _-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
    
    f1 = figure;
    set(f1,'Position',[10,10,1500,1000]);
    poscount1=0;
    poscount2=0;
    for nnn=16*(plotnum)+1:16*(plotnum)+rplotnum
        nn=selectedsi{nnn}(1);
        SFind=selectedsi{nnn}(2);
        TFind=selectedsi{nnn}(3);
        tc=tuning_curve( :,nn,SFind,TFind,2);
        tce=tuning_curve_error( :,nn,SFind,TFind,2);
        dsi=DSI(nn,SFind,TFind,2);
        osi=OSI(nn,SFind,TFind,2);
        mi=modulation_index(find(DIR==gooddir(nnn)),nn,SFind,TFind,1);
        tcg=tuning_curve( :,nn,SFind,TFind,1);
        tceg=tuning_curve_error( :,nn,SFind,TFind,1);
        
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
        hold off;
        set(sb1,'Position',[.00+0.23*(poscount1),.01+0.23*(poscount2),.22,.22]);
        
        sb2=subplot(666,666,666);
        
        tc=tc+0.00001*rand(size(tc));
        sessionname=[listSessions{1,M(nn,1)},'_b',num2str(M(nn,2))];
        sname=strrep(sessionname,'_',' ');
        str = sprintf(['n ',num2str(nn),' o ',num2str(nnn), ' - ', sname,'\n\n','DSI=',num2str(dsi,'%0.1f'),' OSI=',num2str(osi,'%0.1f'),' \n\npeak FR=',num2str(max(tc),'%0.1f'),'+/-',num2str(tce(tc==max(tc)),'%0.1f'),'\n\nSF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),'\n\nZc=',num2str(Zc(nnn),'%0.1f'),' Zp=',num2str(Zp(nnn),'%0.1f')]);
        tx=text(0.5,0.1,str); axis off
        set(sb2,'Position',[.00+0.23*(poscount1)+0.11,.01+0.23*(poscount2)+0.11,0.15,0.15]);
        poscount1=poscount1+1;
        
        if poscount1==4;
            poscount2=poscount2+1;
            poscount1=0;
        end
    end
    %                     h = suptitle(['V1 tuning curves (',type,') ',stimulustype,' SF=',num2str(SF(SFind)),' TF=',num2str(TF(TFind)),' (',num2str(plotnum+1),')']);
    %                     set(gca, 'Visible', 'off');
    %                     set(h, 'Visible', 'on', 'FontSize', 15);
    %                     set(gcf, 'Color', 'w');
    set(gcf, 'PaperPositionMode', 'auto')
    fname=['selected plaid tuning curves ','(',num2str(plotnum+1),') ',target_group];
    fname=strrep(fname,'.','');
    saveas(gcf,fname, 'jpg')
    close all

end
