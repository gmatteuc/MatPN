function [] = plot_tuning_curves( selectedsi )

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

for nnn=1:length(selectedsi)
    
    % retrive neuron number and preferred TF/SF
    nn=selectedsi{nnn}(1);
    SFind=selectedsi{nnn}(2);
    TFind=selectedsi{nnn}(3);
    % do temporary assignations
    Stuning_curve=tuning_curve(:,nn,SFind,TFind,1);
    Sc_tuning_curve=c_tuning_curve(:,nn,SFind,TFind,1);
    Sr_tuning_curve=r_tuning_curve(:,nn,SFind,TFind,1);
    St_tuning_matrix=t_tuning_matrix(:,:,nn,SFind,TFind,1);
    SSF=SF(SFind);
    STF=TF(TFind);
    SOSI=OSI(nn,SFind,TFind,1);
    SZp=Zp(nn,SFind,TFind);
    SDSI=DSI(nn,SFind,TFind,1);
    SZc=Zc(nn,SFind,TFind);
    SPI=PI(nn,SFind,TFind);
    SMI=modulation_index(find(max(Stuning_curve)),nn,SFind,TFind,1);
    Sm_tcorr=m_tcorr(nn,SFind,TFind,1);
    
    % to recycle the figure
    if not(isempty(get(0,'children')))
        f1=get(0,'children');
    else
        f1 = figure;
    end
    set(f1,'Position',[10,10,1500,1000]);
    
    % plot tuning curves
    sb1=subplot(666,666,666);
    par=paruly;
    coffset=10;
    plot(Stuning_curve,'Color',par(end-coffset,:),'LineWidth',3.5); hold on; plot(Sc_tuning_curve,'Color',par(coffset,:),'LineWidth',3.5); plot(Sr_tuning_curve./20,'k','LineWidth',2.5); plot(ones(size(Sr_tuning_curve)),'--k');
    xlim([1,12]); ylim([0,20]); legend('original tc','corrected tc','% responsive tc'); ylabel('firing rate (Hz)'); xlabel('direction (deg)');
    set(gca,'XTickLabel',num2str(round(linspace(0,330,12)')));
    sessionname=[listSessions{1,M(nn,1)},'_b',num2str(M(nn,2))];
    sname=strrep(sessionname,'_',' ');
    title(['n ',num2str(nn),' goodn ',num2str(nnn), ' - ', sname, ' SF=',num2str(SSF),' TF=',num2str(STF),' - grating tc']);
    text(10.4,17,['MI=',num2str(SMI,'%.2f')],'FontSize',10);
    text(10.4,15,['DSI=',num2str(SDSI,'%.2f')],'FontSize',10);
    text(10.4,16,['OSI=',num2str(SOSI,'%.2f')],'FontSize',10);
    text(9,17,['PI=',num2str(SPI,'%.2f')],'FontSize',10);
    text(9,15,['Zp=',num2str(SZp,'%.2f')],'FontSize',10);
    text(9,16,['Zc=',num2str(SZc,'%.2f')],'FontSize',10);
    text(9,14,['mtcorr=',num2str(Sm_tcorr,'%.2f')],'FontSize',10);
    set(sb1,'Position',[.05,.13,.41,.75]);
    
    % plot trial/direction spike count matrix ("t_tuning_matrix")
    sb2=subplot(666,666,666);
    imagesc(flipud(St_tuning_matrix')); colormap('paruly');
    ylabel('trial number'); xlabel('direction (deg)');
    set(gca,'XTickLabel',num2str([30;90;150;210;270;330]));
    set(gca,'YTickLabel',num2str(flipud((2:2:20)')-1));
    hc=colorbar; ylabel(hc,'spike count');
    title(['spike count across trials and directions']);
    set(sb2,'Position',[.51,.13,.41,.75]);
    
    % save figure
    set(gcf, 'PaperPositionMode', 'auto')
    saveas(gcf,strrep(['tuning goodn_',num2str(nnn),' n_',num2str(nn),'_SF',num2str(SSF),'_TF',num2str(STF),'_grating'],'.',''),'jpeg')
    clf(f1)
    
    % ------------------------------------------------------------------- %
    
    % do temporary assignations
    Stuning_curve=tuning_curve(:,nn,SFind,TFind,2);
    Sc_tuning_curve=c_tuning_curve(:,nn,SFind,TFind,2);
    Sr_tuning_curve=r_tuning_curve(:,nn,SFind,TFind,2);
    St_tuning_matrix=t_tuning_matrix(:,:,nn,SFind,TFind,2);
    SSF=SF(SFind);
    STF=TF(TFind);
    SOSI=OSI(nn,SFind,TFind,2);
    SZp=Zp(nn,SFind,TFind);
    SDSI=DSI(nn,SFind,TFind,2);
    SZc=Zc(nn,SFind,TFind);
    SPI=PI(nn,SFind,TFind);
    SMI=modulation_index(find(max(Stuning_curve)),nn,SFind,TFind,2);
    Sm_tcorr=m_tcorr(nn,SFind,TFind,2);
    
    % to recycle the figure
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
    plot(Stuning_curve,'Color',par(end-coffset,:),'LineWidth',3.5); hold on; plot(Sc_tuning_curve,'Color',par(coffset,:),'LineWidth',3.5); plot(Sr_tuning_curve./20,'k','LineWidth',2.5); plot(ones(size(Sr_tuning_curve)),'--k');
    xlim([1,12]); ylim([0,20]); legend('original tc','corrected tc','% responsive tc'); ylabel('firing rate (Hz)'); xlabel('direction (deg)');
    set(gca,'XTickLabel',num2str(round(linspace(0,330,12)')));
    sessionname=[listSessions{1,M(nn,1)},'_b',num2str(M(nn,2))];
    sname=strrep(sessionname,'_',' ');
    title(['n ',num2str(nn),' goodn ',num2str(nnn), ' - ', sname, ' SF=',num2str(SSF),' TF=',num2str(STF),' - plaid tc']);
    text(10.4,17,['MI=',num2str(SMI,'%.2f')],'FontSize',10);
    text(10.4,15,['DSI=',num2str(SDSI,'%.2f')],'FontSize',10);
    text(10.4,16,['OSI=',num2str(SOSI,'%.2f')],'FontSize',10);
    text(9,17,['PI=',num2str(SPI,'%.2f')],'FontSize',10);
    text(9,15,['Zp=',num2str(SZp,'%.2f')],'FontSize',10);
    text(9,16,['Zc=',num2str(SZc,'%.2f')],'FontSize',10);
    text(9,14,['mtcorr=',num2str(Sm_tcorr,'%.2f')],'FontSize',10);
    set(sb1,'Position',[.05,.13,.41,.75]);
    
    % plot trial/direction spike count matrix ("t_tuning_matrix")
    sb2=subplot(666,666,666);
    imagesc(flipud(St_tuning_matrix')); colormap('winter');
    ylabel('trial number'); xlabel('direction (deg)');
    set(gca,'XTickLabel',num2str([30;90;150;210;270;330]));
    set(gca,'YTickLabel',num2str(flipud((2:2:20)')-1));
    hc=colorbar; ylabel(hc,'spike count');
    title(['spike count across trials and directions']);
    set(sb2,'Position',[.51,.13,.41,.75]);
    
    % save figure
    set(gcf, 'PaperPositionMode', 'auto')
    saveas(gcf,strrep(['tuning goodn_',num2str(nnn),' n_',num2str(nn),'_SF',num2str(SSF),'_TF',num2str(STF),'_plaid'],'.',''),'jpeg')
    clf(f1)
    
end

end

