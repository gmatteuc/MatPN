function [] = plot_prediction_results_PN( are )

% [] = plot_prediction_results_PN(  are  )
% produce STA-based prediction population wise plots for a given area
%--------------------------------------------------------------------------

% set pars
pars = set_pars_PN();
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);
stimulustypes=pars.stimPars.stimulustype;

% load tuning analysis results and indexing
contrast=[]; % needed to prevent matlab seeing contrast as a function (bad name choiche)
load('Tuning.mat')
load('Indexing.mat')
load(['prediction_explained_variance_',are,'.mat'])
load(['selected_population_',are,'.mat'])
% sta analysis results
load(['RFs_datastructure_',are])

% select good STAs basing on contrast
contrast_th=pars.contrast_threshold;
contrasts=zeros(size(bestfr_contrast));
for kkk=1:length(bestfr_contrast)
    contrasts(kkk)=contrast(bestfr_contrast(kkk),kkk);
end
goodctr=contrasts>contrast_th;

%% plot prediction goodness distributions

% tuning explained variance  -----------------------------------------------
f1 = figure;
set(f1,'Position',[10,10,1500,1000]);
hold on
for kk=1:numel(stimulustypes)
    if kk==1
        plot(-0.5*ones(size(explvar_T(:,kk)))+1*rand(size(explvar_T(:,kk))),explvar_T(:,kk),'.g','MarkerSize',35);
        xlim([-10,20]); ylim([-0.2,1.2]); hold on; plot(0,mean(explvar_T(:,kk)),'dk','MarkerSize',10,'MarkerFaceColor',[0,0,0]);
    else
        plot(9.5*ones(size(explvar_T(:,kk)))+1*rand(size(explvar_T(:,kk))),explvar_T(:,kk),'.b','MarkerSize',35);
        xlim([-10,20]); ylim([-0.2,1.2]); hold on; plot(10,mean(explvar_T(:,kk)),'dk','MarkerSize',10,'MarkerFaceColor',[0,0,0]);
    end
    title(['tuning prediction explained variance ',are])
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    if kk==1
        te=text(-1.5,0.89*ylimit(2),['mean explained variance G = ',num2str(mean(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        tee=text(-1.5,0.85*ylimit(2),['max explained variance G = ',num2str(max(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        tee=text(-1.5,0.81*ylimit(2),['min explained variance G = ',num2str(min(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        teeee=text(-1.5,0.77*ylimit(2),['median explained variance G = ',num2str(median(explvar_T(:,kk)),'%.2f')],'FontSize',10);
    else
        te=text(8.5,0.89*ylimit(2),['mean explained variance P = ',num2str(mean(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        tee=text(8.5,0.85*ylimit(2),['max explained variance P = ',num2str(max(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        tee=text(8.5,0.81*ylimit(2),['min explained variance P = ',num2str(min(explvar_T(:,kk)),'%.2f')],'FontSize',10);
        teeee=text(8.5,0.77*ylimit(2),['median explained variance P = ',num2str(median(explvar_T(:,kk)),'%.2f')],'FontSize',10);
    end
end
fffname=['tuning prediction explained variance (tuning) ',are];
hold off
%         saveas(gcf,fffname, 'epsc')
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,fffname, 'jpg')

% tuning explained variance  ----------------------------------------------
f2 = figure;
hold on
for kk=1:numel(stimulustypes)
    set(f2,'Position',[10,10,1500,1000]);
    if kk==1
        plot(-0.5*ones(size(explvar_sD(:,kk)))+1*rand(size(explvar_sD(:,kk))),explvar_sD(:,kk),'.g','MarkerSize',35);
        xlim([-10,20]); ylim([-0.2,1.2]); hold on; plot(0,mean(explvar_sD(:,kk)),'dk','MarkerSize',10,'MarkerFaceColor',[0,0,0]);
    else
        plot(9.5*ones(size(explvar_sD(:,kk)))+1*rand(size(explvar_sD(:,kk))),explvar_sD(:,kk),'.b','MarkerSize',35);
        xlim([-10,20]); ylim([-0.2,1.2]); hold on; plot(10,mean(explvar_sD(:,kk)),'dk','MarkerSize',10,'MarkerFaceColor',[0,0,0]);
    end
    title(['dynamics prediction explained variance ',are])
    ylimit=get(gca,'ylim');
    xlimit=get(gca,'xlim');
    if kk==1
        te=text(-1.5,0.89*ylimit(2),['mean explained variance G = ',num2str(mean(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        tee=text(-1.5,0.85*ylimit(2),['max explained variance G = ',num2str(max(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        tee=text(-1.5,0.81*ylimit(2),['min explained variance G = ',num2str(min(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        teeee=text(-1.5,0.77*ylimit(2),['median explained variance G = ',num2str(median(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
    else
        te=text(8.5,0.89*ylimit(2),['mean explained variance P = ',num2str(mean(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        tee=text(8.5,0.85*ylimit(2),['max explained variance P = ',num2str(max(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        tee=text(8.5,0.81*ylimit(2),['min explained variance P = ',num2str(min(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
        teeee=text(8.5,0.77*ylimit(2),['median explained variance P = ',num2str(median(explvar_sD(:,kk)),'%.2f')],'FontSize',10);
    end
end
fffname=['tuning prediction explained variance (dynamics) ',are];
hold off
%         saveas(gcf,fffname, 'epsc')
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,fffname, 'jpg')

% pattern vs. component: bullet plot  -------------------------------------

% get data needed for neuron selection
DSIvec=zeros(length(selectedsi),1);
for nnnn=1:length(selectedsi)
    nidx=selectedsi{nnnn}(1);
    sfidx=selectedsi{nnnn}(2);
    tfidx=selectedsi{nnnn}(3);
    DSIvec(nnnn)=DSI(nidx,sfidx,tfidx,1);
end

% plot predicted partial correlation plot ...

fp=figure;
set(fp,'Position',[10,10,1500,1000]);
% plot every neuron
plot(Zc_P,Zp_P,'.','MarkerSize',45,'Color',[0,0,0]);
hold on
% plot only DS and plaid responsive ones
for nnnn=1:length(selectedsi)
    if goodprds(nnnn) && goodctr(nnnn)
        plot(Zc_P(nnnn),Zp_P(nnnn),'.','MarkerSize',25,'Color',[0,1,0]);
        ylimit=get(gca,'ylim');
        te=text(Zc_P(nnnn),Zp_P(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    elseif goodprds(nnnn) && not(goodctr(nnnn))
        plot(Zc_P(nnnn),Zp_P(nnnn),'.','MarkerSize',25,'Color',[0,0.5,0.2]);
        ylimit=get(gca,'ylim');
        te=text(Zc_P(nnnn),Zp_P(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    end
end
line([1.5 8], [0 6.5],'LineWidth',1.5,'Color','k');
line([0 6.5], [1.5 8],'LineWidth',1.5,'Color','k');
line([1.5 1.5], [-4 0],'LineWidth',1.5,'Color','k');
line([-4 0], [1.5 1.5],'LineWidth',1.5,'Color','k');
xlabel('Zc'); ylabel('Zp'); title('Pattern vs. Component - Predicted');
legend('V1 neurons - predicted','selected neurons - predicted','Location','NorthEast')
% save
set(gcf, 'PaperPositionMode', 'auto'); saveas(fp,'pattern_component_Z_predicted.jpg', 'jpg');

% plot predicted observed partial correlation plot ...

fpp=figure;
set(fpp,'Position',[10,10,1500,1000]);
% plot every neuron
plot(Zc_O,Zp_O,'.','MarkerSize',45,'Color',[0,0,0]);
hold on
% plot only DS and plaid responsive ones
for nnnn=1:length(selectedsi)
    if goodprds(nnnn) && goodctr(nnnn)
        plot(Zc_O(nnnn),Zp_O(nnnn),'.','MarkerSize',25,'Color',[0,1,0]);
        ylimit=get(gca,'ylim');
        te=text(Zc_O(nnnn),Zp_O(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    elseif goodprds(nnnn) && not(goodctr(nnnn))
        plot(Zc_O(nnnn),Zp_O(nnnn),'.','MarkerSize',25,'Color',[0,0.5,0.2]);
        ylimit=get(gca,'ylim');
        te=text(Zc_O(nnnn),Zp_O(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    end
end
line([1.5 8], [0 6.5],'LineWidth',1.5,'Color','k');
line([0 6.5], [1.5 8],'LineWidth',1.5,'Color','k');
line([1.5 1.5], [-4 0],'LineWidth',1.5,'Color','k');
line([-4 0], [1.5 1.5],'LineWidth',1.5,'Color','k');
xlabel('Zc'); ylabel('Zp'); title('Pattern vs. Component - Observed');
legend('V1 neurons - observed','selected neurons - observed','Location','NorthEast')
% save
set(gcf, 'PaperPositionMode', 'auto'); saveas(fpp,'pattern_component_Z_observed.jpg', 'jpg');

% plot observed and predicted partial correlation plot (nice only) ...
fp=figure;
set(fp,'Position',[10,10,1500,1000]);
hold on
nice_idx = [];
% plot only DS and plaid responsive ones
for nnnn=1:length(selectedsi)
    if goodprds(nnnn) && goodctr(nnnn)
        nice_idx = [nice_idx, nnnn];
        ylimit=get(gca,'ylim');
        te=text(Zc_P(nnnn),Zp_P(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',10);
    else
    end
end
plot(Zc_P(nice_idx),Zp_P(nice_idx),'.','MarkerSize',45,'Color',[0,0,1]);
for nnnn=1:length(selectedsi)
    if goodprds(nnnn) && goodctr(nnnn)
        nice_idx = [nice_idx, nnnn];
        ylimit=get(gca,'ylim');
        te=text(Zc_O(nnnn),Zp_O(nnnn)+0.07*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',10);
    else
    end
end
plot(Zc_O(nice_idx),Zp_O(nice_idx),'.','MarkerSize',45,'Color',[1,1,0]);
line([1.5 8], [0 6.5],'LineWidth',1.5,'Color','k');
line([0 6.5], [1.5 8],'LineWidth',1.5,'Color','k');
line([1.5 1.5], [-4 0],'LineWidth',1.5,'Color','k');
line([-4 0], [1.5 1.5],'LineWidth',1.5,'Color','k');
xlabel('Zc'); ylabel('Zp'); title('Pattern vs. Component - Predicted');
legend('predicted','observed','Location','NorthEast')
% save
set(gcf, 'PaperPositionMode', 'auto'); saveas(fp,'pattern_component_Zs.jpg', 'jpg');


% pattern vs. component: PI scatter  --------------------------------------

fppp=figure;
set(fppp,'Position',[10,10,1500,1000]);
hold on
plot(PI_O,PI_P,'.','MarkerSize',45,'Color',[0,0,0]);
xlim([-7,7]);
ylim([-7,7]);
xlabel('observed PI'); ylabel('predicted PI'); title('Patattern Index - observed vs. predicted ');
nice_idx=[];
for nnnn=1:length(selectedsi)
    if goodprds(nnnn) && goodctr(nnnn)
        nice_idx=[nice_idx,nnnn];
        plot(PI_O(nnnn),PI_P(nnnn),'.','MarkerSize',25,'Color',[0,0,1]);
        ylimit=get(gca,'ylim');
        te=text(PI_O(nnnn),PI_P(nnnn)+0.03*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    elseif goodprds(nnnn) && not(goodctr(nnnn))
        plot(PI_O(nnnn),PI_P(nnnn),'.','MarkerSize',25,'Color',[0.7,0.7,1]);
        ylimit=get(gca,'ylim');
        te=text(PI_O(nnnn),PI_P(nnnn)+0.03*diff(ylimit),['n = ',num2str(selectedsi{nnnn}(1))],'FontSize',15);
    end
end
te=text(0.6,0.6,['obs. - pred. PI corr (all) = ',num2str(corr(PI_O,PI_P),'%.2f')],'FontSize',15);
te=text(0.6,1.2,['obs. - pred. PI corr (nice) = ',num2str(corr(PI_O(nice_idx),PI_P(nice_idx)),'%.2f')],'FontSize',15);
legend('all V1 neurons PI','nice V1 neurons PI','Location','NorthEast')
set(gcf, 'PaperPositionMode', 'auto'); saveas(fppp,'PI_prediction_scatter.jpg', 'jpg');

end