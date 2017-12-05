function [ handle] = plot_MI_distribution( goodmidx, target_group )

% [ handle] = plot_MI_distribution( osi_distr, dsi_distr, SF, TF, target_group )
%
% plot MI distribution
%-----------------------------------------------------------------------

binnum=8; %12

ff=figure;
set(ff,'Position',[10,10,1500,1000]);
[v1,c1]=hist(goodmidx,binnum);
bar(c1,v1/sum(v1));
hold on
plot(c1,v1/sum(v1),'--b','LineWidth',3)
goodmidx(isnan(goodmidx))=0;
plot(median(goodmidx).*ones(1,11),0:0.05:0.5,'--k','LineWidth',2.5)
ylimit=get(gca,'ylim');
xlimit=get(gca,'xlim');
te=text(-0.5*xlimit(2),0.40*ylimit(2),['median F1z = ',num2str(median(goodmidx))],'FontSize',10);
xlim([-3,5])
ylim([0,0.55])
title(['modulation index distribution ',target_group])
set(gca, 'box','on');
axis square
hold off
fffname=['modulation index distribution ',target_group];
set(gcf, 'PaperPositionMode', 'auto')
saveas(gcf,fffname, 'jpg')


end

