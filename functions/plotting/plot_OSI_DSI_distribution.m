function [ handle] = plot_OSI_DSI_distribution( osi_distr, dsi_distr, target_group )

% [ handle] = plot_OSI_DSI_distribution( osi_distr, dsi_distr )
%
% plot OSI & DSI distribution
%-----------------------------------------------------------------------

binnum=8;

handle=figure;
set(handle,'Position',[10,10,1500,1000]);
subplot(1,2,1)
[v1,c1]=hist(osi_distr,linspace(0.05,0.95,binnum));
bar(c1,v1/sum(v1));
hold on
plot(c1,v1/sum(v1),'--b','LineWidth',3)
plot(median(osi_distr).*ones(1,11),0:0.05:0.5,'--k','LineWidth',2.5)
ylimit=get(gca,'ylim');
xlimit=get(gca,'xlim');
te=text(0.03*xlimit(2),0.55*ylimit(2),['median OSI = ',num2str(median(osi_distr))],'FontSize',10);
xlim([0,1])
ylim([0,0.55])
title(['OSI distribution ',target_group])
set(gca, 'box','on');
axis square
hold on 
subplot(1,2,2)
[v2,c2]=hist(dsi_distr,linspace(0.05,0.95,binnum));
hold on
bar(c2,v2/sum(v2));
plot(c2,v2/sum(v2),'--b','LineWidth',3)
plot(median(dsi_distr).*ones(1,11),0:0.05:0.5,'--k','LineWidth',2.5)
ylimit=get(gca,'ylim');
xlimit=get(gca,'xlim');
te=text(0.03*xlimit(2),0.55*ylimit(2),['median DSI = ',num2str(median(dsi_distr))],'FontSize',10);
xlim([0,1])
ylim([0,0.55])
axis square
title(['DSI distribution ',target_group])
set(gca, 'box','on');
hold on
handlename=['SI distributions ',target_group];
handlename=strrep(handlename,'.','');
set(gcf, 'PaperPositionMode', 'auto')
saveas(handle,handlename, 'jpg')

end

