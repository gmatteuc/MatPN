function [ handle] = plot_SF_TF_distribution( sf_distr, tf_distr, SF, TF, target_group )

% [ handle] = plot_SF_TF_distribution( osi_distr, dsi_distr, SF, TF, target_group )
%
% plot SF & TF distribution
%-----------------------------------------------------------------------

X=[sf_distr,tf_distr];
ctrs{1}=SF;
ctrs{2}=TF;
[N,C]=hist3(X,ctrs);
handle=figure;
set(handle,'Position',[10,10,1500,1000]);
b=bar3c(N/sum(N(:)),1);
ylabel('SF'); xlabel('TF');
hx=get(b(1),'parent');
set(hx,'yticklabel',num2str(SF'))
set(hx,'xticklabel',num2str(TF'))
colormap('winter');
colorbar;
view(-25,20);
title(['preferred SF and TF distribution ',target_group])
fffname=['st tuning distribution ',target_group];
set(gcf, 'PaperPositionMode', 'auto')
saveas(handle,fffname, 'jpg')

end

