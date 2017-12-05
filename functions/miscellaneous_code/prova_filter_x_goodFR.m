function [bool_goodFR, idx_goodFR] = filter_x_goodFR(S)

%[bool_goodFR, idx_goodFR] = filter_x_goodFR(S)
%
% Returns bool_goodFR in which ther is a 1 for neurons in a given block having a decent FR
% and 0 otherwise and idx_goodFR with the indeces of neurons with decent FR

%--------------------------------------------------------------------------

S=load('SPIKEMAT_12_01_2016_b6.mat');


idx_goodFR=[];
bool_goodFR = zeros(size(M,1),1);
gr_matrix=S.SPIKEmean(100:207,:);
mean_FR=mean(gr_matrix,1);
mean_FR(mean_FR>20)=NaN;
mean_FR(mean_FR<1)=NaN;


std_FR=std(gr_matrix,0,1);
diff_FR=max(gr_matrix,1)-min(gr_matrix,1);
diff_FR=diff_FR(1,:);
modulation_FR=max(gr_matrix,1)./min(gr_matrix,1);%diff_FR./mean_FR;
modulation_FR=modulation_FR(1,:);
figure; plot(mean_FR,'r'); hold on; plot(std_FR,'*b'); hold on; plot(diff_FR,'--k');
figure; plot(modulation_FR,'--k')
figure; plot(mean_FR,'-*b') 
figure; imagesc(gr_matrix);colormap('paruly'); caxis([0,30]); colorbar;
sm=nanstd(mean_FR)
mm=nanmean(mean_FR)

%     for i=1:size(S.SPIKEmean,2);
%         %     if decision
%         %   bool_area(i)=1;
%         %   idx_area=[idx_area,i];
%         %     else
%         %     end
%     end


end

