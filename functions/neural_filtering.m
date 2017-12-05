function [ prate, pcount, pratetime ] = neural_filtering( filter,stimul )
%NEURAL_FILTERING 

STAf=filter-mean(mean(mean(filter)));
STAf=STAf./(max(STAf(:))-min(STAf(:)));

% convolve to compute filter response
fout=zeros(1,size(stimul,3)-size(STAf,3)+1);
for kk=1:size(stimul,3)-size(STAf,3)+1
    fout(kk)=sum(sum(sum(flipdim(STAf,3).*stimul(:,:,(1:size(STAf,3))+(kk-1)))));
end
% prediction results
prate=[fout(1)*ones(size(STAf,3)-1,1);fout'];
pcount=sum(max(prate,zeros(size(prate))));
pratetime=0.033.*(1:length(squeeze(stimul(5,5,:))))-0.033*25;

end

