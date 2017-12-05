function [ bestfr, goodneu, ngoodneu, nneu ] = select_best_frames( RFs )

%[ bestfr, goodneu, ngoodneu, nneu ] = select_best_frames( RFs )
%
%Select good neurons (goodneu), their number (ngoodneu), best frames (bestfr)
%and total neuron number (bestfr)

%--------------------------------------------------------------------------

maxv=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
minv=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
contr=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
integr=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
vcontr=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
vintegr=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
stdv=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
vstd=zeros(size(RFs.DZstafr,3),size(RFs.DZstafr,4));
imaxcontr=zeros(1,size(RFs.DZstafr,4));
maxcontr=zeros(1,size(RFs.DZstafr,4));
imaxintegr=zeros(1,size(RFs.DZstafr,4));
maxintegr=zeros(1,size(RFs.DZstafr,4));
imaxstd=zeros(1,size(RFs.DZstafr,4));
maxstd=zeros(1,size(RFs.DZstafr,4));

% loop over neurons
for nn=1:size(RFs.DZstafr,4);
    % find index corresponding to current neuron number (nn)
    kk=find(RFs.Dneuronum==nn);
    if isempty(kk)
    else
        % loop over frames
        for jj=1:size(RFs.DZstafr,3);
            fram=RFs.DZstafr(:,:,jj,kk);
            % maxv(jj,kk)=max(fram(:)-mean(fram(:)));
            % minv(jj,kk)=min(fram(:)-mean(fram(:)));
            stdv(jj,nn)=(std(fram(:)))^2;
            maxv(jj,nn)=max(fram(:));
            minv(jj,nn)=abs(min(fram(:)));
            contr(jj,nn)=maxv(jj,nn)+minv(jj,nn);
            integr(jj,nn)=mean(abs(fram(:)));
        end
        vcontr(:,nn)=squeeze(contr(:,nn));
        imaxcontr(nn)=find(vcontr(:,nn)==max(vcontr(:,nn)));
        maxcontr(nn)=vcontr(imaxcontr(nn),nn);
        
        vstd(:,nn)=squeeze(stdv(:,nn));
        imaxstd(nn)=find(vstd(:,nn)==max(vstd(:,nn)));
        maxstd(nn)=vcontr(imaxstd(nn),nn);
        
        vintegr(:,nn)=squeeze(integr(:,nn));
        imaxintegr(nn)=find(vintegr(:,nn)==max(vintegr(:,nn)));
        maxintegr(nn)=vintegr(imaxintegr(nn),nn);
%         imaxv(nn)=find(maxv(:,nn)==max(maxv(:,nn)));
    end
end

ngoodneu=sum(maxcontr>8);
goodneu=find(maxcontr>8);
nneu=numel(maxcontr);
bestfr=imaxcontr;

% plot selection results
 figure; 
 subplot(3,1,1)
 plot(imaxcontr,'-*b','LineWidth',2);  hold on; plot(5*maxcontr./max(maxcontr),'--b','LineWidth',1); legend('contrast max frame','contrast max nvalue');
 subplot(3,1,2)
 plot(imaxintegr,'-*k','LineWidth',2); hold on; plot(5*maxintegr./max(maxintegr),'--k','LineWidth',1); legend('luminance','luminance max nvalue');
 subplot(3,1,3)
 plot(imaxstd,'-*g','LineWidth',2); hold on; plot(5*maxstd./max(maxstd),'--g','LineWidth',1); legend('variance','variance max nvalue');
 

 
%   title('max value frame');
end

