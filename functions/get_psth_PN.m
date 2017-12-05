function [ psth, raster, hedges ] = get_psth_PN(  n, SF, TF, DIR, stimulustype  )

% set pars
pars = set_pars_PN();
listSessions = pars.listSessions;
% load indexing
load('Indexing.mat')
% select bitcodes for the tuning curve
selected_idx = get_indexes_PN( SF, TF, DIR, stimulustype );
sessionname=[listSessions{1,M(n,1)},'_b',num2str(M(n,2))];
% get timestamps data
finame=['SPIKEMAT_',sessionname];
S=load(finame);
% current session indexes
ssidx=find(M(:,1)==M(n,1));
% find first and last index of current block
blidx=intersect(find(M(:,2)==M(n,2)),ssidx);
% within block index
nnum=find(blidx==n);

% build psth -------------------------------------------------------------


% collect spike times
S_ts=[];
for i=1:size(S.SPIKEtimes,3)
    S_ts=[S_ts;S.SPIKEtimes{ selected_idx, nnum, i}];
end
S_ts=S_ts(S_ts<=1.2);
S_ts=S_ts(S_ts>=0)-0.2;
T_s=0.010;
hedges=0:T_s:1;
% produce psth
[psth]=hist(S_ts,hedges);

% build raster -------------------------------------------------------------

% collect spike times
spp=cell(1,20);
for trial=1:20
    if trial<=size(S.SPIKEtimes( selected_idx, nnum, :),3)      
        sp=S.SPIKEtimes{ selected_idx, nnum, trial}-0.2;
    else
        sp=[];
    end
    spp{trial}=sp;
end
raster=spp;

end

