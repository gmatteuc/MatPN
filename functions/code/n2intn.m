function [ intn, sessionname ] = n2intn( n )

% set pars and load indexing
pars=set_pars_PN;
data_folder=pars.processed_data_folder;
load(fullfile(data_folder,'Indexing.mat'));
listSessions = pars.listSessions;
% compose session name
sessionname=[listSessions{1,M(n,1)},'_b',num2str(M(n,2))];
% current session indexes
ssidx=find(M(:,1)==M(n,1));
% find first and last index of current block
blidx=intersect(find(M(:,2)==M(n,2)),ssidx);
% within block index
intn=find(blidx==n);

end

