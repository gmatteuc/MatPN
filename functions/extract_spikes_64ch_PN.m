function [ Spikes_g, Spikes_m, Spikes_n, goodc, muac, noisec, clabel, i_sp, clusters, masksidx ] = extract_spikes_64ch_PN(session_name)

%[ Spikes_g, Spikes_m, Spikes_n, goodc, muac, noisec, clabel ] = extract_spikes(filename)
%
%Extract spike time samples (for good, mua and noise), good, mua and noise cluster lists and cluster
%labels fom a given .kwik file
%
%----------------------------------------------------------------------------------

filename=[session_name,'.kwik'];
t_sp_0=hdf5read(filename, '/channel_groups/0/spikes/time_samples');
i_sp_0=hdf5read(filename, '/channel_groups/0/spikes/clusters/main');
try
    t_sp_1=hdf5read(filename, '/channel_groups/1/spikes/time_samples');
    i_sp_1=hdf5read(filename, '/channel_groups/1/spikes/clusters/main');
    message='\n--- 64ch session ---\n';
    fprintf(message);
    is64=1;
catch
    t_sp_1=[];
    i_sp_1=[];
    message='\n--- 32ch session ---\n';
    fprintf(message);
    is64=0;
end
clusters_0=unique(i_sp_0);
clusters_1=unique(i_sp_1);

%% recover cluster identity labels

clabel_0=NaN(1,length(clusters_0) );
for i=1:length(clusters_0)
    
    % find current cluster info
    adress=['/channel_groups/0/clusters/main/',num2str(clusters_0(i))];
    info1=h5info(filename,adress );
    
    % find index of cluster type label
    anames=cell(1,size(info1.Attributes,1));
    for j=1:size(info1.Attributes,1)
        anames{j}=info1.Attributes(j).Name;
        if strcmp(anames{j},'cluster_group')
            clu_info=j;
        else
        end
    end
    
    % strore cluster label for current cluster (0=Noise, 1=Mua, 2=Good)
    clabel_0(i)=info1.Attributes(clu_info).Value;
    fprintf(['Cluster ',num2str(clusters_0(i)),' label read: ',num2str(clabel_0(i)),'\n'])
end

clabel_1=NaN(1,length(clusters_1) );
for i=1:length(clusters_1)
    
    % find current cluster info
    adress=['/channel_groups/1/clusters/main/',num2str(clusters_1(i))];
    info1=h5info(filename,adress );
    
    % find index of cluster type label
    anames=cell(1,size(info1.Attributes,1));
    for j=1:size(info1.Attributes,1)
        anames{j}=info1.Attributes(j).Name;
        if strcmp(anames{j},'cluster_group')
            clu_info=j;
        else
        end
    end
    
    % strore cluster label for current cluster (0=Noise, 1=Mua, 2=Good)
    clabel_1(i)=info1.Attributes(clu_info).Value;
    fprintf(['Cluster ',num2str(clusters_1(i)),' label read: ',num2str(clabel_1(i)),'\n'])
end

%% get spikes corresponding to selected custers

% collate shanks
clabel=cat(1,clabel_0',clabel_1');
t_sp=cat(1,t_sp_0,t_sp_1);
shank2offset=length(clusters_0)+1;
clusters=cat(1,clusters_0,clusters_1+shank2offset);
i_sp=cat(1,i_sp_0,i_sp_1+shank2offset);

% find good clusters
goodc=clusters(clabel==2);

sf= 2.4414e+04;
Spikes_g=cell(1,numel(goodc));
for k=1:numel(goodc)
    Spikes_g{k}=t_sp(i_sp==goodc(k));
    % convert into milliseconds
    Spikes_g{k}=(double(Spikes_g{k})./sf).*1000;
end

% find mua clusters
muac=clusters(clabel==1);

sf= 2.4414e+04;
Spikes_m=cell(1,numel(muac));
for k=1:numel(muac)
    Spikes_m{k}=t_sp(i_sp==muac(k));
    % convert into milliseconds
    Spikes_m{k}=(double(Spikes_m{k})./sf).*1000;
end

% find noise clusters
noisec=clusters(clabel==0);

sf= 2.4414e+04;
Spikes_n=cell(1,numel(noisec));
for k=1:numel(noisec)
    Spikes_n{k}=t_sp(i_sp==noisec(k));
    % convert into milliseconds
    Spikes_n{k}=(double(Spikes_n{k})./sf).*1000;
end


filenamebis=[session_name,'.kwx'];
mask=hdf5read(filenamebis, '/channel_groups/0/features_masks');
masks=squeeze(mask(2,1:3:end,:));
clear mask
[~,masksidx_0]=max(masks+0.001*rand(size(masks),'single'));
clear masks
if is64
    mask=hdf5read(filenamebis, '/channel_groups/1/features_masks');
    masks=squeeze(mask(2,1:3:end,:));
    clear mask
    [~,masksidx_1]=max(masks+0.001*rand(size(masks),'single'));
    clear masks
else
    masksidx_1=[];
end
masksidx=cat(1,masksidx_0',masksidx_1');

end