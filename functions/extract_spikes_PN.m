function [ Spikes_g, Spikes_m, Spikes_n, goodc, muac, noisec, clabel ] = extract_spikes_PN(filename)

%[ Spikes_g, Spikes_m, Spikes_n, goodc, muac, noisec, clabel ] = extract_spikes(filename)
%
%Extract spike time samples (for good, mua and noise), good, mua and noise cluster lists and cluster
%labels fom a given .kwik file
%
%----------------------------------------------------------------------------------

t_sp=hdf5read(filename, '/channel_groups/0/spikes/time_samples');
i_sp=hdf5read(filename, '/channel_groups/0/spikes/clusters/main');

clusters=unique(i_sp);

%% recover cluster identity labels

clabel=NaN(1,length(clusters) );
for i=1:length(clusters)
    
    % find current cluster info
    adress=['/channel_groups/0/clusters/main/',num2str(clusters(i))];
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
    clabel(i)=info1.Attributes(clu_info).Value;
    fprintf(['Cluster ',num2str(clusters(i)),' label read: ',num2str(clabel(i)),'\n'])
end

%% get spikes corresponding to selected custers

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

end