function [ sp, totspikes_after ] = get_noise_response_PN( n, stimlength, boolgpsth )

% set pars
pars = set_pars_PN();
listSessions = pars.listSessions;
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);

% load tuning analysis results and indexing
load('Tuning.mat')
load('Indexing.mat')

% get noise response path for current neuron
session=listSessions{1,M(n,1)};
block=M(n,2);
neuronum=M(n,10);
session_name=[session,'_b',num2str(block)];
response_data_folder=[data_folder,filesep,'noise_response_',session_name];
file_name=['PSTH_RASTER_',num2str(neuronum)];
neuron_path=[response_data_folder,'/',file_name];

% load neuron's frame-binned PSTHs
load(neuron_path)
sbyf=PsthAndRaster_M.MySpikes.Frames;
clear PsthAndRaster_M

% reshape spike data into a single numspikes per frame vector aligned with stimulus frames timig
sp=zeros(1,size(sbyf{1},2)*size(sbyf,1)-40);
for mm=1:size(sbyf,1)
    for nn=1:size(sbyf{1},2)-1
        if not(isempty(sbyf{mm}))
            sp(nn+(mm-1)*stimlength)=numel(sbyf{mm}{nn});
        else
        end
    end
end
sp = sp';

% generate "global PSTH" cumulating spikes of all trials
hsp=zeros(1,size(sbyf{1},2));
for mm=1:size(sbyf,1)
    for nn=1:size(sbyf{1},2)
        if not(isempty(sbyf{mm}))
            hsp(nn)=hsp(nn)+numel(sbyf{mm}{nn});
        else
        end
    end
end

% count total number of spikes
totspikes=sum(hsp);
totspikes_after=sum(sp);

if boolgpsth
    
    % print info on current neuron being analyzed
    message=['\nAnalyzing neuron ',num2str(n),' session ',session_name,' (goodn = ',num2str(neu),') ...\n'];
    fprintf(message)
    message=['total number of spikes = ',num2str(totspikes),'\n'];
    fprintf(message)
    message=['total number of spikes after initial bump elimination = ',num2str(totspikes_after),'\n'];
    fprintf(message)
    message=[num2str(32*18*50),' spikes ideally required for a good STC reconstruction\n'];
    fprintf(message)
    
    f1=figure;
    set(f1,'Position',[10,10,1500,1000]);
    plot(hsp,'color',[0.1,0.3,0.9],'LineWidth',1.5)
    hold on
    plot(sgolayfilt(hsp, 2, 27),'k','LineWidth',3)
    legend('Mean Noise PSTH','Mean Noise PSTH - Smoothed');
    ylabel('spikes per frame');
    xlabel('frame number');
    title(['Neuron ',num2str(neuronum),': ',num2str(totspikes),' spikes before, ',num2str(totspikes_after),' spikes after'])
    hold off
    filename1='global PSTH.jpg';
    set(gcf, 'PaperPositionMode', 'auto');
    saveas(f1, filename1);
    close all
    
else
end
end

