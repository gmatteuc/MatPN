% ------------------------- DO NOISE ANALYSIS -------------------------

clear
close all
clc

% set analysis parameters
ridgeparam = 3; % Tychonov regularization parameter
maxlag=10; % max pre-spike lag
CriticalSize = 1e8;
np=1;
stimlength=1818;
stimwidth=32;
stimheight=18;

% set pars
pars = set_pars_PN();
listSessions = pars.listSessions;
SF=pars.stimPars.SF;
TF=pars.stimPars.TF;
DIR=pars.stimPars.DIR;
types=pars.neuPars.ctypes;
code_folder=pars.code_folder;
addpath(code_folder);
data_folder=pars.processed_data_folder;
addpath(data_folder);
cd(data_folder);

% load tuning analysis results and indexing
load('Tuning.mat')
load('Indexing.mat')

% get noise movie
Stim=get_noise_movie_PN();

% loop over areas
target_areas={'V1','LM','AL','RL'};

for target_area_idx=1:length(target_areas)
    
    % load neuron selection results
    load(['selected_population_',target_areas{target_area_idx},'.mat'])
    
    % create current area's output folder
    out_folder=[target_areas{target_area_idx},'_STA_results'];
    mkdir(out_folder)
    old_folder=cd(out_folder);
    
    for neu=1:length(selectedsi)
        
        n=selectedsi{neu}(1);
        boolgpsth=0;
        % get current neuron response
        [ sp, totspikes_after ] = get_noise_response_PN( n, stimlength, boolgpsth );
        
        % create output folder and cd in
        fnam=['STA_n',num2str(n),'_results'];
        mkdir(fnam);
        old_folder_bis=cd(fnam);
        
        Dstafr=zeros(18,32,maxlag);
        Dwstafr=zeros(18,32,maxlag);
        DZwstafr=zeros(18,32,maxlag);
        DZstafr=zeros(18,32,maxlag);
        Drawcov=zeros(18*32,18*32);
        
        % loop pre-spike lags
        for nlag = 1:maxlag
            
            % apply lag and zero pad
            sp = circshift(sp,-nlag);
            sp(end-nlag+1:end)=zeros(1,nlag);
            % standardize stimulus
            Stim=(Stim-mean(Stim(:)))*(std(Stim(:))).^-1;
            
            if not(exist('rawcov','var')==1)
                % compute stimulus ensemble covariance matrix
                [~,~,~,rawcov] = simpleSTC( Stim,sp, np, CriticalSize);
                % perform Tikhonov regularized inversion of the covariance matrix
                covInv = inv(rawcov+ridgeparam*eye(size(rawcov)));
                covInvsqrt = sqrtm(covInv);
            end
            
            % perform spike triggered average
            [sta] = simpleSTA(Stim, sp, np,CriticalSize);
            
            % whiten STA
            wsta = covInv*sta; %#ok<*MINV>
            
            message=['Neuron (',target_areas{target_area_idx},')',num2str(n),' STA computed\n'];
            
            %% permutation test
            
            % initialize permuation test variables
            nperm=50;
            staT=zeros(size(sta,1),nperm);
            wstaT=zeros(size(wsta,1),nperm);
            
            % loop over permutations
            for j=1:nperm
                
                % randomly reshuffle spike timestamps in time
                spP=sp(randperm(length(sp)));
                
                % perform ST analysis again
                spP = circshift(spP,-nlag);
                spP(end-nlag+1:end)=zeros(1,nlag);
                [staP] = simpleSTA( Stim, spP, np, CriticalSize);
                
                % whiten permuted STA
                wstaP = covInv*staP;
                staT(:,j)=staP;
                wstaT(:,j)=wstaP;
                
                message=['permutation number ',num2str(j),' completed\n'];
                fprintf(message)
                
            end
            
            % compute permutation means
            musta=mean(staT,2);
            muwsta=mean(wstaT,2);
            
            % compute permuatation standard deviations
            sigsta=std(staT,0,2);
            sigwsta=std(wstaT,0,2);
            
            %% permutation test results visualization
            
            % reshape results
            stafr=reshape(sta,stimheight,stimwidth);
            wstafr=reshape(wsta,stimheight,stimwidth);
            mustafr=reshape(musta,stimheight,stimwidth);
            muwstafr=reshape(muwsta,stimheight,stimwidth);
            
            % compute Z-scored STA, whitened STA and STC
            Zstafr = (stafr-mustafr)./mean(sigsta(:));
            Zwstafr = (wstafr-muwstafr)./mean(sigwsta(:));
            
            % visualize Z-scored STA, whitened STA and STC
            Z1=figure; set(Z1,'Position',[10,10,1500,1000]); set(gca,'dataAspectRatio',[1 1 1]);
            imagesc(Zstafr,[-5,5]); colormap('gray'); colorbar; title('Z-scored STA filter'); set(gcf, 'PaperPositionMode', 'auto');
            Z2=figure; set(Z2,'Position',[10,10,1500,1000]); set(gca,'dataAspectRatio',[1 1 1]);
            imagesc(Zwstafr,[-5,5]); colormap('gray'); colorbar; title('Z-scored whitened STA filter'); set(gcf, 'PaperPositionMode', 'auto');
            filenames_lag={['Zsta_at_lag',num2str(nlag),'.jpeg'],['Zwsta_at_lag',num2str(nlag),'.jpeg'],['Zstc1_at_lag',num2str(nlag),'.jpeg'],['Zstc2_at_lag',num2str(nlag),'.jpeg'],['Zeig_at_lag',num2str(nlag),'.jpeg']};
            filename1=['Z_scored_STA_at_lag',num2str(nlag),'.jpeg'];
            filename2=['Z_scored_wSTA_at_lag',num2str(nlag),'.jpeg'];
            saveas(Z1, filename1);
            saveas(Z2, filename2);
            close all
            
            %% compute significance maps:
            
            N=3; % confidence level for the map in sigmas
            
            % significativity for wSTA
            signifmapwsta=zeros(size(wsta));
            for i=1:size(wsta,1)
                if abs(wsta(i)-muwsta(i))>N*sigwsta(i)
                    signifmapwsta(i)=abs(abs(wsta(i))-abs(muwsta(i)));
                end
            end
            signifmapwstafr=reshape(signifmapwsta,18,32);
            
            % significativity for STA
            signifmapsta=zeros(size(sta));
            for i=1:size(sta,1)
                if abs(sta(i)-musta(i))>N*sigsta(i)
                    signifmapsta(i)=abs(abs(sta(i))-abs(musta(i)));
                end
            end
            signifmapstafr=reshape(signifmapsta,18,32);
            
            %% plot results summarizing figure
            
            imag=figure;
            set(imag,'Position',[10,10,1500,1000]);
            
            % raw results
            subplot(3,2,1)
            imagesc(wstafr,[min(wstafr(:)),max(wstafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('reconstructed wSTA filter', 'FontSize', 7);
            subplot(3,2,2)
            imagesc(stafr,[min(stafr(:)),max(stafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('reconstructed STA filter', 'FontSize', 7);
            % signifmaps
            subplot(3,2,3)
            imagesc(signifmapwstafr,[min(signifmapwstafr(:)),max(signifmapwstafr(:))+0.0001]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('signifmap wSTA filter', 'FontSize', 7);
            subplot(3,2,4)
            imagesc(signifmapstafr,[min(signifmapstafr(:)),max(signifmapstafr(:))+0.0001]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('signifmap STA filter', 'FontSize', 7);
            % raw permutation means
            subplot(3,2,5)
            imagesc(muwstafr,[min(muwstafr(:)),max(muwstafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('wSTA permutation mean', 'FontSize', 7);
            subplot(3,2,6)
            imagesc(mustafr,[min(mustafr(:)),max(mustafr(:))]); colormap('paruly'); set(gca,'dataAspectRatio',[1 1 1]);
            title('STA permutation mean', 'FontSize', 7);
            filenam=['ST analysis result at lag',num2str(nlag),'.jpeg'];
            saveas(imag, filenam);
            close all
            
            %% save results
            
            % store results at current lag
            Dstafr(:,:,nlag)=stafr;
            Dwstafr(:,:,nlag)=wstafr;
            DZwstafr(:,:,nlag)=Zwstafr;
            DZstafr(:,:,nlag)=Zstafr;
            Dtotspikes=totspikes_after;
            
            if nlag==maxlag
                % store also stimulus ensamble covariance matrix
                Drawcov(:,:)=rawcov;
                % save results
                save(['STA_results','_n',num2str(n),'.mat'],'Dstafr','Dwstafr','DZwstafr','DZstafr','Drawcov','Dtotspikes','neu','n')
                % cd back to original folder
                cd(old_folder_bis)
                % clean rawcov variable
                clear rawcov
            else
            end
            
            close all
            
        end
        
    end
    
    cd(old_folder)
    
end