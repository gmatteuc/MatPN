function [Syy,f] = p_BT(y,maxlag,nfft,fs,lagwindow)

% Syy = btpsd(y, maxlag, nfft, lagwindow)
%
% Blackman-Tukey spectral estimation
% Returns power spectrum and corresponding frequencies
% lagwindow: 
%            't' for a triangular window
%            'g' for a Gauss window 
%            'h' for a Hamming window 
%            'b' for a Blackman window. 

ryy = xcorr(y-mean(y),maxlag,'biased');
ryy=ryy';

if lagwindow == 't'
   w = triang(2*maxlag+1);
elseif lagwindow == 'g'
   w = gauss(2*maxlag+1);
elseif lagwindow == 'h'
   w = hamming(2*maxlag+1);
elseif lagwindow == 'b'
   w = blackman(2*maxlag+1);   
end
Syy = abs(fft(w.*ryy, nfft));
deltf=fs/nfft;
f=0:deltf:(fs-deltf);


