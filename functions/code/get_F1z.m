function [ F1z, pow, f, N_f ] = get_F1z( psth, TF, T_s )

    % demean signal
    ds=psth;
    baseline=mean(ds);
    ds = ds-baseline;
    
    % compute power spectrum
    N=length(psth);
    N_f=2^8;
    f_s=1/T_s;
    maxlag=fix(N/3);
    lagwindow='t';
    [pow,f]=p_BT(ds, maxlag, N_f,f_s,lagwindow);
    
    % find stimulus frequency vector index
    fidx=find(f<=TF);
    fidx=fidx(end);
    spect=pow(fidx); sigspect=std(pow); meanspect=mean(pow);
    F1z = (spect-meanspect)/sigspect;

end

