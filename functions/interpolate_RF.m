function [ rf_int ] = interpolate_RF( rf )

%[ RFs ] = interpolate_RFs( RFs )
%
%Add to RFs datastructure the interpolated RFs


%--------------------------------------------------------------------------
    rf_int=zeros((size(rf,1))*30,(size(rf,2))*30,size(rf,3));
    % spline interpolation
    for jj=1:size(rf,3)

        [Xo,Yo] = meshgrid(1:1:size(rf,2), 1:1:size(rf,1));
        [Xp,Yp] = meshgrid(linspace(1,size(rf,2),(size(rf,2))*30), linspace(1,size(rf,1),(size(rf,1))*30));

            fram=squeeze(rf(:,:,jj));
            if sum(isnan(fram(:)))
            fram=zeros(size(rf,1),size(rf,2));
            else
            end
            rf_int(:,:,jj)=interp2(Xo,Yo ,fram,Xp,Yp,'spline');
    end
end

