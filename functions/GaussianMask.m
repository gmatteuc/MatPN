function  M  = GaussianMask( a,x0,i0,sigma1,sigma2,theta,hdscalefactor )

% M  = GaussianMask( a,x0,i0,sigma1,sigma2,theta,hdscalefactor )

sx=1080/hdscalefactor;
sy=1920/hdscalefactor;
M=zeros(sx,sy);
[X,Y]=meshgrid(1:sy,1:sx);
% % MODIFICA !!!
Y=flipud(Y);
for xx=1:sx
    for yy=1:sy
        
xi=X(xx,yy);
yi=Y(xx,yy);
M(xx,yy)=a*exp(-((((xi-x0)*cos(theta) + (yi-i0)*sin(theta)).^2)/2/sigma1^2)-((-(xi-x0)*sin(theta) + (yi-i0)*cos(theta)).^2)/2/sigma2^2);

    end
end

end

