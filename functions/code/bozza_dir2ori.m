
        x = (0:1:12)';
        y1 = gaussmf(x, [1.5 10]);
        y1_ori=(y1(1:length(y1)/2)+y1(1+length(y1)/2:length(y1)))./2;
        figure; plot(x,y1); hold on; plot(x(1:6),y1_ori,'r');