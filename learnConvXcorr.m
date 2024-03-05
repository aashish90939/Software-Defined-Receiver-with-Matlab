% corrvsconv .m: ” c o r r e l a t i o n ” vs ” convolution”
x=[1 -1 2 -2 3 -3]; % de f ine sequence x [ k ]
h=[1 2 3 4 5 6 -5 -4 -3 -2 -1];
% de f ine sequence h [ k ]
yconv=conv(x,h) ;% convolve x [ k ] * h [ k ]
%ycorr=xcorr (h , f l i p l r (x ) ) % c o r r e l a t i o n of f l ipped h and x
ycorr = xcorr(x,fliplr(h));
check=max(abs([yconv,zeros(1,length(ycorr)-length(yconv))]-ycorr))

%%
k=[-1,1,1,-1,-1,1,1,1]    %8 bit
l=[-1,1,1,-1,-1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,-1,1,1,-1,-1,1,1,1]  %header+data+header =48
lengthheadernew=length(k)
lengthdatawithheader=length(l)
lengthofdataonly=lengthdatawithheader-2*lengthheadernew
corr=xcorr(k,l);
convulv= conv(flip(k),l);

figure(1)
subplot(3,1,1),stem(corr);
title('xcorr k and l')
subplot(3,1,2),stem(convulv);
title('using conv flip')
subplot(3,1,3),stem(l);
title('head8+data32+head8');
[pks,locs]=findpeaks(abs(corr),MinPeakHeight=6);


% angle(fft)