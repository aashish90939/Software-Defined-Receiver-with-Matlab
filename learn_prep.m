clc 
clear all
close all

M=7;         %no of sample per width
L=4;    % 
beta=0.33;   %roll of factor  , for sinc beta is 0
adv=0;    

%% creating srrc pulses and ploting spectrum ..understadning
fs=700e3;
ts=1/fs;
ps=srrc(L,beta,M,adv);
ps1=srrc(L,beta,M,1);

pls=repmat(ps,1,1000);   %create pulse series
figure(10),plotspec(pls,ts)
title('series of srrc')
% figure(99)
% periodogram(pls)
% title('psd without noise test')

%% % Generate Random Noise
noiselevel=2.5;
noise = noiselevel * randn(1,length(pls));

noisypulse=pls+noise;
figure(11),plotspec(noisypulse,ts)
title('noise plus pulse')

figure(98)
periodogram(noisypulse)
title('psd with noise test')


%%lpf to cut noise -7e4 ,7e4 
beq=[0 7e4/(fs/2) 7e4/(fs/2)+0.1 1];
lpf_coef=firpm(100,beq,[1 1 0 0]);
lpfnoisesignal=filter(lpf_coef,1,noisypulse);
figure(97),plotspec(lpfnoisesignal,ts) ,title('noisey lpfed')

figure(96)
periodogram(lpfnoisesignal)
title('psd with less noise test')

%%pulse shaping
mp=srrc(L,beta,M,adv);
r=filter(mp,1,pls);
figure(12),plotspec(r(1:100),ts)
title('pulse shaped with srrcc itself ')

%%
t=0:ts:length(pls)*ts-ts;
carrier=cos(2*pi*100000*t);
figure(13)
plot(t,carrier)

carrier2=cos(2*pi*100000*t)+1.5*cos(2*pi*100000*t);


modr=noisypulse.*carrier;
figure(14)
plotspec(modr,ts);

modr2=noisypulse.*carrier2;
figure(15)
plotspec(modr2,ts);



%%

