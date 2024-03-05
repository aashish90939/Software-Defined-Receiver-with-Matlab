%% lpf ra matchedfilter

noise=abs(matchedfiltersignal)-abs(demodr);
snr=(mean(abs(matchedfiltersignal)).^2)/(mean(noise).^2)
return
%% lpf

noise2=abs(baseband_signal)-abs(demodr);
snr2=(mean(abs(baseband_signal)).^2)/(mean(noise2).^2)

% %%
% r=signal+noise_a+noise_b
% 
% 
% a-b-c
% 
% (c-a)-(b-a)