clc, clear, close all

userInput = input('Enter signal  ', 's');
userInput = upper(userInput);

switch userInput
    case 'A'
            srrclen=4;    
            beta=0.33;
            T_t=8.9e-6;
            f_s=700e3;  
            f_if=1.6e6;   
            T_s = 1/f_s;   
            tres=500;   %threshold of finding peaks in correlation
            load("mysteryA.mat");

    case 'B'
            srrclen=5;
            beta=0.4;
            T_t=7.5e-6;
            f_s=950e3; 
            f_if=1.2e6;   
            T_s = 1/f_s;
            tres=500;
            load("mysteryB.mat");

    case 'C'
            srrclen=3;
            beta=0.14;
            T_t=8.14e-6;
            f_if=2.2e6;  
            f_s=819e3;   
            T_s = 1/f_s; 
            tres=450;
            load("mysteryC.mat");
          
    otherwise
        disp('no such signals');
           
end
 
% plot the received signal
figure(1), plotspec(r, T_s)
title('Received  signal')

%sub nyquist sampling condition check
fc = abs(f_if - round(f_if/f_s)*f_s);  %aliased freq at fc centered 
% return

%% carrierrecovery with dualcostas
t= 0:T_s:length(r)*T_s;  %timevector
muc1=0.03;             %0.03
muc2=0.0003;              %0.0003
fo=fc;              %assumed at receiver
%forlpf
fl=100;                   
ff=[0 0.01 0.02 1]; fa=[1 1 0 0];   
h=firpm(fl,ff,fa);
% figure(11),freqz(h)

%costasloop
theta1=zeros(1,length(t)); theta1(1)=0;   %initialize estimate vector
theta2=zeros(1,length(t)); theta2(1)=0;

zs1=zeros( 1 , fl+1); zc1=zeros (1 , fl +1);  % buffer for lpfs
zs2=zeros( 1 , fl+1); zc2=zeros (1 , fl +1);

for k=1:length(t)-1        %zs contain past fl+1 inputs
    zs1=[zs1(2:fl+1), 2*r(k)*sin(2*pi*fo*t(k)+ theta1(k))] ;
    zc1=[zc1(2:fl+1), 2*r(k)*cos(2*pi*fo*t(k)+ theta1(k))] ;
    lpfs1=fliplr(h)*zs1' ;  lpfc1=fliplr(h)*zc1';               %new output of filters
    theta1(k+1)=theta1(k)-muc1*lpfs1*lpfc1;                     % algorithm update

    zs2=[zs2(2:fl+1), 2*r(k)*sin(2*pi*fo*t(k)+ theta1(k)+theta2(k))] ;
    zc2=[zc2(2:fl+1), 2*r(k)*cos(2*pi*fo*t(k)+ theta1(k)+theta2(k))] ;
    lpfs2=fliplr(h)*zs2' ;  lpfc2=fliplr(h)*zc2';               %new output of filters
    theta2(k+1)=theta2(k)-muc1*lpfs2*lpfc2;                     % algorithm update

    carest(k)= cos(2*pi*fo*t(k)+theta2(k)+theta1(k));

end


figure(2)
subplot(2,1,1),plot(t,theta1)
title('frequency offset')
subplot(2,1,2),plot(t,theta2)
title('phase offset')

% return
figure(3),plotspec(carest,T_s)
title('recovered carrier')


%demodulating
demodr= r'.*carest;
figure(4),plotspec(demodr,T_s) 
title('demodulated with carrier..')

%% lpf design to remove replicas 
fn=f_s/2;
low_cutoff = (1.3e4/(fn)) + 0.15;   %normalized freq with nyquist freq
high_cutoff = (1.3e4/(fn)) + 0.20;

h_demod = firpm(100,[0 low_cutoff high_cutoff 1],[1 1 0 0]);                  
baseband_signal = filter(h_demod, 1, demodr); 

figure(5),plotspec(baseband_signal,T_s)
title('name','basebandsignal ..low passed');

%% matched filter



M=T_t*f_s;      % upsampling factor
raised_cosine =srrc(srrclen, beta, M, 0);  %matched filtering
matchedfiltersignal=filter(raised_cosine, 1, baseband_signal);  

figure(6),plotspec(matchedfiltersignal,T_s);
title('Matched Filtered Signal')


% return

% neye=2;
% sl=round(srrclen*M+1);    % skip length
% % sl=0;
% M=round(M); 
% ul=floor((length(matchedfiltersignal(1:1000))-sl)/(neye*M)); 
% xeye=reshape(matchedfiltersignal(sl+1:ul*neye*M+sl),neye*M,ul);
% figure(14)
% plot(xeye)

% return

%% downsampling with timing recovery
%opm - timing recovery algorithm
tnow        = srrclen*M+1;     
tau=0;
len        = round(length(r)/M);
xs         = zeros(1,len);       %downsampled signal
i           = 0;
mut        = 0.05;               
delta_opm    = 0.7;          %time for derivative
tausave    = zeros(1,len);
while tnow < length(matchedfiltersignal)-srrclen*M
    i        = i+1;
    xs(i)   = interpsinc(matchedfiltersignal,tnow+tau,srrclen);       % interpolated at tnow +tau
    x_deltap = interpsinc(matchedfiltersignal,tnow+tau+delta_opm,srrclen); % value to right
    x_deltam = interpsinc(matchedfiltersignal,tnow+tau-delta_opm,srrclen); % value to left
    dx       = x_deltap-x_deltam;                                      % numerical derivative
    tau     = tau+mut*dx*xs(i);                                       % alg update (energy)
    tausave(i) = tau;                                                  % save for plotting
    tnow    = tnow+M;       % update
end

figure(7), 
subplot(2,1,1), plot(xs(1:i),'b.')
title('Constellation diagram');
subplot(2,1,2), plot(tausave(1:i))
title('timing offset \tau ');
ylabel('offset estimates'), xlabel('iterations')

% return
%% dd equalizer as ps

taps_eq = 7;
b_eq=zeros(taps_eq,1);
b_eq(round(taps_eq/2))=1;
% b_eq      = [0 0 1 0 0]';          % center spike initialization

mu_eq     = 0.00001;    
for i  = taps_eq+1:length(xs)
    rr = xs(i:-1:i-taps_eq+1)';                      %input
    e  = quantalph(b_eq'*rr,[-3,-1,1,3])-b_eq'*rr;   %desired..output-input
    b_eq = b_eq+mu_eq*e*rr;
end
equalizedsignal = filter(b_eq,1,xs);
figure(8), 
plot(equalizedsignal(1:i-2),'b.')
title('Equalizer Output: soft decisions');

% return
 
%% decision device / quantization 
quantized=quantalph(equalizedsignal,[-3,-1,1,3])';     

%frame sync ..
header= '0x0 This is is the Frame Header 1y1';
header = letters2pam(header);
y=conv(flip(header),quantized);

[pks,locs]=findpeaks(abs(y),'MinPeakHeight',tres); % finding peaks above threshold tres
figure(9)
subplot(3,1,1), stem(header)             % plot header
title('Header')
subplot(3,1,2), stem(quantized)             % plot data sequence
title('Data with embedded header')
subplot(3,1,3), stem(y)                % plot correlation
title('Correlation of header with data')


% figure(10)
% stem(y(1:1500))                % plot correlation
% title('Correlation of header with data')

frame_data=[];
for i=1:length(locs)-1
    data_start=locs(i)+1;                      %start index of data
    datalen=locs(i+1)-locs(i)-length(header);  %length of data
    data_end =locs(i)+datalen;                 %end index of data
    if(y(locs(i)) > 0)
        frame_data = [frame_data, quantized(data_start:data_end)];
    else
        frame_data = [frame_data,-quantized((data_start:data_end))];
    end
end

%decodepamsymbol to letters
reconstructed_message=pam2letters(frame_data)




