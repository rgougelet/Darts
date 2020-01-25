% Define basic items of EEG structure.
EEG = eeg_emptyset();
length_in_sec = 20;
srate = 1000;
t =  linspace(0,length_in_sec,srate*length_in_sec);
% EEG.data  = cumsum(randn(1,length_in_sec*srate));
% EEG.data = EEG.data+std(EEG.data)*sin(2*pi*60*t);
EEG.data  = 5*sin(2*pi*60*t);
EEG.times  = t;
EEG.xmin   = EEG.times(1);
EEG.xmax   = EEG.times(end);
EEG.srate  = round(1/((EEG.xmax-EEG.xmin)/length(EEG.times))); % Rounded actual sampling rate. Note that the unit of the time must be in second.
EEG.nbchan = size(EEG.data,1);
EEG.pnts   = size(EEG.data,2);
nyq = EEG.srate/2;

close all;
figure;
hold on;
res = 5;
nfft = res*length_in_sec*srate+1;
fx = fft(EEG.data,nfft);
afx = 2*abs(res*fx/(nfft+2));
f = -nyq:(1/(res*length_in_sec)):nyq;
subplot(2,1,1)
plot(f,fftshift(afx),'b'); hold on
xlim([0 120]);
subplot(2,1,2)
plot(f,fftshift(afx),'b'); hold on
xlim([55 65]);
cl_EEG = pop_cleanline(EEG, 'bandwidth',2,'chanlist',1:EEG.nbchan ,'computepower',1,'linefreqs', 60:60:500 ,'normSpectrum',1,'p',0.01,'pad',1,'plotfigures',0,'scanforlines',0,'sigtype','Channels','tau',100,'verb',1,'winsize',1,'winstep',0.25);
fx = fft(cl_EEG.data,nfft);
afx = 2*abs(res*fx/nfft);
f = -nyq:(1/(res*length_in_sec)):nyq;
subplot(2,1,1)
title('CleanLine')
plot(f,fftshift(afx), '.-r'); hold on
xlim([0 120]);
subplot(2,1,2)
plot(f,fftshift(afx),'.-r'); hold on
xlim([55 65]);


figure;
hold on;
res = 5;
nfft = res*length_in_sec*srate+1;
fx = fft(EEG.data,nfft);
afx = 2*abs(res*fx/(nfft+2));
f = -nyq:(1/(res*length_in_sec)):nyq;
subplot(2,1,1)
plot(f,fftshift(afx), 'b'); hold on
xlim([0 120]);
subplot(2,1,2)
plot(f,fftshift(afx),'b'); hold on
xlim([55 65]);
fl_EEG = pop_eegfiltnew(EEG, 59.5,60.5, [],true);
fx = fft(fl_EEG.data,nfft);
afx = 2*abs(res*fx/nfft);
f = -nyq:(1/(res*length_in_sec)):nyq;
subplot(2,1,1)
title('Notch')

plot(f,fftshift(afx), '.-r'); hold on
xlim([0 120]);
subplot(2,1,2)
plot(f,fftshift(afx),'.-r'); hold on
xlim([55 65]);
