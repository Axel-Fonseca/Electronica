% -------------------------------------------------------------------------
%  Programa para eliminar artefactos musculares usando filtros pasa-bajas.
%  Adem�s, es posible usar filtros pasa-banda para analizar las diferentes
%  bandas de frecuencia de los ritmos cerebrales. 
%  Los filtros pueden ser realizados en 4 diferentes aproximaciones
%  (Butterworth, Chebyshev I, Chebyshev II, El�ptico). De acuerdo a las 
%  necesidades que se tengan, ser� posible definir la mejor aproximaci�n.
%  Los datos utilizados en este ejemplo fueron obtenidos con el Epoc+ de 
%  Emotiv (R), con una frecuencia de muestreo de 128 Hz.
% -------------------------------------------------------------------------

close all;
clear all;

% -------------------------------------------------------------------------
%                Cargar se�ales, inicializaci�n de variables
% -------------------------------------------------------------------------
%cd eeglab14_1_2b
%eeglab
%cd ..

load('data.mat')
data2Use = recordings;

% -------------------------------------------------------------------------
%              Initializing variables
% -------------------------------------------------------------------------
[N, T]     = size(data2Use);
fs         = 128;        % Sampling Frequency
fc         = 30;          % Cutoff Frequency
order      = 4; 
offset     = 1;          % Starting point
len        = 20;
interval   = offset:offset+len*fs-1;%length(data2Use);%90*fs;91*fs:(300*fs)-1;
channels   = 1:14;
L          = fs;         % No. of points to compute the spectra
data2Use   = data2Use(channels,interval);
[NN, TT]   = size(data2Use);
timebase   = (offset/fs):(1/fs):(TT+offset-1)/fs;
chann2Filt = 1;

% -------------------------------------------------------------------------
%               Plotting signals in time
% -------------------------------------------------------------------------
%eegplot(recordings,'srate',fs,'winlength',20);

% -------------------------------------------------------------------------
%            Zero-mean, Unit-variance data
% -------------------------------------------------------------------------
data2Use = (data2Use - mean(data2Use,2) * ones(1,TT))./(std(data2Use,0,2)  * ones(1,TT));    
%eegplot(data2Use,'srate',fs,'winlength',20);

% -------------------------------------------------------------------------
%            PSD of Original Data
% -------------------------------------------------------------------------
[ao, bo] = specL(data2Use(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Original Channel');
figure(200); plot(ao,10*log10(bo)); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;
[ao, bo] = specL(data2Use(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Original Channel');
figure(300); plot(ao,10*log10(bo)); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;
[ao, bo] = specL(data2Use(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Original Channel');
figure(400); plot(ao,10*log10(bo)); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;   
% -------------------------------------------------------------------------
%           Lowpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------


% ------------------------ Lowpass Butterworth ----------------------------
[bl, al]  = butter(order, fc/(fs/2));
fvtool(bl,al)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(bl,al,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af,10*log10(bf),'r')
legend('Original','Filtered Butter')

% ------------------------ Lowpass Cheby1 ---------------------------------
[b2, a2] = cheby1(order, 3, fc/(fs/2));
fvtool(b2,a2)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b2,a2,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'g')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af,10*log10(bf),'g')
legend('Original','Filtered Butter','Filtered Cheby1')

% ------------------------ Lowpass Cheby2 ---------------------------------
[b3, a3]  = cheby2(order, 10, fc/(fs/2));
fvtool(b3,a3)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b3,a3,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'m')
        legend('Original','Filtered Cheby2')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af,10*log10(bf),'m')
legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2')

% ------------------------ Lowpass Ellip  ---------------------------------
[b4,a4] = ellip(order, 3, 60, fc/(fs/2));
fvtool(b4,a4)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b4,a4,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'k')
        legend('Original','Filtered Ellip')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af,10*log10(bf),'k')
legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2','Filtered Ellip')


% -------------------------------------------------------------------------
%           Highpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------


% ------------------------ Highpass Butterworth ---------------------------
[b5, a5]  = butter(order, fc/(fs/2), 'high');
fvtool(b5,a5)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b5,a5,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af,10*log10(bf),'r')
legend('Original','Filtered Butter')

% ------------------------ Highpass Cheby1 --------------------------------
[b6, a6] = cheby1(order, 3, fc/(fs/2), 'high');
fvtool(b6,a6)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b6,a6,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'g')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af,10*log10(bf),'g')
legend('Original','Filtered Butter','Filtered Cheby1')

% ------------------------ Highpass Cheby2 --------------------------------
[b7, a7]  = cheby2(order, 10, fc/(fs/2), 'high');
fvtool(b7,a7)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b7,a7,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'m')
        legend('Original','Filtered Cheby2')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af,10*log10(bf),'m')
legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2')

% ------------------------ Highpass Ellip ---------------------------------
[b8,a8] = ellip(order, 3, 60, fc/(fs/2), 'high');
fvtool(b8,a8)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(b8,a8,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'k')
        legend('Original','Filtered Ellip')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af,10*log10(bf),'k')
legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2','Filtered Ellip')


% -------------------------------------------------------------------------
%           Bandpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------


% ------------------------ Bandpass Butterworth ---------------------------
[b9, a9] = butter(order,[10/(fs/2) 30/(fs/2)]);
fvtool(b9,a9)
bpf_signal = filtfilt(b9,a9,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal,'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af,   bf] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(300); plot(af,[10*log10(bo) 10*log10(bf)]);hold on;
             legend('Original','Filtered Butter')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')

% ------------------------ Bandpass Cheby1 --------------------------------
[b10,a10] = cheby1(order, 3, [10/(fs/2) 30/(fs/2)]);
fvtool(b10,a10)
bpf_signal = filtfilt(b10,a10,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal,'g')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af,   bf] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(300); plot(af,[10*log10(bo) 10*log10(bf)]);hold on;
legend('Original','Filtered Butter', 'Filtered Cheby1')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             
% ------------------------ Bandpass Cheby2 --------------------------------
[b11, a11]  = cheby2(order, 10, [10/(fs/2) 30/(fs/2)]);
fvtool(b11,a11)
bpf_signal = filtfilt(b11,a11,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal,'m')
        legend('Original','Filtered Cheby2')
        xlabel('Time (s)')

[af,   bf] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(300); plot(af,[10*log10(bo) 10*log10(bf)]);hold on;
legend('Original','Filtered Butter', 'Filtered Cheby1','Filtered Cheby2')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             
% ------------------------ Bandpass Ellip ---------------------------------
[b12, a12]  = ellip(order, 3, 60, [10/(fs/2) 30/(fs/2)]); 
fvtool(b12,a12)
bpf_signal = filtfilt(b12,a12,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal,'b')
        legend('Original','Filtered Ellip')
        xlabel('Time (s)')

[af,   bf] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(300); plot(af,[10*log10(bo) 10*log10(bf)]);hold on;
legend('Original','Filtered Butter', 'Filtered Cheby1','Filtered Cheby2', 'Filtered Ellip')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')