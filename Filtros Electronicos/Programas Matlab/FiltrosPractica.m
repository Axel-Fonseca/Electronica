% -------------------------------------------------------------------------
%  Programa para eliminar artefactos musculares usando filtros pasa-bajas.
%  Además, es posible usar filtros pasa-banda para analizar las diferentes
%  bandas de frecuencia de los ritmos cerebrales. 
%  Los filtros pueden ser realizados en 4 diferentes aproximaciones
%  (Butterworth, Chebyshev I, Chebyshev II, Elíptico). De acuerdo a las 
%  necesidades que se tengan, será posible definir la mejor aproximación.
%  Los datos utilizados en este ejemplo fueron obtenidos con el Epoc+ de 
%  Emotiv (R), con una frecuencia de muestreo de 128 Hz.
% -------------------------------------------------------------------------

close all;
clear all;
                    %PROGRAMA REALIZADO POR GERARDO Y JORGE 
% -------------------------------------------------------------------------
%                Cargar señales, inicialización de variables
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
%eegplot(recordings,'srate',fs,'winlength',20)

% -------------------------------------------------------------------------
%            Zero-mean, Unit-variance data
% -------------------------------------------------------------------------
data2Use = (data2Use - mean(data2Use,2) * ones(1,TT))./(std(data2Use,0,2)  * ones(1,TT));    
%eegplot(data2Use,'srate',fs,'winlength',20)

                        %PROGRAMA REALIZADO POR GERARDO Y JORGE 
% -------------------------------------------------------------------------
%            PSD of Original Data
% -------------------------------------------------------------------------
[ao, bo] = specL(data2Use(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Original Channel');
figure(200); plot(ao,10*log10(bo),'b'); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;
figure(400); plot(ao,10*log10(bo),'b'); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;
figure(600); plot(ao,10*log10(bo),'b'); xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             hold on;

% -------------------------------------------------------------------------
%           Lowpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------

% ------------------------ Lowpass Butterworth ----------------------------
[bl, al]   = butter(order, fc/(fs/2));
fvtool(bl,al)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(bl,al,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af,   bf] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af,10*log10(bf),'r')
             legend('Original','Filtered')

% ------------------------ Lowpass Cheby1 ---------------------------------
[bl2, al2]   = cheby1(order,0.5,fc/(fs/2));
fvtool(bl2,al2)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(bl2,al2,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af2,   bf2] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af2,10*log10(bf2),'g')
             legend('Original','Filtered Butter','Filtered Cheby1')

                    %PROGRAMA REALIZADO POR GERARDO Y JORGE 
                    
% ------------------------ Lowpass Cheby2 ---------------------------------
[bl3, al3]   = cheby2(order,20,fc/(fs/2));
fvtool(bl3,al3)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(bl3,al3,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af3,   bf3] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af3,10*log10(bf3),'y')
             legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2')
             
% ------------------------ Lowpass Ellip  ---------------------------------
[bl4, al4]   = ellip(order,0.5,40,fc/(fs/2));
fvtool(bl4,al4)
for n=chann2Filt
    lpf_signal(n,:) = filtfilt(bl4,al4,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,lpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Cheby1')
        xlabel('Time (s)')

[af4,   bf4] = specL(lpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(200); plot(af4,10*log10(bf4),'c')
             legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2','Ellip')
             title('Lowpass Filters: Butter, Cheby1, Cheby2, Ellip')
             
% -------------------------------------------------------------------------
%           Highpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------
                        %PROGRAMA REALIZADO POR GERARDO Y JORGE 
                        
% ------------------------ Highpass Butterworth ---------------------------
[b2, a2]   = butter(order, fc/(fs/2),'high');
fvtool(b2,a2)
for n=chann2Filt
    hpf_signal(n,:) = filtfilt(b2,a2,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,hpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af20,   bf20] = specL(hpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af20,10*log10(bf20),'m')
             legend('Original','Filtered Butterworth')
             

% ------------------------ Highpass Cheby1 --------------------------------
[b3, a3]   = cheby1(order,0.5,fc/(fs/2),'high');
fvtool(b3,a3)
for n=chann2Filt
    hpf_signal(n,:) = filtfilt(b3,a3,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,hpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af21,   bf21] = specL(hpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af21,10*log10(bf21),'r')
             legend('Original','Filtered Butterworth','Filtered Cheby1')

% ------------------------ Highpass Cheby2 --------------------------------
[b4, a4]   = cheby2(order,20,fc/(fs/2),'high');
fvtool(b4,a4)
for n=chann2Filt
    hpf_signal(n,:) = filtfilt(b4,a4,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,hpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af22,   bf22] = specL(hpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af22,10*log10(bf22),'g')
             legend('Original','Filtered Butterworth','Filtered Cheby1','Filtered Cheby2')
             
% ------------------------ Highpass Ellip ---------------------------------
[b5, a5]   = ellip(order,0.5,40,fc/(fs/2),'high');
fvtool(b5,a5)
for n=chann2Filt
    hpf_signal(n,:) = filtfilt(b5,a5,data2Use(n,:));
end
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,hpf_signal(chann2Filt,:),'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af23,   bf23] = specL(hpf_signal(chann2Filt,:)',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(400); plot(af23,10*log10(bf23),'b')
             legend('Original','Filtered Butterworth','Filtered Cheby1','Filtered Cheby2','Filtered Ellip')
             title('Highpass Filters: Butter, Cheby1, Cheby2, Ellip')

% -------------------------------------------------------------------------
%           Bandpass Filters: Butter, Cheby1, Cheby2, Ellip
% -------------------------------------------------------------------------
                        %PROGRAMA REALIZADO POR GERARDO Y JORGE 

% ------------------------ Bandpass Butterworth ---------------------------
[bb, ab] = butter(order,[10/(fs/2) 30/(fs/2)]);
fvtool(bb,ab)
bpf_signal = filtfilt(bb,ab,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal,'r')
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af,   bf] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(600); plot(af,10*log10(bf),'g'); hold on;
             legend('Original','Filtered Butter')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')

% ------------------------ Bandpass Cheby1 --------------------------------
[bb2, ab2] = cheby1(order,20,[10/(fs/2) 30/(fs/2)]);
fvtool(bb2,ab2)
bpf_signal = filtfilt(bb2,ab2,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal)
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af2,   bf2] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(600); plot(af2,10*log10(bf2),'k'); hold on;
             legend('Original','Filtered Butter','Filtered Cheby1')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')

% ------------------------ Bandpass Cheby2 --------------------------------
[bb3, ab3] = cheby2(order,20,[10/(fs/2) 30/(fs/2)]);
fvtool(bb3,ab3)
bpf_signal = filtfilt(bb3,ab3,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal)
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af3,   bf3] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(600); plot(af3,10*log10(bf3),'m'); hold on;
             legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')

% ------------------------ Bandpass Ellip ---------------------------------
[bb4, ab4] = ellip(order,0.5,20,[10/(fs/2) 30/(fs/2)]);
fvtool(bb4,ab4)
bpf_signal = filtfilt(bb4,ab4,data2Use(chann2Filt,:));
figure; plot(timebase,data2Use(chann2Filt,:)); hold on; plot(timebase,bpf_signal)
        legend('Original','Filtered Butter')
        xlabel('Time (s)')

[af4,   bf4] = specL(bpf_signal',fs,1,fs/2,L,0,'Power Spectrum Filtered Channel');
figure(600); plot(af4,10*log10(bf4),'r'); hold on;
             legend('Original','Filtered Butter','Filtered Cheby1','Filtered Cheby2','Filtered Ellip')
             xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
             title('Bandpass Filters: Butter, Cheby1, Cheby2, Ellip')
