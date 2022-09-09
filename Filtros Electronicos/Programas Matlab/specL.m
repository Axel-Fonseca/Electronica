% This functions computes the average power spectrum with a 
% window length L = 1000
% spec input arguments:
% x = time signal to compute the average power spectrum;
% sf =  sampling frequency;
% mif = lower limit of frequency to be plotted;
% maf = upper limit of frequency to be plotted;
% flag = 1 if the figure should be displayed;
%        0 if the figure should not be displayed;
% nameplot = string with the title of the plot;


function [afr,aco]=specL(x,sf,mif,maf,L,flag,nameplot)

[nChan T] = size(x);
if nChan<T
    x = x';
end
%L=1000;
lx=length(x);sp1=zeros(L,1);
if mod(lx,L) == 0
    lx = lx;
else
    lx = lx - mod(lx,L);
end
m=0;
for i=1:L:lx
    m=m+1;
    x1=x(i:i+L-1,1);
    x1=(x1-mean(x1));%/std(x1);
    f1=fft(x1)/length(x1);
    p1=f1.*conj(f1);
    sp1=sp1+p1;
end
sp1=sp1/m;
fr=(0:L-1)/L*sf;
q1=find(fr>=mif);
q2=find(fr>=maf);
aco=sp1(q1(1):q2(1));
afr=fr(q1(1):q2(1));
cf=0.01^(1/(m-1));% sinignificance of zero coherence 
con_lim=(1-cf);
% disp('Take down the confidence limit')
% disp(con_lim)
if(flag)
    figure; plot(afr,10*log10(aco)); title([nameplot],'FontSize',14,'FontWeight','Bold')
    xlabel('Frequency (Hz)'); ylabel('PSD (dB)')
end
