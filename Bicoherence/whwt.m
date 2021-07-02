% % harmonic wavelet of windowed rectangular shape in frequency domain
function [x_hwt,W,Fs,Ts]=whwt(x,f1,f2,fs,N)
% extract desired frequency band component of signal
% x is the signal to be decomposed,size of N*1
% fs is the sampling frequency
% f1,f2: low- and high-est frequency
% fs: sampling frequency
% N:number of FFT points, if absent,will be defined in code.
% x_hwt:the wavelet coefficients,every line correspond to a
% subband
% W: wavelet filter used in frequency domain
% Fs: frequency spread
% Ts: time spread
% prepare the signal
xsize=size(x);
if xsize(2)==1
    x=x';
end
L=length(x);

df=fs/L;
if nargin<5
    bf=f2-f1;% bandwidth
    if bf<1
        df=min(bf/10,df);    
    end
    N=ceil(fs/df);% Nfft    
end

X=fft(x,N); % N points FFT of the signal
A=1/sqrt(3);
B=1/sqrt(3);% parameter for window
Fs=f2-f1;
s=pi*sqrt(1/3+(B^2-16*A*B)/(4*A^2+2*B^2)/pi^2)/Fs;

w0=(f2+f1)/2*s;% central frequency in Hz
Ts=s/sqrt(3);
W=zeros(1,N);
w=[w0/s-1/2/s+df:df:w0/s+1/2/s];
% fl>=0
if w0/s-1/2/s<0
    w=[df:df:w0/s+1/2/s];
end

for i=round(w/df)+1
    W(i)=A+B*cos(2*pi*s*((i-1)*df-w0/s));
%     W(i)=W(i)/(2*A);% normalization
end
XW=X.*conj(W);
x_hwt0=ifft(XW);
x_hwt=x_hwt0(1:L);

% % plot wavelet
% figure,subplot(211),plot([0:df:fs/2-df],W([1:fs/2/df])),
% xlabel('frequency(Hz)'),title('Harmonic wavelet');
% t=[-4+1/fs:1/fs:4];
% Lt=length(t);
% w=ifft(W,Lt);% time signal
% w2=[conj(w(4*fs:-1:1)),w(1:4*fs)];
% subplot(212),plot(t,real(w2),'r'),hold on,plot(t,imag(w2)),
% xlabel('time(s)'),legend('real','imag');



