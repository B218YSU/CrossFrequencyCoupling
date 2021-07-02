
function [P,MOIN] = CalPAC(DataL,DataH)

%  Input:  
%         DataL : The lower frequency components of Data, filtered by FIR filter.
%         DataH : The high frequency components of Data, filtered by FIR filter.

%  Output:
%         P : The normalized mean amplitude over phase-bins
%         MOIN :  Modulation Index (MI)


% Hilbert transform
hl = hilbert(DataL);
hh = hilbert(DataH);

% instantaneous phase & amplitude envelope
pha_l = angle(hl);
amp_h = abs(hh);

N = length(pha_l);

Abin = zeros(18,1);
num = zeros(18,1);
for j = 1:N

    for k = 1:18

        if pha_l(j)>=(-pi+(k-1)*pi/9) && pha_l(j)<(-pi+k*pi/9)

            Abin(k,1) = Abin(k,1) + amp_h(j);
            num(k,1) = num(k,1) + 1;

        end  

    end

end

% normalized mean amplitude over phase-bins
mAbin = Abin./num;
P = mAbin/sum(mAbin);

% MI calculated by KL divergence 
MOIN = 1+sum(P.*log(P))/log(18); 

end