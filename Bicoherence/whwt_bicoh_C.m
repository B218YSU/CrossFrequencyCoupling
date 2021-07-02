
function [bicoh,bispect,ST_bicoh,ST_bispect,bic0]=whwt_bicoh_C(x,y,f1,f2,fre_wid,sf,ST,NumSurrogate)

% This function is to calculate the bicoherence between two series x, and y
% at the frequency f1 and f2 at given frequency width;

% If x=y, means auto bicoherence.

% Inoput:

%       x, y    : two time series;
%       fs      : sampling frequency;
%       f1, f2, : central frequency for x and y
%       fre_wid : frequency width

%       St: statistic test

% Output:

%       bispect   : bispectrum;
%       bicoh     : bicoherence;

%       St_bicoh   : statisitic test
%       St_bispect : statisitic test

% ST: ST=1: means it conduct the statistc analysis;
%     ST=0: means it does not conduct the statistic analysis


if ~exist('NumSurrogate')
    NumSurrogate=19;
end

if ~exist('ST')
    ST=0;
end

bispect=0; bicoh=0;
St_bicoh=0; St_bispect=0;

% WHWT

Wx1=whwt(x,f1-fre_wid/2,f1+fre_wid/2,sf);
Wx2=whwt(y,f2-fre_wid/2,f2+fre_wid/2,sf);
Wc=whwt(y,(f1+f2)-fre_wid/2,(f1+f2)+fre_wid/2,sf);

% Bispectrum / bicoherence:
[bispect,bicoh,ST_bispect,ST_bicoh,bic0]=BiSC(Wx1,Wx2,Wc,sf,ST,NumSurrogate);

return;


function [bis,bic,ST_bis,ST_bic,bic0]= BiSC(Wx1,Wx2,Wc,sf,ST,NumSurrogate)

% Mean method is aplied to reduce the noise:
% The series is divided into SN segments,

% SN=10, in fact we may revise it based on the observisions on practical data
% Liang revised 10
% SN  = 20;
SN=10;
Len = length(Wx1);
WN = floor(Len/SN);

% Segment number:

for jj=1:SN
    NN=((jj-1)*WN+1):jj*WN;
    [bs(jj),bc(jj),ST_bs(jj),ST_bc(jj),bc0(jj)]= BiSCP(Wx1(NN),Wx2(NN),Wc(NN),sf,ST,NumSurrogate);
end

% trimmend mean:
% K=10;
% bis=trimmean(bs,K);
% bic=trimmean(bc,K);

% mean:
bis=mean(bs);
bic=mean(bc);

ST_bis=mean(ST_bs);
ST_bic=mean(ST_bc);

bic0=mean(bc0);

% bic=trimmean(bc,10);

return;


function [bis,bic,ST_bis,ST_bic,bic0]= BiSCP(Wx1,Wx2,Wc,sf,ST,NumSurrogate)

%
bis=0;
bic=0;
ST_bis=0;
ST_bic=0;
bic0=0;

% Cross - spectrum

Wxxy=(Wx1.*Wx2.*conj(Wc));
% Bi-spectrum: version 1
bis=abs(sum(Wxxy));

s=sum(abs(Wc).^2);
ss=sum(abs(Wx1.*Wx2).^2);

% Bi-coherence
biphase=angle(Wxxy);

R=pi*randn(length(biphase),1)';
Rx=exp(sqrt(-1)*(R.*biphase));

% Calcualtion of bicoherence
bxxy=sum(abs(Wxxy).*Rx);
bic=sqrt((abs(bxxy).^2)/[s*ss+0.005]);

bxxy=sum(abs(Wxxy));
bic0=sqrt((abs(bxxy).^2)/[s*ss+0.005]);


if ST>0
    kk=NumSurrogate;
    S=0; W=0;

    for i=1:kk

        R=pi*randn(length(Wx1),1)';
        bip=biphase+R;
        Rx=exp(sqrt(-1)*(R.*bip));

        Sbxxy=sum(abs(Wxxy).*Rx);
        Sbic=sqrt((abs(Sbxxy).^2)/[s*ss+0.005]);

        Wbis(i)=abs(Sbxxy);
        Wbic(i)=Sbic;
    end


    g=2;
    if bic>(mean(Wbic)+g*std(Wbic))
        ST_bic=bic;
    else
        ST_bic=0; end

    abis=abs(Wbis);

    if bis>(mean(abis)+g*std(abis))
        ST_bis=bis;
    else
        ST_bis=0;
    end

end


% % Bi-spectrum: version 2: Aug. 12 2008,
% % Ref. : Anesthesiolgy, v 100 No. 4, 2004
%
%
%     TP=(Wx1.*Wx2.*conj(Wc));
%     bis=abs(sum(TP));           %    bispecrum
%
%     % Bi-coherence
%     biphase=angle(TP);
%     R=pi*randn(length(TP),1)';
%     Rx=exp(sqrt(-1)*(R.*biphase));
%
%     % Calcualtion of bicoherence
%     R_bis=abs(sum(abs(TP).*Rx));
%
%     Nm=sum(abs(TP));
%     bic=R_bis/Nm;
%
%
%     if ST>0
%        kk=NumSurrogate;
%        S=0; W=0;
%
%        for i=1:kk
%
%            R=pi*randn(length(TP),1)';
%            bip=biphase+R;
%            Rx=exp(sqrt(-1)*(R.*bip));
% %            Rx=exp(sqrt(-1)*(bip));
%
%            Sbxxy=abs(sum(abs(TP).*Rx));
%            Sbic=Sbxxy/Nm;
%
%            Wbis(i)=Sbxxy;
%            Wbic(i)=Sbic;
%        end
%
%
%      g=2;
%      if bic>(mean(Wbic)+g*std(Wbic))
%          ST_bic=bic;
%      else
%          ST_bic=0; end
%
%      abis=abs(Wbis);
%
%      if abs(bis)>(mean(abis)+g*std(abis))
%          ST_bis=bis;
%      else
%          ST_bis=0; end
%
%
%

return;



