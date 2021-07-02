
function QPC=QPC_whwt_bicoh(X,Y,FreqBand1,FreqBand2,FreqWidth,sf,ST,SurrNum)

% Bicoeherence/ Bispectrum 

% In this program we will have two part: part 1 for X-X-X (X series); 
% another part is for X-X-Y (X, Y series)

% input:
% X,Y: signal to be analysed
% FreqBand1,FreqBand2: cental frequencies to be intersted 
% FreqWidth: band width
% sf: sampling frequency
% ST: ST=1: means it conduct the statistc analysis;
%     ST=0: means it does not conduct the statistic analysis 
% SurrNum: number of surrogate data, when ST=1

% output:
% QPC.BS,QPC.BC:bispectrum and bicoherence
% QPC.ST_BS,QPC.ST_BC:statisitic test

if ~exist('SurrNum')
    SurrNum=19;
end

if ~exist('ST')
    ST=0;
end

N=length(FreqBand1);
M=length(FreqBand2);

BS=zeros(N,M);
BC=zeros(N,M);
ST_BS=zeros(N,M);
ST_BC=zeros(N,M);

if mean(X-Y)==0
    Part=1; 
else
    Part=2;
end

% Part 1:  one channel

if Part==1
    for ii=1:N
        f1=FreqBand1(ii); 
        for hh=ii:M
            f2=FreqBand2(hh); 
            % call whwt_bich 
%             [bicoh,bispect,St_bicoh,St_bispect]=whwt_bicoh_A(X,Y,f1,f2,FreqWidth,sf,ST);
            
            if ii==hh
               [bicoh,bispect,St_bicoh,St_bispect,bic0]=whwt_bicoh_C(X,Y,f1,f2,FreqWidth,sf,ST,SurrNum);
               bicoh=bicoh/2; bispect=bispect/2;St_bicoh=St_bicoh/2; St_bispect=St_bispect/2;
               bic0=bic0/2;
            else
               [bicoh,bispect,St_bicoh,St_bispect,bic0]=whwt_bicoh_C(X,Y,f1,f2,FreqWidth,sf,ST,SurrNum);
            end
            
                        
            BS(ii,hh)=bispect;
            BC(ii,hh)=bicoh;
            ST_BS(ii,hh)=St_bispect;
            ST_BC(ii,hh)=St_bicoh;
            BC0(ii,hh)=bic0;
        end
    end

QPC.BS=[BS+BS'];
QPC.BC=[BC+BC'];

QPC.ST_BS=[ST_BS+ST_BS'];
QPC.ST_BC=[ST_BC+ST_BC'];

QPC.BIC0=[BC0+BC0'];

end

% Part 2: two channels

if Part==2
    for ii=1:N
        f1=FreqBand1(ii); 
        for hh=1:M
            f2=FreqBand2(hh); 
            % call whwt_bich 
            [bicoh,bispect,St_bicoh,St_bispect,bic0]=whwt_bicoh_C(X,Y,f1,f2,FreqWidth,sf,ST,SurrNum);
                        
            BS(ii,hh)=bispect;
            BC(ii,hh)=bicoh;
            ST_BS(ii,hh)=St_bispect;
            ST_BC(ii,hh)=St_bicoh;
            BC0(ii,hh)=bic0;
        end
    end

QPC.BS=BS;
QPC.BC=BC;

QPC.ST_BS=ST_BS;
QPC.ST_BC=ST_BC;

QPC.BIC0=BC0;


end




return;  

