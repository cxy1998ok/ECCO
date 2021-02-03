clearvars -except pc
%close all;
%get SAM Index
SAM=load('SAM_Since1992.mat');
SAM_1=SAM.n';
ttotal=288;
SAM=SAM_1(1:ttotal);
%get PC1
pc1=pc(1,:);
pc2=pc(2,:);
%get Psis
load('PSIs.mat');
%Transfer Function

%Divide into subintervals x_k & y_k
blocklength=12; %length of each subinterval in months

%% no overlap
%
xk=reshape(SAM,blocklength,numel(SAM)/blocklength)'; %each row is the time series of blocklength months
Nblocks=size(xk,1);
yk=reshape(psi2,blocklength,numel(SAM)/blocklength)';
%}
%% overlapping subintervals
%{
overlaprate=0.5;
startmonth=[1:blocklength*overlaprate:ttotal];
Nblocks=length(startmonth)-1;
xk=[];
yk=[];
for i=1:Nblocks
    if startmonth(i)+blocklength-1<=ttotal
    xk=[xk; SAM(startmonth(i):startmonth(i)+blocklength-1)];
    yk=[yk; psi2(startmonth(i):startmonth(i)+blocklength-1)];   
    end
end
%redefine Nblocks for rounding down
Nblocks=size(xk,1);
%}
%% windowing
%detrend
xk_detrend=detrend(xk');
yk_detrend=detrend(yk');
xk_detrend=xk_detrend';
yk_detrend=yk_detrend';
%Hann Window
hann=(sin(pi*[0:blocklength-1]/(blocklength-1))).^2; 
hamming=0.54-0.46*cos(2*pi*[0:blocklength-1]/(blocklength-1));
xk=xk_detrend.*hamming;
yk=yk_detrend.*hamming;
%% FFT
Nyears=24;
Fs=1/3600/24/30; %month in a sec/sampling frequency
N=12*Nyears;
f=[0:blocklength/2]*Fs/blocklength; %frequency for fft
Nf=size(f,2); %number of different freq
xk_f=fft(xk,[],2); %do fft across rows 
yk_f=fft(yk,[],2);

%% Transfer Functions
xk_fs=zeros(Nblocks,blocklength/2+1);
yk_fs=zeros(Nblocks,blocklength/2+1);
Sxx=zeros(1,blocklength/2+1);
Sxy=zeros(1,blocklength/2+1);
Syy=zeros(1,blocklength/2+1);
for i=1:Nblocks %number of subintervals
%Make single-sided
P2=xk_f(i,:)/blocklength;
P1=P2(:,1:blocklength/2+1);
P1(2:end-1)=2*P1(2:end-1);
xk_fs(i,:)=P1;
P2=(yk_f(i,:)/blocklength);
P1=P2(1:blocklength/2+1);
P1(2:end-1)=2*P1(2:end-1);
yk_fs(i,:)=P1;
for j=1:blocklength/2+1
Sxy(j)=Sxy(j)+xk_fs(i,j)'*yk_fs(i,j);
Sxx(j)=Sxx(j)+xk_fs(i,j)'*xk_fs(i,j);
Syy(j)=Syy(j)+yk_fs(i,j)'*yk_fs(i,j);
Txy(j)=Sxy(j)/Sxx(j);
Coherence(j)=abs(Sxy(j))/(sqrt(Sxx(j)*Syy(j)));
Terrors(j)=1/2/Nblocks*(1-(Coherence(j))^2)/(Coherence(j)^2)*abs(Txy(j))^2;
Aerrors(j)=1/2/Nblocks*(1-(Coherence(j))^2)/(Coherence(j)^2);
end
end
Sxx=Sxx/Nblocks;
Sxy=Sxy/Nblocks;


%scaling due to windowing?

%% plotting
f_in_year=f*60*60*24*360; %12*30=360
period=1./f_in_year*12;
%Magnitude
MTxy=abs(Txy);
ATxy=angle(Txy);
ATxy=rad2deg(ATxy);
Terrors(Terrors<=0.01)=0;
Terrors=sqrt(Terrors);
Aerrors(Aerrors<=0.01)=0;
Aerrors=atan(sqrt(Aerrors));
Aerrors=rad2deg(Aerrors);
figure;
yyaxis left
errorbar(period,MTxy,2*Terrors);
xlabel('Period (Months)');
ylabel('TF Magnitude');
yyaxis right
plot(period,ATxy);
ylabel('TF Phase (deg)');
title(['      Transfer Function of Psi2 and SAM in ' num2str(blocklength) '-Month Blocks']);
%xlabel('Frequency')
%ylabel('Phase Angle of TF');
