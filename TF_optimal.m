%clear all;
close all;
%get SAM Index
SAM=load('SAM_Since1992.mat');
SAM_1=SAM.n';
ttotal=288;
SAM=SAM_1(1:ttotal);
%get PC1
pc=load('pc.mat');
pc=pc.pc;
pc1=pc(1,:);
pc2=pc(2,:);
%get Psis
load('PSIs.mat');

%Transfer Function

PeriodList=[12:4:12*24]; % period in months
Nperiod=length(PeriodList);
Txy=zeros(1,Nperiod);
Coherence=zeros(1,Nperiod);
Terrors=zeros(1,Nperiod);
Aerrors=zeros(1,Nperiod);
alpha=1;
for b=1:Nperiod
    blocklength=alpha*PeriodList(b);
    if blocklength>ttotal  %if blocklength is too big, continue
    continue;
    end
    %% overlapping subintervals
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

%% windowing
%detrend
xk_detrend=detrend(xk');
yk_detrend=detrend(yk');
xk_detrend=xk_detrend';
yk_detrend=yk_detrend';
%Hann Window
hann=(sin(pi*[0:blocklength-1]/blocklength)).^2; 
hamming=0.54-0.46*cos(2*pi*[0:blocklength-1]/blocklength);
xk=xk_detrend.*hamming;
yk=yk_detrend.*hamming;
%% FFT
Nyears=24;
Fs=1/3600/24/30; %month in a sec/sampling frequency
N=12*Nyears;
f=[0:blocklength/2]*Fs/blocklength; %frequency for fft


Nf=size(f,2); %number of different freq
xk_f=fft(xk,[],2); %do fft across rows with f points
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
end
end
Sxx=Sxx/Nblocks;
Sxy=Sxy/Nblocks;
Syy=Syy/Nblocks;
periodindex=2; %only want the frequency corresponds to the PeriodList(b) should be 2 for now
Txy(b)=Sxy(periodindex)/Sxx(periodindex); 
Coherence(b)=abs(Sxy(periodindex))/(sqrt(Sxx(periodindex)*Syy(periodindex)));
Terrors(b)=1/2/Nblocks*(1-(Coherence(b))^2)/(Coherence(b)^2)*abs(Txy(b))^2;
Aerrors(b)=1/2/Nblocks*(1-(Coherence(b))^2)/(Coherence(b)^2);
end

%% plotting
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
errorbar(PeriodList,MTxy,2*Terrors);
xlabel('Period (Months)');
ylabel('TF Magnitude');
yyaxis right
%errorbar(PeriodList,ATxy,2*Aerrors);
ylabel('TF Phase (deg)');
title(['      TF of Psi2 and SAM in blocklength=' num2str(alpha) '*period']);
%xlabel('Frequency')
%ylabel('Phase Angle of TF');

