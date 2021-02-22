close all;
load 'PSIs.mat';
load 'pcolorxy.mat';
addpath 'C:\Users\user\Downloads' %where we have calcalphabeta, dens and sw_cp
p = genpath('C:\cygwin64\home\user\ECCO\gcmfaces\'); addpath(p); %Change to whatever folder you put gcmfaces in (Andrew sent a scp/download link to the folder on Caolila)
p = genpath('m_map/'); addpath(p); %I do not think you need m_map in this particular case, but at some point it may become useful to download. It makes useful plots using geographic projections
dirv4r3 = 'C:\cygwin64\home\user\ECCO\Version4\Release3\'; %Change to you equivalent data folder
dirGrid = 'C:\cygwin64\home\user\ECCO\Version4\Release3\nctiles_grid';

%MONTH = 3; %March, for example


%'straight','cube','compact'
grid_load; %Load the grid.
gcmfaces_global; % Define global variables
r=nan(1080,360);
pval=nan(1080,360);
Maxlag=nan(1080,360);
for i=1:1080
    for j=1:360
        %Correlation
        %[M,P]=corrcoef(squeeze(alldensflux(i,j,:)),psi3');
        %r(i,j)=M(2);
        %pval(i,j)=P(2);
        
        
        %lagged correlation
        
        [r1,lags]=xcorr(squeeze(allFLD(i,j,:)),psi2','normalized');
        [Maxval,MaxIndex]=max(r1);
        MaxMonth=lags(MaxIndex);
        
       % optimallag=circshift(psi3',MaxMonth);
        %{
        if MaxMonth>0
            optimallag(1:MaxMonth)=0;
        end
        if MaxMonth<0
            optimallag(end+MaxMonth,end)=0;
        end
        %}
       %[M,P]=corrcoef(squeeze(alldensflux(i,j,:)),optimallag);
         
        r(i,j)=Maxval;
        Maxlag(i,j)=MaxMonth;

    end
end

%% Plot - with gcmfaces
%This is much better:
figure;
pcolor(X,Y,r);
if ~isempty(find(X>359)); axis([0 360 -90 90]); else; axis([-180 180 -90 90]); end;
dd1 = 1;
cc=[-1:0.1:1]*dd1; % color bar set to -1 to 1 N/m2
shading flat; cb=gcmfaces_cmap_cbar(cc);
%Add labels and title
xlabel('longitude'); ylabel('latitude');
title('Max Correlation between ETAN and PSI2 time series');

figure;
pcolor(X,Y,Maxlag);
if ~isempty(find(X>359)); axis([0 360 -90 90]); else; axis([-180 180 -90 90]); end;
dd1 = 10;
cc=[-1:0.1:1]*dd1; % color bar set to -1 to 1 N/m2
shading flat; cb=gcmfaces_cmap_cbar(cc);
%Add labels and title
xlabel('longitude'); ylabel('latitude');
title('Max Lag Time between ETAN and PSI2 time series');
%the first argument is t time before second if its -t