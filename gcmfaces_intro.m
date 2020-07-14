%Based on the scripts at the following webpage: https://ecco.jpl.nasa.gov/products/latest/user-guide/
%Add paths
addpath 'C:\Users\user\Downloads' %where we have calcalphabeta, dens and sw_cp
p = genpath('C:\cygwin64\home\user\ECCO\gcmfaces\'); addpath(p); %Change to whatever folder you put gcmfaces in (Andrew sent a scp/download link to the folder on Caolila)
p = genpath('m_map/'); addpath(p); %I do not think you need m_map in this particular case, but at some point it may become useful to download. It makes useful plots using geographic projections
dirv4r3 = 'C:\cygwin64\home\user\ECCO\Version4\Release3\'; %Change to you equivalent data folder
dirGrid = 'C:\cygwin64\home\user\ECCO\Version4\Release3\nctiles_grid';

%MONTH = 3; %March, for example

%% "global" loads
%'straight','cube','compact'
grid_load; %Load the grid.
gcmfaces_global; % Define global variables
%nFaces,fileFormat,memoryLimit,omitNativeGrid
%You will be prompted for the path of the nctiles_grid directory. Type it in
alldensflux=nan(1080,360,288);
allFLD2=nan(1080,360,288);
%% Reading 2D NC fields (ETAN is the liquid sea surface height)
for MONTH=1:288
%fileName= [dirv4r3 'nctiles_monthly\ETAN\' 'ETAN'];
%fldName='ETAN';
%fileName= [dirv4r3 'nctiles_monthly\PHIBOT\' 'PHIBOT'];
%fldName='PHIBOT';
fileName= [dirv4r3 'nctiles_monthly\surfSALT\' 'surfSALT'];
fldName='SALT';

fld1=read_nctiles(fileName,fldName,MONTH,1); %need to plug in 1 for surface depth if its a 4d field
[X,Y,FLD1]=convert2pcol(mygrid.XC,mygrid.YC,fld1); 

fileName= [dirv4r3 'nctiles_monthly\surfTemp\' 'surfTemp'];
fldName='THETA';
fld2=read_nctiles(fileName,fldName,MONTH,1); %Read in the MONTH-th monthly record of ETAN
[X,Y,FLD2]=convert2pcol(mygrid.XC,mygrid.YC,fld2); 

fileName= [dirv4r3 'nctiles_monthly\SFLUX\' 'SFLUX'];
fldName='SFLUX';
fld3=read_nctiles(fileName,fldName,MONTH); %need to plug in 1 for surface depth if its a 4d field
[X,Y,FLD3]=convert2pcol(mygrid.XC,mygrid.YC,fld3); 
fileName= [dirv4r3 'nctiles_monthly\TFLUX\' 'TFLUX'];
fldName='TFLUX';
fld4=read_nctiles(fileName,fldName,MONTH); %need to plug in 1 for surface depth if its a 4d field
[X,Y,FLD4]=convert2pcol(mygrid.XC,mygrid.YC,fld4); 
%Vectors use:
% fileName= [dirv4r3 'nctiles_monthly\oceTAUX\' 'oceTAUX'];
 %fldName='oceTAUX';
%fldx=read_nctiles(fileName,fldName,MONTH); %Read in the MONTH-th monthly record of ETAN
%fileName= [dirv4r3 'nctiles_monthly\oceTAUY\' 'oceTAUY'];
%fldName='oceTAUY';
%fldy=read_nctiles(fileName,fldName,MONTH); %Read in the MONTH-th monthly record of ETAN
%[fldUe,fldVn]=calc_UEVNfromUXVY(fldx,fldy);
%[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fldUe); %Ue for east Vn for North

dens=densmdjwf(FLD1,FLD2,10*ones(1080,360));
cp=sw_cp(FLD1,FLD2,10*ones(1080,360));

[alpha,beta]=calcAlphaBeta(FLD1,FLD2,5*ones(1080,360));
densflux=beta.*FLD3-alpha.*FLD4./cp; %change this to sflux and tflux
alldensflux(:,:,MONTH)=densflux;
end




%for .mat file use:

%UVELtot= load('C:\cygwin64\home\user\ECCO\Version4\Release3\myproducts_monthly\UVELtot_mean.mat');
%UVELtot=UVELtot.UVELtot_mean;
%[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,UVELtot); 

%[allFLDeof,allFLDpc,allFLDexpvar]=eof(allFLD);
%allFLDpc1=allFLDpc(1,:);



%AMOC
latrange=[-40 -60];
index=[find(lat==latrange(end)):find(lat==latrange(1))] ;
for t=1:288
    depthline=find(dens_bnds>1036,1);
    psi1(t)=min(squeeze(PSI_notrend(find(lat==-50),depthline:end,t)));
    psi2(t)=min(min(PSI_notrend(index,depthline:end,t),[],2));
    psi3(t)=max(min(PSI_notrend(index,depthline:end,t),[],2));
end


r=nan(1080,360);
pval=nan(1080,360);
for i=1:1080
    for j=1:360
        [M,P]=corrcoef(squeeze(alldensflux(i,j,:)),psi3');
        r(i,j)=M(2);
        pval(i,j)=P(2);
    end
end
%[M,P]=corrcoef(SAM,pc1'); %corrcoef for SAM with pc1
%[M,P]=corrcoef(SAM,allFLDpc1'); %correcoef for SAM with field pc1




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
title('Correlation between Densflux and PSI3 time series');










%{
if 0

    %Read in binary files
    %To read in and display the 15th 6-hourly record of the E-W wind stress for year 1992
    %Read in data
    dirv4r3 = '/mydir/v4r3/';
    ff= [dirv4r3 'input_forcing/' 'eccov4r3_ustr_1992'];
    fld = read_bin(ff,15,0); %read in the 15th 2-D record

    %Display
    figure;
    [X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fld); pcolor(X,Y,FLD);
    if ~isempty(find(X>359)); axis([0 360 -90 90]); else; axis([-180 180 -90 90]); end;
    dd1 = 1;
    cc=[-1:0.1:1]*dd1; % color bar set to -1 to 1 N/m2
    shading flat; cb=gcmfaces_cmap_cbar(cc);
    %Add labels and title
    xlabel('longitude'); ylabel('latitude');
    title('display using convert2pcol');


end
%}
