%Based on the scripts at the following webpage: https://ecco.jpl.nasa.gov/products/latest/user-guide/
%Add paths
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
allFLD=nan(1080,360,288);
%% Reading 2D fields (ETAN is the liquid sea surface height)
for MONTH=1:288
%fileName= [dirv4r3 'nctiles_monthly\ETAN\' 'ETAN'];
%fldName='ETAN';
%fileName= [dirv4r3 'nctiles_monthly\PHIBOT\' 'PHIBOT'];
%fldName='PHIBOT';
%fld=read_nctiles(fileName,fldName,MONTH); %Read in the MONTH-th monthly record of ETAN
%[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fld); 

%Last two use:
 fileName= [dirv4r3 'nctiles_monthly\oceTAUX\' 'oceTAUX'];
 fldName='oceTAUX';
fldx=read_nctiles(fileName,fldName,MONTH); %Read in the MONTH-th monthly record of ETAN
fileName= [dirv4r3 'nctiles_monthly\oceTAUY\' 'oceTAUY'];
fldName='oceTAUY';
fldy=read_nctiles(fileName,fldName,MONTH); %Read in the MONTH-th monthly record of ETAN
[fldUe,fldVn]=calc_UEVNfromUXVY(fldx,fldy);
[X,Y,FLD]=convert2pcol(mygrid.XC,mygrid.YC,fldUe); %Ue for east Vn for North

allFLD(:,:,MONTH)=FLD;
end
%[allFLDeof,allFLDpc,allFLDexpvar]=eof(allFLD);
%allFLDpc1=allFLDpc(1,:);
latrange=[-40 -60];
index=[find(lat==latrange(end)):find(lat==latrange(1))] ;
for t=1:288
    depthline=find(dens_bnds>1036,1);
    psi1(t)=min(squeeze(PSItot(find(lat==-60),depthline:end,t)));
    psi2(t)=min(min(PSItot(index,depthline:end,t),[],2));
    psi3(t)=max(min(PSItot(index,depthline:end,t),[],2));
end


r=nan(1080,360);
pval=nan(1080,360);
for i=1:1080
    for j=1:360
        [M,P]=corrcoef(squeeze(allFLD(i,j,:)),psi3');
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
title('Correlation between ETAN and PSI1 time series');










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
