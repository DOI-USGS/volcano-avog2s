%#!/usr/bin/octave
%  or to use matlab
% cat plot_Nby2D_tloss.m | matlab -nodesktop -nosplash

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Edit only the fields below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmethid= 1;
srcnam ='Clev.';
srclon =-169.945;
srclat = 52.822;
srcz   = 1.73;
Az     = 41.0;
freq   = 0.1;
na=361;
nr=1000;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nmethid == 1
  NMeth  = 'Strat-Modess';  pref='Nby2D_tloss';      sufx='nm';
elseif nmethid == 2
  NMeth  = 'Strat-CModess'; pref='Nby2D_tloss';      sufx='cnm';
elseif nmethid == 3
  NMeth  = 'Strat-WMod';    pref='Nby2D_wtloss';     sufx='nm';
else
  quit;
end

dim1='_1d.';

file1=sprintf('%s%s%s',pref,dim1,sufx);
file3=sprintf('%s%slossless.%s',pref,dim1,sufx);
dat1=load(file1);
if exist(file3)==2
  HasLossless = 1;
  dat3=load(file3);
else
  HasLossless = 0;
end

r1=dat1(:,1);
s1=r1/(2.0*6371.0*pi)*(360.0);
az1=dat1(:,2);
x1=r1.*sin(az1.*pi/180.);
y1=r1.*cos(az1.*pi/180.);
tl1_Re=dat1(:,3);
tl1_Im=dat1(:,4);
tl1=10*log10(tl1_Re.^2 + tl1_Im.^2);
inco=10*log10(dat1(:,5).^2);

r2=dat3(:,1);
az2=dat3(:,2);
tl2_Re=dat3(:,3);
tl2_Im=dat3(:,4);
tl2=10*log10(tl2_Re.^2 + tl2_Im.^2);

R1=reshape(r1,nr,na);
S1=reshape(s1,nr,na);
AZ1=reshape(az1,nr,na);
X1=reshape(x1,nr,na);
Y1=reshape(y1,nr,na);
TL1=reshape(tl1,nr,na);

figure
xmin=min(x1);
xmax=max(x1);
ymin=min(y1);
ymax=max(y1);
surf(X1,Y1,TL1,'edgecolor','none');
view(2)
axis equal;
h=colorbar;
ylabel(h, 'Tloss (Db)')
colormap(flipud(hot));
caxis([-150 -80]);
titstr=sprintf('%s : 2d transmission loss: %.2f Hz\nSource=%s',NMeth,freq,srcnam);
title(titstr);

%print -dpng DLL.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This section resamples the azimuth and range data onto a lon/lat grid
%  writing out a netcdf
%  This section requires the mapping toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%outdat=zeros(nr*na,3);
%for ia = 1:na
% for ir = 1:nr
%   ii = nr*(ia-1)+ir;
%   [outdat(ii,1) outdat(ii,2)] = reckon(srclat,srclon,S1(ir,ia),AZ1(ir,ia));
%   outdat(ii,3) = TL1(ir,ia);
% end
%end
%
%[xq,yq] = meshgrid(-180:0.02:-150, 45:0.02:65);
%vq=griddata(outdat(:,2),outdat(:,1),outdat(:,3),xq,yq);  
%
%surf(xq,yq,vq,'edgecolor','none')
%view(2)
%colormap(flipud(hot));
%
%a=size(xq);
%ncid=netcdf.create('tloss.nc','CLOBBER');
%londimid=netcdf.defDim(ncid,'lon',a(2));
%latdimid=netcdf.defDim(ncid,'lat',a(1));
%lonvarid=netcdf.defVar(ncid,'lon','NC_DOUBLE',londimid);
%latvarid=netcdf.defVar(ncid,'lat','NC_DOUBLE',latdimid);
%tlsvarid=netcdf.defVar(ncid,'tloss','NC_DOUBLE',[londimid latdimid]);
%netcdf.endDef(ncid);
%netcdf.putVar(ncid,lonvarid,xq(1,:));
%netcdf.putVar(ncid,latvarid,yq(:,1));
%netcdf.putVar(ncid,tlsvarid,vq');
%netcdf.close(ncid);


