%#!/usr/bin/octave

clear all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Edit only the fields below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nmethid= 5;
srcnam ='Clev.';
srclon =-169.945;
srclat = 52.822;
srcz   = 1.73;
Az     = 41.0;
freq   = 0.1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nmethid == 1
  NMeth  = 'Strat-Modess';  pref='tloss';      sufx='nm';
elseif nmethid == 2
  NMeth  = 'Strat-CModess'; pref='tloss';      sufx='cnm';
elseif nmethid == 3
  NMeth  = 'Strat-WMod';    pref='wtloss';     sufx='nm';
elseif nmethid == 4
  NMeth  = 'RngDep-Modess'; pref='tloss_rd';   sufx='nm';
elseif nmethid == 5
  NMeth  = 'RngDep-pape';   pref='tloss';      sufx='pe';
end

dim1='_1d.';
dim2='_2d.';

file1=sprintf('%s%s%s',pref,dim1,sufx);
file2=sprintf('%s%s%s',pref,dim2,sufx);
file3=sprintf('%s%slossless.%s',pref,dim1,sufx);
dat1=load(file1);
dat2=load(file2);
if exist(file3)==2
  HasLossless = 1;
  dat3=load(file3);
else
  HasLossless = 0;
end

r1=dat1(:,1);
nr=max(size(dat1));
nc=min(size(dat1));
tl1_Re=dat1(:,2);
tl1_Im=dat1(:,3);
tl1=10*log10(tl1_Re.^2 + tl1_Im.^2);
if nc>3
  IsStrat=1;
  inco=10*log10(dat1(:,4).^2);
else
  IsStrat=0;
end

r2=dat2(:,1);
z2=dat2(:,2);
tl2_Re=dat2(:,3);
tl2_Im=dat2(:,4);
tl2=10*log10(tl2_Re.^2 + tl2_Im.^2);

if HasLossless == 1
  r3=dat3(:,1);
  tl3_Re=dat3(:,2);
  tl3_Im=dat3(:,3);
  tl3=10*log10(tl3_Re.^2 + tl3_Im.^2);
end

nz=max(size(dat2))/nr;
R2 =reshape(r2,nz,nr);
Z2 =reshape(z2,nz,nr);
TL2=reshape(tl2,nz,nr);

fig1=figure;
subplot(10,1,2:6),imagesc(R2(1,:),Z2(:,1),TL2,[-130 -90]);
ylabel('Altitude (km)')
set(gca,'XTick',[])
titstr=sprintf('%s : 2d transmission loss: %.2f Hz\nSource=%s, Az=%0.1f',NMeth,freq,srcnam,Az);
title(titstr);
axis xy;
h=colorbar('SouthOutside');
xlabel(h, 'Tloss (Db)')
colormap(flipud(hot));
subplot(10,1,8:10),plot(r1,tl1,'r-'); hold on;
if HasLossless == 1
  subplot(10,1,8:10),plot(r3,tl3,'b-');
end
if IsStrat==1
  subplot(10,1,8:10),plot(r1,inco,'g-');
end
% Now sort out legend
if HasLossless == 1 && IsStrat==1
  legend('lossy','lossless','incoher.');
elseif HasLossless == 1 && IsStrat==0
  legend('lossy','lossless');
else
  legend('lossy');
end


xlabel('Range (km)')
ylabel('Tloss (Db)');
hold off;

%print "-S750,250" -dpng DLL.png
%print -dpng DLL.png


