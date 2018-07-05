%#!/usr/bin/octave
%  or to use matlab
% cat plot_art2d.m | matlab -nodesktop -nosplash

clear all;

vref=0.33;  % reducing velocity
xoff=2651;
xoff2=2500;
%close all


file1   = 'art2d_xsec.bin';
filein  = 'art2d_data.in';
filell  = 'ProfLonLat01.dat';
filetopo= 'topo.in';

  %da = 180.0/(nprofile);
  %azim = 0.0 + (iprof-1)*da;
azim = 41.0;
%ttitle = sprintf('%s: Rays above Ceff, (Topo/Reflections, Az=%.1f)',DetectorName,azim);

% do binary read on file1/2
fid=fopen(file1,'rb');
buffer=fread(fid,inf,'double');
ANG=buffer(1:6:end);
X=buffer(2:6:end);
Y=buffer(3:6:end);
Z=buffer(4:6:end);
T=buffer(5:6:end);
R=buffer(6:6:end);
fclose(fid);
hi=[ANG X-xoff2 Z T R Y];

A=load('dims.txt');
nx = A(1);
nz = A(2);

A=load(filein);
x10  = A(:,1);
y10  = A(:,2);
z10  = A(:,3);
vx10 = A(:,4);
vy10 = A(:,5);

c=reshape(z10,[nz nx]);
vx=reshape(vx10,[nz nx]);
vy=reshape(vy10,[nz nx]);
x=reshape(x10,[nz nx]);
z=reshape(y10,[nz nx]);

figure;

vmin=-100.0; vmax=100.0;
cmin=200.0; cmax=400.0;
Tmin=150.0; Tmax=400.0;
zmin=0.0; zmax = 140.0;
rmin=0.0; rmax = 1000.0;
clf;
subplot(1,10,1),plot(1000.0*vx(:,2385),z(:,2385),'k'); hold on;
axis([vmin,vmax,zmin,zmax]);
set(gca,'XTick',[]);
ylabel('Altitude (km)');
xlabel('Vel');
hold off;

subplot(1,10,2),plot(1000.0*c(:,2385),z(:,2385),'k'); hold on;
axis([cmin,cmax,zmin,zmax]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlabel('C');
hold off;

subplot(1,10,3),plot(c(:,2385),z(:,2385),'k'); hold on;
axis([Tmin,Tmax,zmin,zmax]);
set(gca,'YTick',[]);
set(gca,'XTick',[])
xlabel('T');
hold off;


%hold on;
for i=-180:1:180
  j=find(hi(:,1)==i); subplot(1,10,4:10),plot(hi(j,2),hi(j,3),'k-','linewidth',0.05); hold on;
end
axis([0 1000.0 0 140]);
set(gca,'YTick',[])
%ylabel('Altitude (km)'); xlabel('Range (km)'); 
%hold off
%print "-S750,350" -dpng DLL.png

