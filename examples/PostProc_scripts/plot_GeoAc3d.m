%#!/usr/bin/octave
%  or to use matlab
% cat plot_GeoAc3d.m | matlab -nodesktop -nosplash

clear all;

theta_min  =   0.0;
theta_max  =  10.0;
theta_step =   0.5;
phi_min    =  75.0;
phi_max    = 230.0;
phi_step   =   5.0;

xvec=load('Bogo.xloc');
xmin=min(xvec);
xmax=max(xvec);
yvec=load('Bogo.yloc');
ymin=min(yvec);
ymax=max(yvec);
zmin=0.0;
zmax=140.0;

ntheta = (theta_max-theta_min)/theta_step +1;
theta  = linspace(theta_min,theta_max,ntheta);
nphi   = (phi_max-phi_min)/phi_step+1;
phi    = linspace(phi_min,phi_max,nphi);

% Loop through file and find the longest ray
fid=fopen('Bogo_raypaths.dat');
tline = fgetl(fid);
iphi   = 1;
itheta = 1;
it = 0;
tmax = 0;
while 1
  tline = fgetl(fid);
  check=min(size(tline));
  if ~ischar(tline), break, end
  if check == 1
    A = sscanf(tline,'%f %f %f %f %f %f');
    it = it + 1;
    if it > tmax
      tmax = it;
    end
  else
    it = 0;
  end
end
fclose(fid);

%Initialize arrays and loop back through the file
x = NaN([tmax ntheta nphi]);
y = NaN([tmax ntheta nphi]);
z = NaN([tmax ntheta nphi]);
t = NaN([tmax ntheta nphi]);
fid=fopen('Bogo_raypaths.dat');
tline = fgetl(fid);
iphi   = 1;
itheta = 1;
it     = 1;
while 1
  tline = fgetl(fid);
  check=min(size(tline));
  if ~ischar(tline), break, end
  if check == 1
    A = sscanf(tline,'%f %f %f %f %f %f');
    x(it,itheta,iphi) = A(1);
    y(it,itheta,iphi) = A(2);
    z(it,itheta,iphi) = A(3);
    t(it,itheta,iphi) = A(3);
    it = it + 1;
  else
    itheta=itheta+1;
    it = 1;
    if itheta > ntheta
      itheta = 1;
      iphi = iphi + 1;
    end
  end
  %disp(tline);
end
fclose(fid);

 %Bogoslof : -168.0381 53.9280
Bogox = -1198.70; Bogoy = -3680.87;

 %Okmok : -168.1750 53.4681
Okx = -1223.97; Oky = -3728.21;

 %Makusin : -166.9311 53.8864
Makx = -1128.75; Maky = -3707.91;


figure;
hold on;
for ip = 1:nphi
  for it = 1:ntheta
    %scatter3(x(:,it,1),y(:,it,1),z(:,it,1),1,'k');
    plot3(x(:,it,ip),y(:,it,ip),z(:,it,ip),'k-')
  end
end
plot3(Bogox,Bogoy,0,'m^');
plot3(Okx,Oky,0,'m^');
plot3(Makx,Maky,0,'m^');
hold off;
axis([xmin xmax ymin ymax zmin zmax])
grid on
axis equal;



%raypath_data = load('Bogo_raypaths.dat');

%# x [km]       y [km]  z [km]  Geo. Atten. [dB]        Atmo. Atten. [dB]       Travel Time [s]
%-1198.67        -3680.87        1.36764e-05     0       -8.98466e-10    0.0749126
%-1198.65        -3680.87        5.48002e-05     0       -1.7985e-09     0.149957
%x   = raypath_data(:,1); 
%y   = raypath_data(:,2);
%z   = raypath_data(:,3);         % z, altitude [km]
%dB  = raypath_data(:,4);         % Geo. Atten. [dB]
%Aa  = raypath_data(:,5);         % Atmo. Atten. [dB]
%Tt  = raypath_data(:,6);         % Travel Time [s]

%figure;
%scatter3(x,y,z,1,'k')

