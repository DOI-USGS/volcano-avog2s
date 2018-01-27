clear all;

theta_min  =   0.0;
%theta_max  =  10.0;
theta_max  =   0.5;
theta_step =   0.5;
phi_min    =  75.0;
phi_max    = 230.0;
%phi_step   =   5.0;
phi_step   =   1.0;

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

% Get Topography
topox=ncread('3d_tephra.nc','x');
topoy=ncread('3d_tephra.nc','y');
topodata=ncread('3d_tephra.nc','Topography');
%topodata(topodata<0)=0;

 %Bogoslof : -168.0381 53.9280
Bogox = -1198.70; Bogoy = -3680.87;

 %Okmok : -168.1750 53.4681
Okx = -1223.97; Oky = -3728.21;

 %Makusin : -166.9311 53.8864
Makx = -1128.75; Maky = -3707.91;

figure;
hold on;
nplotphi   = nphi;
nplottheta = 1;
for ip = 1:nplotphi
  for it = 1:nplottheta
    %scatter3(x(:,it,1),y(:,it,1),z(:,it,1),1,'k');
    plot3(x(:,it,ip),y(:,it,ip),z(:,it,ip),'k-')
  end
end

surface(topox,topoy,topodata')
shading interp
demcmap(topodata)
plot3(Bogox,Bogoy,3,'m^');
plot3(Okx,Oky,3,'m^');
plot3(Makx,Maky,3,'m^');
hold off;
axis([xmin xmax ymin ymax zmin zmax])
grid on
axis equal;

