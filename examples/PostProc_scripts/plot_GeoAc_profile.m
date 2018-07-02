%#!/usr/bin/octave
%  or to use matlab
% cat plot_GeoAc_profile.m | matlab -nodesktop -nosplash


%Pref='Clev0';
Pref='Bering';


% Source location
%  Note: the source lon/lat is needed for calculating ranges and azimuths to stations
%        the source z is needed for calculating sound_speed at source in load_GeoAc_atmo.m
Srcx=-169.945;  Srcy=52.822; Srcz=1.73;

% Get the location of the infrasound arrays
load_StationInfo

% Here's where you need to know a little bit about how the GeoAc job was run
theta_min  = -30.0; theta_max  =  55.0; theta_step =   1.0;
%phi_min    =  40.0; phi_max    =  90.0; phi_step   =   5.0;
% For a single azimuth, use this
phi_min    =  41.0; phi_max    =  41.0; phi_step   =   0.0;


ntheta = (theta_max-theta_min)/theta_step +1;
theta  = linspace(theta_min,theta_max,ntheta);
if phi_step>0.0  % for azimuthal sweeps
  nphi = (phi_max-phi_min)/phi_step+1;
else             % for a single azimuth
  nphi = 1;
end
phi    = linspace(phi_min,phi_max,nphi);

% Now load the raypath file and double-check that the number of rays equal the phi/theta info
% from above
load_GeoAc_raypaths
if nrays ~= nphi*ntheta
  nrays
  nphi
  ntheta
  error('ERROR: nrays does not equal nphi*ntheta')
end

Z=reshape(z,tmax,ntheta,nphi);
Y=reshape(y,tmax,ntheta,nphi);
X=reshape(x,tmax,ntheta,nphi);
T=reshape(t,tmax,ntheta,nphi);
R=reshape(r,tmax,ntheta,nphi);

% Now load atmosphere file
load_GeoAc_atmo


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dillingham plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% We need to know which phi() is closest to the azimuth to DLL.
% If it is outside some tolerance, then reject
%
tol=5.0;
[dum ip] = min(abs(phi-Dil_az));
if abs(dum) > tol
  error('ERROR: swath of profiles > 5-degrees from requested azimuth')
end

vmin=-100.0; vmax=100.0;
cmin=200.0; cmax=400.0;
Tmin=150.0; Tmax=400.0;
zmin=0.0; zmax = 140.0;
rmin=0.0; rmax = 1000.0;
clf;
subplot(1,10,1),plot(a_u,a_z,'r',a_v,a_z,'b',wind_dir_prop,a_z,'k'); hold on;
axis([vmin,vmax,zmin,zmax]);
set(gca,'XTick',[]);
ylabel('Altitude (km)');
xlabel('Vel');
hold off;

subplot(1,10,2),plot(a_cth,a_z,'k',a_cef,a_z,'k--',[sound_speed sound_speed],[zmin zmax],'g-'); hold on;
axis([cmin,cmax,zmin,zmax]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);
xlabel('C');
hold off;

subplot(1,10,3),plot(a_T,a_z,'k'); hold on;
axis([Tmin,Tmax,zmin,zmax]);
set(gca,'YTick',[]);
set(gca,'XTick',[])
xlabel('T');
hold off;

for it = 1:ntheta
  subplot(1,10,4:10),plot(r(:,it,ip),z(:,it,ip),'k-'); hold on;
end
axis([rmin,rmax,zmin,zmax]);
subplot(1,10,4:10),plot(0.0,0.0,'r^','MarkerSize',5,'MarkerFaceColor','r');
subplot(1,10,4:10),plot(Dil_r,0.0,'gs','MarkerSize',5,'MarkerFaceColor','g');
set(gca,'YTick',[])
xlabel('Range');
hold off
%print "-S750,350" -dpng Clev_DLL_GeoAc.png
%print('Clev_DLL_GeoAc.png','-dpng')


