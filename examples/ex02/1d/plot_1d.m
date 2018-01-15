clear all;

raypath_data = load('Bogo0_raypaths.dat');
r  = raypath_data(:,1);          % r [km]
z  = raypath_data(:,2);          % z, altitude [km]
dB = raypath_data(:,3);          % Geo. Atten. [dB]
Aa = raypath_data(:,4);          % Atmo. Atten. [dB]
Tt = raypath_data(:,5);          % Travel Time [s]

atmo_data = load('atmo.dat');
z2  = atmo_data(:,1);            % z, altitude [km]
u2  = atmo_data(:,3);            % u(z) [m/s]
v2  = atmo_data(:,4);            % v(z) [m/s]
c2  = atmo_data(:,7)*1000.0;     % c, sound speed based on temperature [m/s]
ceff2  = atmo_data(:,8)*1000.0;  % ceff, effective sound speed [m/s]

% Wind speed in the direction of propagation
wind_dir_prop = ceff2 - c2;

subplot(1,6,1),plot(u2,z2,'r',v2,z2,'b',wind_dir_prop,z2,'k')
subplot(1,6,2:6),scatter(r,z,1,'k');
axis([0 1000.0 0 140]);

