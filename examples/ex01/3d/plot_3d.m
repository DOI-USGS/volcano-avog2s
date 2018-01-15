clear all;

raypath_data = load('Clev_raypaths.dat');
ncol = min(size(raypath_data));
if ncol==5
 r  = raypath_data(:,1);          % r [km]
 z  = raypath_data(:,2);          % z, altitude [km]
 dB = raypath_data(:,3);          % Geo. Atten. [dB]
 Aa = raypath_data(:,4);          % Atmo. Atten. [dB]
 Tt = raypath_data(:,5);          % Travel Time [s]
elseif ncol==6
 z   = raypath_data(:,1);         % z, altitude [km]
 lat = raypath_data(:,2);         % Latitude [deg]
 lon = raypath_data(:,3);         % Longitude [deg]
 dB  = raypath_data(:,4);         % Geo. Atten. [dB]
 Aa  = raypath_data(:,5);         % Atmo. Atten. [dB]
 Tt  = raypath_data(:,6);         % Travel Time [s]
end

figure;
scatter3(lat,lon,z,1,'k')

