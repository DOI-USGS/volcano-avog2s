
try
 atmofilename=sprintf('%s_atmo.dat',Pref);
 atmo=load(atmofilename);
catch
 atmofilename=sprintf('atmo.dat');
 atmo=load(atmofilename);
end

a_z   = atmo(:,1)-6370.0 ; % z [km]
a_T   = atmo(:,2);         % Temperature [K}
a_u   = atmo(:,3);         % Zonal wind [m/s]
a_v   = atmo(:,4);         % Meridional wind [m/s]
a_rho = atmo(:,5);         % Density [g/cm3]
a_p   = atmo(:,6);         % presure [mbar]
a_cth = atmo(:,7)*1000.0;  % c therm. [m/s]
a_cef = atmo(:,8)*1000.0;  % c eff. [m/s]
% Wind speed in the direction of propagation
wind_dir_prop = a_cef - a_cth;
[dum indx] = min(abs(a_z-Srcz));
sound_speed = a_cth(indx);


