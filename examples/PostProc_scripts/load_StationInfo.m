
% Stations locations
Adkx  =-176.658; Adky   = 51.88;
Clevx =-169.721922; Clevy  = 52.786357;
Okmx  =-168.175; Okmy   = 53.4681;
Aktx  =-165.99;  Akty   = 54.1356;
SndPtx=-160.49;  SndPty = 55.3397;
Dilx  =-158.4575;Dily   = 59.0397;

% Now get the range and azimuth from the src to the stations
phi1=Srcy * pi/180.0; phi2=Adky * pi/180.0; dellam = (Srcx-Adkx) * pi/180.0;
 Adk_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 Adk_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;
phi1=Srcy * pi/180.0; phi2=Clevy * pi/180.0; dellam = (Srcx-Clevx) * pi/180.0;
 Clev_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 Clev_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;
phi1=Srcy * pi/180.0; phi2=Okmy * pi/180.0; dellam = (Srcx-Okmx) * pi/180.0;
 Okm_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 Okm_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;
phi1=Srcy * pi/180.0; phi2=Akty * pi/180.0; dellam = (Srcx-Aktx) * pi/180.0;
 Akt_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 Akt_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;
phi1=Srcy * pi/180.0; phi2=SndPty * pi/180.0; dellam = (Srcx-SndPtx) * pi/180.0;
 SndPt_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 SndPt_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;
phi1=Srcy * pi/180.0; phi2=Dily * pi/180.0; dellam = (Srcx-Dilx) * pi/180.0;
 Dil_r = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
 Dil_az=-atan2( sin(dellam)*cos(phi2),cos(phi1)*sin(phi2)-sin(phi1)*cos(phi2)*cos(dellam))*180.0/pi;

