
 rayfilename=sprintf('%s_raypaths.dat',Pref);
 fid=fopen(rayfilename,'r');

%rayfilename=sprintf('raypaths.dat');

% First order of business is to determine if we are dealing with a 2d (5-column)
% of a 3d (6-column) case
fid=fopen(rayfilename);
tline = fgetl(fid); % the first line is a comment line
tline = fgetl(fid);
A = sscanf(tline,'%f %f %f %f %f %f');
ncols=max(size(A));
fclose(fid);
if ncols==5
 is2d=1;
else
 is2d=0;
end

% Next loop through file again, count rays and find the longest
fid=fopen(rayfilename);
tline = fgetl(fid);
it = 0;
tmax = 0;
nrays=0;
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
    % If the line has no length, then reset 'it' to start a new ray
    nrays=nrays+1;
    it = 0;
  end
end
fclose(fid);

% Initialize arrays and loop back through the file
r = NaN([tmax nrays]); % Range [km]
x = NaN([tmax nrays]); % Long [deg]
y = NaN([tmax nrays]); % Lat [deg]
z = NaN([tmax nrays]); % z [km]
t = NaN([tmax nrays]); % Travel Time [s]
Ga= NaN([tmax nrays]); % Atmospheric Atten. [dB]
Aa= NaN([tmax nrays]); % Geometric Atten. [dB]

% Finally, read the file once more filling the variables
fid=fopen(rayfilename);
tline = fgetl(fid);
it = 0;
ir = 1;
while 1
  tline = fgetl(fid);
  check=min(size(tline));
  if ~ischar(tline), break, end
  if check == 1
    A = sscanf(tline,'%f %f %f %f %f %f');
    it = it + 1;
    if is2d
      % 2-d case
      r(it,ir) = A(1);
      z(it,ir) = A(2);
      t(it,ir) = A(5);
    else
      % 3-d case
      z(it,ir) = A(1);
      y(it,ir) = A(2);
      x(it,ir) = A(3);
      t(it,ir) = A(6);
       % get range using haversine formula (https://en.wikipedia.org/wiki/Great-circle_distance)
      phi1=Srcy * pi/180.0;
      phi2=A(2) * pi/180.0;
      dellam = (Srcx-A(3)) * pi/180.0;
      r(it,ir) = 2.0*asin(sqrt(sin(0.5*(phi1-phi2)).^2.0 + cos(phi1)*cos(phi2)*sin(dellam*0.5).^2.0))*6371.0;
    end
  else
    % If the line has no length, then reset 'it' to start a new ray
    it = 0;
    ir = ir + 1;
  end
end
fclose(fid);


%nrays=0;
%while 1
%  tline = fgetl(fid);
%  check=min(size(tline));
%  if ~ischar(tline), break, end
%  if check == 1
%    A = sscanf(tline,'%f %f %f %f %f %f');
%    z(it,itheta,iphi) = A(1);
%    y(it,itheta,iphi) = A(2);
%    x(it,itheta,iphi) = A(3);
%    t(it,itheta,iphi) = A(6);
%    it = it + 1;
%  else
%    % If the line has no length, then reset 'it' to start a new ray
%    nrays=nrays+1;
%    itheta=itheta+1;
%    it = 1;
%    if itheta > ntheta
%      itheta = 1;
%      iphi = iphi + 1;
%    end
%  end
%  %disp(tline);
%end
%fclose(fid);
%
%zbounce = 1.0;
%nb = 0;
%for ip = 1:nphi
%  for it = 1:ntheta
%    smax = max(size(z(:,it,ip)));
%    for is = 1:smax
%      if z(is,it,ip)<=zbounce
%        nb = nb+1;
%        xb(nb) = x(is,it,ip);
%        yb(nb) = y(is,it,ip);
%        tb(nb) = t(is,it,ip);
%      end
%    end
%  end
%end
%zb = xb*0.0 + zbounce;
%
%
%Clevx=-169.945;
%Dilx= -158.4575;
%Okmx= -168.175;
%SndPtx=-160.49;
%Aktx=-165.99;
%Adkx=-176.6581;
%
%% Dillingham plot
%clf;
%hold on;
%nplotphi   = nphi;
%nplottheta = ntheta;
%ip=1; % to DLL (Az=42.0 )
%for it = 1:nplottheta
%  plot(x(:,it,ip),z(:,it,ip),'k-')
%end
%axis([-177.0,-158,0,180])
%plot(Clevx,0.0,'r^','MarkerSize',10,'MarkerFaceColor','r');
%plot(Dilx,0.0,'gs','MarkerSize',10,'MarkerFaceColor','g')
%ylabel('Altitude (km)')
%titlestr=sprintf('Clev. to Dillingham; Azimuth = 41.0');
%title(titlestr);
%hold off
%print "-S750,250" -dpng DLL.png
