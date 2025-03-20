%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.02.26                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***             okLLtoLocalXY.m             ***             %
%             ***********************************************             %
%                                                                         %
% Convert Lon Lat to local XY cartesian coordinate system                 %
% Translate from Python to Matlab function                                %
% Translate from:(http://blog.ez2learn.com/2009/08/15/lat-lon-to-twd97/)  %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement        %
% 3. lon_origin: Numeric. The conversion longitude origin                 %
%                                                                         %
% Example:                                                                %
% DataStruct = okLLtoLocalXY(DataStruct,'Original',120)                   %
% DataStruct = okLLtoLocalXY(DataStruct,'Subset',120)                     %
%                                                                         %
% Output:                                                                 %
% 1. DataLocal: Structure. Converted coordinates will be stored as        %
%    LocalX: (Longitude)                                                  %
%    LocalY: (Latitude)                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataLocal = okLLtoLocalXY(DataStruct,Dataset,lon_origin)

disp('*** Converting from Lon Lat to local XY')
disp(strcat('*** Local longitude origin:',32,num2str(lon_origin)))
disp(' ')
% Retrieve Longitude and Latitude
lon = DataStruct.(Dataset).Longitude;
lat = DataStruct.(Dataset).Latitude;


a = 6378137.0; %Equatorial radius
b = 6356752.314245; %Polar radius
long0 = lon_origin*pi/180; %Central longitude (radian)
k0 = 0.9999; %Scale factor
dx = 250000; 

lon = lon.*(pi/180); %Convert to radian
lat = lat.*(pi/180); %Convert to radian

e = (1-b.^2./a.^2).^0.5;
e2 = e.^2./(1-e.^2);
n = (a-b)./(a+b);
nu = a./(1-(e.^2).*(sin(lat).^2)).^0.5;
p = lon-long0;

A = a.*(1 - n + (5/4.0).*(n.^2 - n.^3) + (81/64.0).*(n.^4  - n.^5));
B = (3.*a.*n./2.0).*(1 - n + (7/8.0).*(n.^2 - n.^3) + (55/64.0).*(n.^4 - n.^5));
C = (15.*a.*(n.^2)./16.0).*(1 - n + (3./4.0).*(n.^2 - n.^3));
D = (35.*a.*(n.^3)./48.0).*(1 - n + (11/16.0).*(n.^2 - n.^3));
E = (315.*a.*(n.^4)/51.0).*(1 - n);

S = A.*lat - B.*sin(2.*lat) + C.*sin(4.*lat) - D.*sin(6.*lat) + E.*sin(8.*lat);

K1 = S.*k0;
K2 = k0.*nu.*sin(2.*lat)/4.0;
K3 = (k0.*nu.*sin(lat).*(cos(lat).^3)/24.0) .* (5 - tan(lat).^2 + 9.*e2.*(cos(lat).^2) + 4.*(e2.^2).*(cos(lat).^4));

y = K1 + K2.*(p.^2) + K3.*(p.^4);

K4 = k0.*nu.*cos(lat);
K5 = (k0.*nu.*(cos(lat).^3)./6.0) .* (1 - tan(lat).^2 + e2.*(cos(lat).^2));

x = K4.*p + K5.*(p.^3) + dx;

lon = x;
lat = y;

% Ouptut
DataLocal = DataStruct;
DataLocal.(Dataset).LocalX = lon;
DataLocal.(Dataset).LocalY = lat;
DataLocal.(Dataset).LongitudeOrigin = lon_origin;
end