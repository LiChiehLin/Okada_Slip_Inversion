%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.14                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***             okLocalXYtoLL.m             ***             %
%             ***********************************************             %
%                                                                         %
% Convert local XY to Lon Lat cartesian coordinate system                 %
% Translate from php to Matlab function                                   %
% Translate from:(https://gist.github.com/pingyen/1346895)                %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. lon_origin: Numeric. The conversion longitude origin                 %
%                                                                         %
% Example:                                                                %
% Coord_LL = okLocalXYtoLL(120,'localx',DataStruct.Subset.LocalX, ...     %
%       'localy',DataStruct.Subset.LocalY)                                %
% DataStruct = okLocalXYtoLL(120,'localx',DataStruct.Subset.LocalX, ...   %
%       'localy',DataStruct.Subset.LocalY,DataStruct.Subset)              %
%                                                                         %
% Output:                                                                 %
% 1. DataLL: Structure or matrix depending on whether or not input data   %
%    structure. Converted coordinates will be stored as                   %
%    Lon: (Longitude)                                                     %
%    Lat: (Latitude)                                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataLL = okLocalXYtoLL(lon_origin,varargin)
p = inputParser;
default_localx = 0;
default_localy = 0;
default_datastruct = 0;
addParameter(p, 'localx', default_localx, @(x) isnumeric(x));
addParameter(p, 'localy', default_localy, @(x) isnumeric(x));
addParameter(p, 'datastruct', default_datastruct, @(x) isstruct(x));

parse(p, varargin{:});
x = p.Results.localx;
y = p.Results.localy;
DataStruct = p.Results.datastruct;


disp('*** Converting from Local XY back to Lon Lat')
disp(strcat('*** Local longitude origin:',32,num2str(lon_origin)))
disp(' ')


a = 6378137.0; %Equatorial radius
b = 6356752.314245; %Polar radius
long0 = lon_origin*pi/180; %Central longitude (radian)
k0 = 0.9999; %Scale factor
dx = 250000; 

dy = 0;
e = (1 - b.^2./a.^2).^0.5;

lon = x-dx;
lat = y-dy;

M = lat ./ k0;

mu = M ./ (a .* (1.0 - e.^2 ./ 4.0 - 3 .* e.^4 ./ 64.0 - 5 .* e.^6 ./ 256.0));
e1 = (1.0 - (1.0 - e.^2).^0.5) / (1.0 + (1.0 - e.^2).^0.5);

J1 = (3 .* e1 ./ 2 - 27 .* e1.^3 ./ 32.0);
J2 = (21 .* e1.^2 / 16 - 55 .* e1.^4 ./ 32.0);
J3 = (151 .*e1.^3 ./ 96.0);
J4 = (1097 .* e1.^4 ./ 512.0);

fp = mu + J1 .* sin(2 .* mu) + J2 .* sin(4 .* mu) + J3 .* sin(6 .* mu) + J4 .* sin(8 .* mu);

e2 = (e .* a ./ b).^2;
C1 = e2 .* cos(fp).^2;
T1 = tan(fp).^2;
R1 = a .* (1 - e.^2) ./ (1 - e.^2 * sin(fp).^2).^1.5; %%%%%%
N1 = a ./ (1 - e.^2 .* sin(fp).^2).^0.5; %%%%%%%%

D = lon ./ (N1 .* k0);

Q1 = N1 .* tan(fp) ./ R1;
Q2 = D.^2 ./ 2.0;
Q3 = (5 + 3 .* T1 + 10 .* C1 - 4 .* C1.^2 - 9 .* e2) .* D.^4 ./ 24.0;
Q4 = (61 + 90 .* T1 + 298 .* C1 + 45 .* T1.^2 - 3 .* C1.^2 - 252 .* e2) .* D.^6 ./ 720.0;
lat = fp - Q1 .* (Q2 - Q3 + Q4);

Q5 = D;
Q6 = (1 + 2 .* T1 + C1) .* D.^3 ./ 6;
Q7 = (5 - 2 * C1 + 28 * T1 - 3 * C1.^2 + 8 * e2 + 24 * T1.^2) .* D.^5 / 120.0;
lon = long0 + (Q5 - Q6 + Q7) ./ cos(fp);

lat = (lat .* 180) ./ pi;
lon = (lon .* 180) ./ pi;

% Ouptut
if isstruct(DataStruct)
    DataLL = DataStruct;
    DataLL.Lon = lon;
    DataLL.Lat = lat;
else
    DataLL = [lon,lat];
end


end