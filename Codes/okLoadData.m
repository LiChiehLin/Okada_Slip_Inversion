%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.09                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okLoadData.m               ***             %
%             ***********************************************             %
%                                                                         %
% (Update: 2025.05.13)                                                    %
%   Supports GNSS input                                                   %
%                                                                         %
% Read in data and make the data strucuture that complies with later      %
% processing                                                              %
% Currently only supports output from ISCE                                %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% Case 1: InSAR data                                                      %
% 1. Displacement: Character or matrix. The InSAR LOS or Azi displacement %
% 2. Azimuth: Character or matrix. The azimuth direction                  %
% 3. Inicidence: Character or matrix. The Incidence direction             %
% 4. Coherence: Character or matrix. The coherence                        %
% 5. Processor: Character. Specify InSAR processor                        %
%    ['ISCE','GMTSAR' etc.]                                               %
% 6. LookDirection: Character. InSAR looking direction                    %
%    ['r','l']                                                            %
%                                                                         %
% Example (InSAR-ISCE):                                                   %
% DataStruct = okLoadData('Displacement','UnwrapPhase.tif', ...           %
%    'Azimuth','Azimuth.tif', ...                                         %
%    'Incidence','Incidence.tif', ...                                     %
%    'Coherence','Coherence.tif', ...                                     %
%    'Processor','ISCE', ...                                              %
%    'LookDir','r');                                                      %
%                                                                         %
% Case 2: GNSS data (Not supported)                                       %
% 1. Longitude: Vector. The longitude of the GNSS stations                %
% 2. Latitude: Vector. The latitude of the GNSS stations                  %
% 3. DisplacementEW: Vector. The EW displacement                          %
% 4. DisplacementNS: Vector. The NS displacement                          %
% 5. DisplacementUD: Vector. The UD displacement                          %
% 6. UncertaintyEW: Vector. The EW displacement uncertainty               %
% 7. UncertaintyNS: Vector. The NS displacement uncertainty               %
% 8. UncertaintyUD: Vector. The UD displacement uncertainty               %
% 9. StationName: Character. The station names                            %
%                                                                         %
% Any input left blank will be output as NaN                              %
%                                                                         %
% Note that all read-in data will be stored in Field 'Original'           %
%                                                                         %
% Output:                                                                 %
% Case 1 (InSAR-ISCE):                                                    %
%    DataStruct.Original.Displ = Displacement                             %
%    DataStruct.Original.Azimuth = Azimuth                                %
%    DataStruct.Original.Incidence = Incidence                            %
%    DataStruct.Original.Coherence = Coherence                            %
%    DataStruct.Original.Longitude = Longitude                            %
%    DataStruct.Original.Latitude = Latitude                              %
% Case 2 (GNSS):                                                          %
%    DataStruct.Original.Longitude = Longitude                            %
%    DataStruct.Original.Latitude = Latitude                              %
%    DataStruct.Original.DisplEW = DisplacementEW                         %
%    DataStruct.Original.DisplNS = DisplacementNS                         %
%    DataStruct.Original.DisplUD = DisplacementUD                         %
%    DataStruct.Original.UncertEW = UncertaintyEW                         %
%    DataStruct.Original.UncertNS = UncertaintyNS                         %
%    DataStruct.Original.UncertUD = UncertaintyUD                         %
%    DataStruct.Original.StationName = StationName                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataStruct = okLoadData(varargin)
% Parse the parameters
p = inputParser;
% InSAR
default_displ = 0;
default_azi = 0;
default_inc = 0;
default_coh = 0;
default_processor = 'ISCE';
default_lookdir = 'r';
% GNSS
default_lon = nan;
default_lat = nan;
default_displEW = nan;
default_displNS = nan;
default_displUD = nan;
default_uncertEW = nan;
default_uncertNS = nan;
default_uncertUD = nan;
default_staName = nan;

% InSAR input parameter parser
addParameter(p, 'Displacement', default_displ, @(x) ischar(x) || ismatrix(x));
addParameter(p, 'Azimuth', default_azi, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Incidence', default_inc, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Coherence', default_coh, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Processor', default_processor, @(x) ischar(x));
addParameter(p, 'LookDir', default_lookdir, @(x) ischar(x));

% GNSS input parameter parser
addParameter(p, 'Longitude', default_lon, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'Latitude', default_lat, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'DisplacementEW', default_displEW, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'DisplacementNS', default_displNS, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'DisplacementUD', default_displUD, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'UncertaintyEW', default_uncertEW, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'UncertaintyNS', default_uncertNS, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'UncertaintyUD', default_uncertUD, @(x) ismatrix(x) && isnumeric(x));
addParameter(p, 'StationName', default_staName, @(x) iscell(x));
parse(p, varargin{:});

% InSAR
Displarg = p.Results.Displacement;
Aziarg = p.Results.Azimuth;
Incarg = p.Results.Incidence;
Coharg = p.Results.Coherence;
Processor = p.Results.Processor;
LookDirection = p.Results.LookDir;

% GNSS
Lonarg = p.Results.Longitude;
Latarg = p.Results.Latitude;
DisplEWarg = p.Results.DisplacementEW;
DisplNSarg = p.Results.DisplacementNS;
DisplUDarg = p.Results.DisplacementUD;
UncertEWarg = p.Results.UncertaintyEW;
UncertNSarg = p.Results.UncertaintyNS;
UncertUDarg = p.Results.UncertaintyUD;
StaNamearg = p.Results.StationName;

if ischar(Displarg) || ismatrix(Displarg) && Displarg ~= 0
    % InSAR inputs
    flag = ones(4,1);
    disp(' ')
    disp(strcat('InSAR processor:',Processor,32,'Looking direction:',LookDirection))
    % Displacement
    if ischar(Displarg)
        disp(strcat('*** Reading displacement from',32,Displarg))
        [Displ,A] = readgeoraster(Displarg);
        DataStruct.Original.Attributes = A;
    else
        disp('*** Already a read-in matrix (Displacement)')
        Displ = Displarg;
        flag(1) = 0;
    end
    % Azimuth
    if ischar(Aziarg)
        disp(strcat('*** Reading Azimuth from',32,Aziarg))
        [Azimuth,A] = readgeoraster(Aziarg);
        Azimuth(Azimuth==0) = nan;
        % Convert Azimuth angle
        if strcmp(Processor,'ISCE') && strcmp(LookDirection,'r')
            Azimuth = -Azimuth + 90;
        elseif strcmp(Processor,'ISCE') && strcmp(LookDirection,'l')
            Azimuth = -Azimuth - 90;
        end
        DataStruct.Original.Attributes = A;
    else
        disp('*** Already a read-in matrix (Azimuth)')
        Azimuth = Aziarg;
        flag(2) = 0;
    end
    % Incidence
    if ischar(Incarg)
        disp(strcat('*** Reading Incidence from',32,Incarg))
        [Incidence,A] = readgeoraster(Incarg);
        Incidence(Incidence==0) = nan;
        DataStruct.Original.Attributes = A;
    else
        disp('*** Already a read-in matrix (Incidence)')
        Incidence = Incarg;
        flag(3) = 0;
    end
    % Coherence (This is optional)
    if ischar(Coharg)
        disp(strcat('*** Reading Coherence from',32,Coharg))
        [Coherence,A] = readgeoraster(Coharg);
        DataStruct.Original.Attributes = A;
        Cohflag = 1;
    elseif ismatrix(Coharg) && Coharg ~= 0
        disp('*** Already a read-in matrix (Coherence)')
        Coherence = Coharg;
        flag(4) = 0;
        Cohflag = 1;
    elseif Coharg == 0
        disp('*** No coherence data will be read in')
        Cohflag = 0;
        Coherence = nan;
    end
    
    
    % Get Lon Lat if read from Tiff
    if any(flag==1)
        LonMin = A.LongitudeLimits(1);
        LonMax = A.LongitudeLimits(2);
        LatMin = A.LatitudeLimits(1);
        LatMax = A.LatitudeLimits(2);
        LatInc = A.CellExtentInLatitude;
        LonInc = A.CellExtentInLongitude;
        Row = A.RasterSize(1);
        Col = A.RasterSize(2);
    
        Lon = LonMin:LonInc:LonMax-LonInc;
        LonMat = repmat(Lon,Row,1);
        Lat = LatMin:LatInc:LatMax-LatInc;
        LatMat = repmat(flipud(transpose(Lat)),1,Col);
    else
        warning('No Lon Lat information was read!! Output the row and col indices as Longitude and Latitude.')
        [Row,Col] = size(Displ);
        Lon = 1:Col;
        LonMat = repmat(Lon,Row,1);
        Lat = 1:Row;
        LatMat = repmat(flipud(transpose(Lat)),1,Col);
    end
    
    % Output data structure
    DataStruct.Original.Displ = Displ;
    DataStruct.Original.Azimuth = Azimuth;
    DataStruct.Original.Incidence = Incidence;
    if Cohflag == 1
        DataStruct.Original.Coherence = Coherence;
    end
    DataStruct.Original.Longitude = LonMat;
    DataStruct.Original.Latitude = LatMat;
    DataStruct.Processor = Processor;
    DataStruct.LookDirection = LookDirection;
    DataStruct.DataType = 'InSAR';
    % Report a bit of attributes
    disp('*******************************************')
    disp(strcat('Size of the input matrix:',32,num2str(size(Displ))))
    disp(strcat('Azimuth angle:',32,num2str(mean(Azimuth(:),'omitnan'))))
    disp(' ')

elseif Displarg == 0 
    % GNSS inputs
    disp('*** Loading GNSS data')
    DataStruct.Original.Longitude = Lonarg;
    DataStruct.Original.Latitude = Latarg;
    DataStruct.Original.DisplEW = DisplEWarg;
    DataStruct.Original.DisplNS = DisplNSarg;
    DataStruct.Original.DisplUD = DisplUDarg;
    DataStruct.Original.UncertEW = UncertEWarg;
    DataStruct.Original.UncertNS = UncertNSarg;
    DataStruct.Original.UncertUD = UncertUDarg;
    DataStruct.Original.StaNamearg = StaNamearg;
    DataStruct.DataType = 'GNSS';

    % Make some report attributes
    if isnan(Lonarg)
        LonCheck = 'Not loaded';
    else
        LonCheck = 'Check';
    end
    if isnan(Latarg)
        LatCheck = 'Not loaded';
    else
        LatCheck = 'Check';
    end
    if isnan(DisplEWarg)
        DisplEWCheck = 'Not loaded';
    else
        DisplEWCheck = 'Check';
    end
    if isnan(DisplNSarg)
        DisplNSCheck = 'Not loaded';
    else
        DisplNSCheck = 'Check';
    end
    if isnan(DisplUDarg)
        DisplUDCheck = 'Not loaded';
    else
        DisplUDCheck = 'Check';
    end
    if isnan(UncertEWarg)
        UncertEWCheck = 'Not loaded';
    else
        UncertEWCheck = 'Check';
    end
    if isnan(UncertNSarg)
        UncertNSCheck = 'Not loaded';
    else
        UncertNSCheck = 'Check';
    end
    if isnan(UncertUDarg)
        UncertUDCheck = 'Not loaded';
    else
        UncertUDCheck = 'Check';
    end
    if ~iscell(StaNamearg) && isnan(StaNamearg)
        StaNameCheck = 'Not loaded';
    else
        StaNameCheck = 'Check';
    end

    % Report a bit of attributes
    disp('*******************************************')
    disp(strcat('* GNSS data count:',32,num2str(length(Lonarg))))
    disp(strcat('* Longitude:',32,LonCheck))
    disp(strcat('* Latitude:',32,LatCheck))
    disp(strcat('* EW Displacement:',32,DisplEWCheck))
    disp(strcat('* NS Displacement:',32,DisplNSCheck))
    disp(strcat('* UD Displacement:',32,DisplUDCheck))
    disp(strcat('* EW Uncertainty:',32,UncertEWCheck))
    disp(strcat('* NS Uncertainty:',32,UncertNSCheck))
    disp(strcat('* UD Uncertainty:',32,UncertUDCheck))
    disp(strcat('* Station names:',32,StaNameCheck))
else

    error_message = ['Input arguments should be either 5, 6 or 1: \n' ...
        '  5 or 6 input arguments (InSAR arguments): \n' ...
        '    "Displacement", "UnwrapPhase.tif" \n' ...
        '    "Azimuth", "Azimuth.tif" \n'
        '    "Coherence", "Coherence.tif" (This is optional) \n'...
        '    "Incidence", "Incidence.tif" \n' ...
        '    "Processor", "ISCE" \n' ...
        '    "LookDir", "r" \n' ...
        '  9 input arguments (GNSS arguments): \n' ...
        '    "Longitude", Lon \n' ...
        '    "Latitude", Lat \n' ...
        '    "DisplacementEW", EWdispl \n' ...
        '    "DisplacementNS", NSdispl \n' ...
        '    "DisplacementUD", UDdispl \n' ...
        '    "UncertaintyEW",  EWuncert \n' ...
        '    "UncertaintyNS",  NSuncert \n' ...
        '    "UncertaintyUD",  UDuncert \n' ...
        '    "StationName",  StaName \n'];
    error('u:stuffed:it', error_message);
end

end
