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
% 4. Processor: Character. Specify InSAR processor                        %
%    ['ISCE','GMTSAR' etc.]                                               %
% 5. LookDirection: Character. InSAR looking direction                    %
%    ['r','l']                                                            %
%                                                                         %
% Case 2: GNSS data (Not supported)                                       %
% 1. Displacement: Matrix. The GNSS displacements and its error estimates %
%                                                                         %
% Example (InSAR-ISCE):                                                   %
% DataStruct = okLoadData('filt_topophase_unw.geo','Azimuth.tif', ...     %
%       'Incidence.tif', 'ISCE', 'r')                                     %
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
%    DataStruct.Original.Displ = Displacement                             %
%    DataStruct.Original.Longitude = Longitude                            %
%    DataStruct.Original.Latitude = Latitude                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataStruct = okLoadData(varargin)
% Parse the parameters
p = inputParser;
default_displ = 0;
default_azi = 0;
default_inc = 0;
default_coh = 0;
default_processor = 'ISCE';
default_lookdir = 'r';

addParameter(p, 'Displacement', default_displ, @(x) ischar(x) || ismatrix(x));
addParameter(p, 'Azimuth', default_azi, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Incidence', default_inc, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Coherence', default_coh, @(x)  ischar(x) || ismatrix(x));
addParameter(p, 'Processor', default_processor, @(x) ischar(x));
addParameter(p, 'LookDir', default_lookdir, @(x) ischar(x));
parse(p, varargin{:});

Displarg = p.Results.Displacement;
Aziarg = p.Results.Azimuth;
Incarg = p.Results.Incidence;
Coharg = p.Results.Coherence;
Processor = p.Results.Processor;
LookDirection = p.Results.LookDir;


if nargin >= 5
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
    % Report a bit of attributes
    disp('*******************************************')
    disp(strcat('Size of the input matrix:',32,num2str(size(Displ))))
    disp(strcat('Azimuth angle:',32,num2str(mean(Azimuth(:),'omitnan'))))
    disp(' ')

else
    Displarg = varargin{1};
    flag = ones(1,2);
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
    DataStruct.Original.Longitude = LonMat;
    DataStruct.Original.Latitude = LatMat;
    % Report a bit of attributes
    disp('*******************************************')
    disp(strcat('Size of the input matrix:',32,num2str(size(Displ))))
    disp(' ')
    
    
    error_message = ['Input arguments should be either 5, 6 or 1: \n' ...
        '  5 input arguments: \n' ...
        '    "Displacement", "UnwrapPhase.tif" \n' ...
        '    "Azimuth", "Azimuth.tif" \n'
        '    "Coherence", "Coherence.tif" (This is optional) \n'...
        '    "Incidence", "Incidence.tif" \n' ...
        '    "Processor", "ISCE" \n' ...
        '    "LookDir", "r" \n' ...
        '  1 input arguments (This is for GNSS observations. However, currently not supported):\n' ...
        '    Displacement \n'];
    error('u:stuffed:it', error_message);
end

end
