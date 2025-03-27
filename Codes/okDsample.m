%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.06                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okDsample.m                ***             %
%             ***********************************************             %
%                                                                         %
% (Update 2025.03.08)                                                     %
%   Change function output. Also support downsample the viewing geometry  %
% (Update 2025.03.26)                                                     %
%   Populate InSAR auto-correlation result to the output and some minor   %
%   changes                                                               %
%                                                                         %
% Perform InSAR data downsampling based on the input downsampling function%
% to prepare for fault slip inversion                                     %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement, and   %
%    the local coordinates                                                %
% 3. FuncParam: Structure. Containing the downsampling parameters. See:   %
%    okMakeDsampleParam.m                                                 %
%                                                                         %
% Example:                                                                %
% DataStruct = okDsample(DataStruct,'Original',FuncParam)                 %
% The output data strucuture will be updated with the downsampled data    % 
%                                                                         %
% Output:                                                                 %
% 1. DsampleDataStruct: Structure. Downsampled displacement and its       %
%    attributes. Will be stored in field 'Dsample'                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DsampleDataStruct = okDsample(DataStruct,Dataset,FuncParam)
% These three must exist, otherwise report error and escape downsampling
Data = DataStruct.(Dataset).Displ;
LocalXMat = DataStruct.(Dataset).LocalX;
LocalYMat = DataStruct.(Dataset).LocalY;

if any(strcmp(fieldnames(DataStruct.(Dataset)),'LongitudeOrigin'))
    LonOrigin = DataStruct.(Dataset).LongitudeOrigin;
else
    LonOrigin = [];
end

% To fit for okMakeGreenFunc.m
% Case 1:
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Azimuth'))
    Azimuth = DataStruct.(Dataset).Azimuth;
    Aziflag = 1;
else
    Aziflag = 0;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Incidence'))
    Incidence = DataStruct.(Dataset).Incidence;
    Incflag = 1;
else
    Incflag = 0;
end
% Case 2:
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Ecoeff'))
    Ecoeff = DataStruct.(Dataset).Ecoeff;
    Eflag = 1;
else
    Eflag = 0;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Ncoeff'))
    Ncoeff = DataStruct.(Dataset).Ncoeff;
    Nflag = 1;
else
    Nflag = 0;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Ucoeff'))
    Ucoeff = DataStruct.(Dataset).Ucoeff;
    Uflag = 1;
else
    Uflag = 0;
end

% Also downsample longitude, latitude and coherence(optional)
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Longitude'))
    Longitude = DataStruct.(Dataset).Longitude;
    Lonflag = 1;
else
    Lonflag = 0;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Latitude'))
    Latitude = DataStruct.(Dataset).Latitude;
    Latflag = 1;
else
    Latflag = 0;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Coherence'))
    Coherence = DataStruct.(Dataset).Coherence;
    Cohflag = 1;
else
    Cohflag = 0;
end

% Parse the function type and retrieve corresponding parameters
Dfunc = FuncParam.FunctionType;

if strcmp(Dfunc,'quadtree') || strcmp(Dfunc,'Quadtree')
    % Function parameters
    Criterion = FuncParam.Criterion;
    Tolerance = FuncParam.Tolerance;
    LeafMinSize = FuncParam.LeafMinSize;
    LeafMaxSize = FuncParam.LeafMaxSize;
    NaNPixelAllow = FuncParam.NaNPixelAllow;

    % Report function inputs
    disp(' ')
    disp('******* Downsampling data okDsample.m *******')
    disp(strcat('*** Downsample function:',32,Dfunc))
    disp(strcat('*** Downsample criterion:',32,Criterion))
    disp(strcat('*** Threshold tolerance:',32,num2str(Tolerance)))
    disp(strcat('*** Minimum leaf allowed:',32,num2str(LeafMinSize)))
    disp(strcat('*** Maximum leaf allowed:',32,num2str(LeafMaxSize)))
    disp(strcat('*** NaN pixel allowed per leaf:',32,num2str(NaNPixelAllow)))
    disp(' ')


    % Matrix boundaries, downsample the whole matrix
    rstart = 1;
    rend = size(Data,1);
    cstart = 1;
    cend = size(Data,2);
    LeafCoord = okQuadtreeD(Data,rstart,rend,cstart,cend,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow);

    % Downsample the data based on the quadtree result
    LeafCount = size(LeafCoord,1);
    disp('*********************************')
    disp(strcat('*** Total leaf count:',32,num2str(LeafCount)));
    disp(' ')

    DsampleGeom = zeros(size(LeafCoord,1),8);
    DsampleData = zeros(size(LeafCoord,1),3);
    for i = 1:size(LeafCoord,1)
        rmin = LeafCoord(i,1);
        rmax = LeafCoord(i,2);
        cmin = LeafCoord(i,3);
        cmax = LeafCoord(i,4);
        rmid = floor((rmin+rmax)/2);
        cmid = floor((cmin+cmax)/2);

        % Displacements
        xCoord = mean(LocalXMat(rmin:rmax,cmin:cmax),'all','omitnan');
        yCoord = mean(LocalYMat(rmin:rmax,cmin:cmax),'all','omitnan');
        % Calculate the NaN pixel counts
        NaNPixelRatio = sum(sum(isnan(Data(rmin:rmax,cmin:cmax)))) / (size(Data(rmin:rmax,cmin:cmax),1)*size(Data(rmin:rmax,cmin:cmax),2));

        % If larger than the allowed, then return NaN
        if NaNPixelRatio >= NaNPixelAllow
            DataVal = nan;
        else
            DataVal = mean(Data(rmin:rmax,cmin:cmax),'all','omitnan');
        end
        DsampleData(i,:) = [yCoord,xCoord,DataVal];

        % Geometry and Lon Lat Coherence if any
        if Aziflag == 1
            tmpAzi = mean(Azimuth(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,1) = tmpAzi;
        end
        if Incflag == 1
            tmpInc = mean(Incidence(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,2) = tmpInc;
        end
        if Eflag == 1
            tmpEcoeff = mean(Ecoeff(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,3) = tmpEcoeff;
        end
        if Nflag == 1
            tmpNcoeff = mean(Ncoeff(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,4) = tmpNcoeff;
        end
        if Uflag == 1
            tmpUcoeff = mean(Ucoeff(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,5) = tmpUcoeff;
        end
        if Lonflag == 1
            tmpLoncoeff = mean(Longitude(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,6) = tmpLoncoeff;
        end
        if Latflag == 1
            tmpLatcoeff = mean(Latitude(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,7) = tmpLatcoeff;
        end
        if Cohflag == 1
            tmpCohcoeff = mean(Coherence(rmin:rmax,cmin:cmax),'all','omitnan');
            DsampleGeom(i,8) = tmpCohcoeff;
        end

    end

    % Find the NaN values and remove them before output
    % This is to ensure the generation of Green's function and inversion
    % are correct
    rmNaN = find(~isnan(DsampleData(:,3)));

    % Populate the input to the output
    DsampleDataStruct = DataStruct;
    DsampleDataStruct.Dsample = DataStruct.(Dataset);

    % Output downsampled data and downsample parameters strucuture 
    DsampleDataStruct.Dsample.Displ = DsampleData(rmNaN,3);
    DsampleDataStruct.Dsample.LocalX = DsampleData(rmNaN,2);
    DsampleDataStruct.Dsample.LocalY = DsampleData(rmNaN,1);
    DsampleDataStruct.Dsample.LeafCoord = LeafCoord;

    % Output downsampled geometry and Lon Lat Coherence
    if Aziflag == 1
        DsampleDataStruct.Dsample.Azimuth = DsampleGeom(rmNaN,1);
    end
    if Incflag == 1
        DsampleDataStruct.Dsample.Incidence = DsampleGeom(rmNaN,2);
    end
    if Eflag == 1
        DsampleDataStruct.Dsample.Ecoeff = DsampleGeom(rmNaN,3);
    end
    if Nflag == 1
        DsampleDataStruct.Dsample.Ncoeff = DsampleGeom(rmNaN,4);
    end
    if Uflag == 1
        DsampleDataStruct.Dsample.Ucoeff = DsampleGeom(rmNaN,5);
    end
    if Lonflag == 1
        DsampleDataStruct.Dsample.Longitude = DsampleGeom(rmNaN,6);
    end
    if Latflag == 1
        DsampleDataStruct.Dsample.Latitude = DsampleGeom(rmNaN,7);
    end
    if Cohflag == 1
        DsampleDataStruct.Dsample.Coherence = DsampleGeom(rmNaN,8);
    end

    % Output the downsample parameters
    DsampleDataStruct.Dsample.Criterion = Criterion;
    DsampleDataStruct.Dsample.Tolerance = Tolerance;
    DsampleDataStruct.Dsample.LeafMinSize = LeafMinSize;
    DsampleDataStruct.Dsample.LeafMaxSize = LeafMaxSize;
    DsampleDataStruct.Dsample.NaNPixelAllow = NaNPixelAllow;
    DsampleDataStruct.Dsample.DsampleFlag = 1;
    DsampleDataStruct.Dsample.SourceDataset = Dataset;
    if ~isempty(LonOrigin)
        DsampleDataStruct.Dsample.LongitudeOrigin = LonOrigin;
    end



elseif strcmp(Dfunc,'uniform') || strcmp(Dfunc,'Uniform')
    % To-be updated
    % Make a okUniformD.m

elseif strcmp(Dfunc,'rosampling') || strcmp(Dfunc,'Rosampling') || strcmp(Dfunc,'R-based')
    % To-be updated

else
    error_message = ['Cannot find corrsponding downsampler. For different downsampling algorithm, see below: \n', ...
        '  Quadtree: put quadtree or Quadtree \n', ...
        '  Uniform sampling: put uniform or Uniform \n', ...
        '  R-based sampling: put rosampling or Rosampling or R-based (Lohman and Simons, 2005) \n', ...
        'Or see okMakeDsampleParam.m to make a compatible one for it'];

    error('u:stuffed:it', error_message);

end



end
