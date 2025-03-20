%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.09                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***           okMakeDataSubset.m            ***             %
%             ***********************************************             %
%                                                                         %
% Make subset of the input data structure from okLoadData.m               %
% Cropping the data to make it smaller                                    %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement        %
% 3. rmin: Numeric. Minimum row index                                     %
% 4. rmax: Numeric. Maximum row index                                     %
% 5. cmin: Numeric. Minimum col index                                     %
% 6. cmax: Numeric. Maximum col index                                     %
%                                                                         %
% Example:                                                                %
% DataStruct = okMakeDataSubset(DataStruct,'Original',1,3000,1,5000)      %
%                                                                         %
% Output:                                                                 %
% 1. DataSubset: Structure. Cropped data will be stored in Field 'Subset' %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DataSubset = okMakeDataSubset(DataStruct,Dataset,rmin,rmax,cmin,cmax)

[Row,Col] = size(DataStruct.(Dataset).Displ);
disp('*** Making subset of the input data structure')
disp(strcat('*** Orig size:',32,num2str(Row),32,num2str(Col)))
disp('**************************************')
DataSubset = DataStruct;
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Displ'))
    Displ = DataStruct.(Dataset).Displ;
    DisplSub = Displ(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Displacement')
    DataSubset.Subset.Displ = DisplSub;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Azimuth'))
    Azimuth = DataStruct.(Dataset).Azimuth;
    AzimuthSub = Azimuth(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Azimuth')
    DataSubset.Subset.Azimuth = AzimuthSub;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Incidence'))
    Incidence = DataStruct.(Dataset).Incidence;
    IncidenceSub = Incidence(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Incidence')
    DataSubset.Subset.Incidence = IncidenceSub;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Longitude'))
    Longitude = DataStruct.(Dataset).Longitude;
    LongitudeSub = Longitude(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Longitude')
    DataSubset.Subset.Longitude = LongitudeSub;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Latitude'))
    Latitude = DataStruct.(Dataset).Latitude;
    LatitudeSub = Latitude(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Latitude')
    DataSubset.Subset.Latitude = LatitudeSub;
end
if any(strcmp(fieldnames(DataStruct.(Dataset)),'Geometry'))
    Geometry = DataStruct.(Dataset).Geometry;
    GeometrySub = Geometry(rmin:rmax,cmin:cmax);
    disp('*** Making subset of Geometry')
    DataSubset.Subset.Geometry = GeometrySub;
end

% Store the subset index boundaries
DataSubset.Subset.Boundaries = [rmin,rmax,cmin,cmax];

% Report a bit of attributes
disp('**************************************')
disp(strcat('Subset size:',32,num2str(rmax-rmin+1),32,num2str(cmax-cmin+1)))
disp(' ')
end


