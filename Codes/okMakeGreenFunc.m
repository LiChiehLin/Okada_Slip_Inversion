%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.02.27                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***            okMakeGreenFunc.m            ***             %
%             ***********************************************             %
%                                                                         %
% (Update 2025.10.08)                                                     %
% Support inputs of different Poisson's ratio                             %
% (Update 2025.03.08)                                                     %
% Support converting to LOS or Azimuth displacement vector                %
%                                                                         %
% Make Green's function from the fault model generated from               %
% okMakeFaultModel.m leveraging okada85.m.                                %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement, and   %
%    the local coordinates                                                %
% 3. Direction: Character. Green's function conversion direction          %
%    'LOS': Convert to LOS direction                                      %
%    'Azi': Convert to Azimuth direction                                  %
%    '3D': Remain in EW, NS and UD direction                              %
% 4. FaultModel: Structure. The fault geometry from okMakeFaultModel.m    %
% 5. Rake: Numeric. Rake angle                                            %
% 6. Slip: Numeric. How much the patch slips                              %
% 7. Opening: Numeric. How much the patch tensile opens                   %
% 8. GreenFuncName: Character. Output field name of the Green's function  %
%                                                                         % 
% Example (Strike-slip):                                                  %
% FaultModel = okMakeGreenFunc(DataStruct,'Dsample', ...                  %
%       'LOS', FaultModel, 180, 1, 0, 'GreenLOS'                          %           
%                                                                         %
% Output:                                                                 %
% 1. GreenFunc: Structure. The Green's function will be stored in the     %
%    field name put in GreenFuncName                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function GreenFunc = okMakeGreenFunc(DataStruct,Dataset,Direction,FaultModel,Rake,Slip,Opening,GreenFuncName,varargin)
p = inputParser;
default_nu = 0.25;
addParameter(p,'nu',default_nu, @(x) isnumeric(x));

parse(p, varargin{:});
nu = p.Results.method;

GreenFunc = FaultModel;
if strcmp(Direction,'LOS')
    flag = 1;
elseif strcmp(Direction,'Azi')
    flag = 1;
elseif strcmp(Direction,'3D')
    flag = 0;
end

% Retrieve data (From Displacement)
ObsE = DataStruct.(Dataset).LocalX;
ObsN = DataStruct.(Dataset).LocalY;
if flag == 1
    Azimuth = DataStruct.(Dataset).Azimuth;
    Incidence = DataStruct.(Dataset).Incidence;
    LookDirection = DataStruct.LookDirection;
    % Calculate the quadrant the satellite is flying at
    if (mean(Azimuth(:),'omitnan')) > 0 && (mean(Azimuth(:),'omitnan')) < 90
        Quadrant = '1';
    elseif (mean(Azimuth(:),'omitnan')) > 90 && (mean(Azimuth(:),'omitnan')) < 180
        Quadrant = '4';
    elseif (mean(Azimuth(:),'omitnan')) > 180 && (mean(Azimuth(:),'omitnan')) < 270
        Quadrant = '3';
    elseif (mean(Azimuth(:),'omitnan')) > 270 && (mean(Azimuth(:),'omitnan')) < 360
        Quadrant = '2';
    end
end
% From FaultModel
okFault = FaultModel.okFault;

% Report the attributes
PatchCount = size(okFault,1);
disp(' ')
disp("******* Constructing Green's function okMakeGreenFunc.m *******")
if flag == 1
    disp(strcat('*** Convert Green func to:',32,Direction,32,'direction'))
end
disp(strcat('*** Patch count:',32,num2str(PatchCount)))
disp(strcat('*** Rake angle:',32,num2str(Rake)))
disp(strcat('*** Unit slip:',32,num2str(Slip)))
disp(strcat('*** Unit opening:',32,num2str(Opening)))
disp(strcat("*** Poisson's ratio:",32,num2str(nu)))
disp('*********************************')

% Start forwarding Okada
[Row,Col] = size(ObsE);
uE = zeros(Row*Col,PatchCount);
uN = zeros(Row*Col,PatchCount);
uZ = zeros(Row*Col,PatchCount);
digit = length(num2str(PatchCount));
progressStr = strcat('*** Processing patch: #',repmat('0',1,digit));
fprintf(progressStr)
for i = 1:PatchCount
    Xct = okFault(i,1);
    Yct = okFault(i,2);
    Zct = -okFault(i,3);
    Str = okFault(i,4);
    Dip = okFault(i,5);
    ASlength = okFault(i,6);
    ADwidth = okFault(i,7);

    % Execute okada85.m
    Num = strcat(repmat('0',1,digit-length(num2str(i))),num2str(i));
    progressStr = strcat('*** Processing patch: #',Num);
    fprintf([repmat('\b', 1, length(progressStr)), progressStr]);

    [Etmp,Ntmp,Ztmp] = okada85(ObsE-Xct,ObsN-Yct,Zct,Str,Dip,ASlength,ADwidth,Rake,Slip,Opening,nu);
    uE(:,i) = Etmp(:);
    uN(:,i) = Ntmp(:);
    uZ(:,i) = Ztmp(:);

end
fprintf('\n')

if flag == 1
    GreenLOS = zeros(Row*Col,PatchCount);
    GreenAzi = zeros(Row*Col,PatchCount);
    % Convert to LOS or Azi
    % Account for all types of LOS or Azimuth displacement directions
    if strcmp(LookDirection,'l') && strcmp(Quadrant,'1')
        GreenLOS = uE.*(-1).*cosd(Azimuth).*sind(Incidence) + uN.*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'r') && strcmp(Quadrant,'1')
        GreenLOS = uE.*cosd(Azimuth).*sind(Incidence) + uN.*(-1).*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'l') && strcmp(Quadrant,'2')
        GreenLOS = uE.*(-1).*cosd(Azimuth).*sind(Incidence) + uN.*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'r') && strcmp(Quadrant,'2')
        GreenLOS = uE.*cosd(Azimuth).*sind(Incidence) + uN.*(-1).*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'l') && strcmp(Quadrant,'3')
        GreenLOS = uE.*(-1).*cosd(Azimuth).*sind(Incidence) + uN.*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'r') && strcmp(Quadrant,'3')
        GreenLOS = uE.*cosd(Azimuth).*sind(Incidence) + uN.*(-1).*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'l') && strcmp(Quadrant,'4')
        GreenLOS = uE.*(-1).*cosd(Azimuth).*sind(Incidence) + uN.*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    elseif strcmp(LookDirection,'r') && strcmp(Quadrant,'4')
        GreenLOS = uE.*cosd(Azimuth).*sind(Incidence) + uN.*(-1).*sind(Azimuth).*sind(Incidence) + uZ.*(-1).*cosd(Incidence);
        GreenAzi = uE.*sind(Azimuth) + uN.*cosd(Azimuth);
    end
    if strcmp(Direction,'LOS')
        GreenFunc.(GreenFuncName) = GreenLOS;
    elseif strcmp(Direction,'Azi')
        GreenFunc.(GreenFuncName) = GreenAzi;
    end

else
    % No converting
    GreenFunc.(GreenFuncName).E = uE;
    GreenFunc.(GreenFuncName).N = uN;
    GreenFunc.(GreenFuncName).U = uZ;

end


end




