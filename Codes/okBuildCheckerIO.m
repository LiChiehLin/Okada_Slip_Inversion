%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2026.01.12                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***           okBuildCheckerIO.m            ***             %
%             ***********************************************             %
%                                                                         %
% Construct the checker board style of slip/opening I/O for generating    %
% Green's function or searching for optimal fault model                   %
% Routines that use this sub-routine:                                     %
%   - okMakeGreenFunc.m                                                   %
%   - okSearchOptimalFaultModel.m                                         %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
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
function Checker = okBuildCheckerIO(FaultModel)
% Retrieve some fault model attributes
PatchCountStrike = FaultModel.PatchCountStrike;
PatchCountDip = FaultModel.PatchCountDip;
TotalPatchCount = sum(FaultModel.TotalPatchCount);

% Get the Starting and ending indices of each layer
EndI = cumsum(PatchCountStrike.*PatchCountDip);
StartI = [1;EndI+1];
StartI = StartI(1:end-1);

% Determine the checker for each layer
Checker = zeros(TotalPatchCount,1);
for j = 1:length(StartI)
    SInd = StartI(j);
    EInd = EndI(j);

    if j == 1
        CheckInd = SInd:2:EInd;
        Checker(CheckInd) = 1;

        CheckIndpre = CheckInd;
        SIndpre = SInd;
    else
        if any(CheckIndpre == SIndpre)
            % Push the checker board 1 patch ahead
            CheckInd = (SInd+1):2:EInd;
            Checker(CheckInd) = 1;
        else
            CheckInd = SInd:2:EInd;
            Checker(CheckInd) = 1;
        end
        CheckIndpre = CheckInd;
        SIndpre = SInd;
    end
end

end