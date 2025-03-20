%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.12                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***         okCombineFaultModel.m           ***             %
%             ***********************************************             %
%                                                                         %
% Combine connected fault segments or make multiple fault models into 1   %
% fault data structure. Later processing only requires 1 input of fault   %
% data strucuture, so use this routine to comebine them all               %
%                                                                         %
% Current version only supports connecting faults at depth (along-dip) as %
% along-strike fault segments might intersect each other and gets tricky  %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. FaultModels. Cell. Contain the fault data structure                  %
% 2. Direction. Character. Specify which direction the fault is connected %
%    'strike': Along-strike (Currently not supported)                     %
%    'dip': Along-dip (e.g. Listric fault model)                          %
%    'separate': A secondary fault that is connected with another         %
%                                                                         %
% Example 1: Listric fault                                                %
% CombineFault = okCombineFaultModel({FaultModel1,FaultModel2},'dip')     %
% Example 2: Combine a non-connected secondary strucuture                 %
% CombineFault = okCombineFaultModel({FaultModel1,FaultModel2},'separate')%
%                                                                         %
% Output:                                                                 %
% 1. CombinedFaultModels: Structure. Combined fault data structure        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function CombinedFaultModels = okCombineFaultModel(FaultModels,Direction)
% Either connecting along-strike or along-dip
% Right now, only along-dip
N = length(FaultModels);
if strcmp(Direction,'strike')
    % To-be updated
    % Mind the intersecting fault segments, raise warning

elseif strcmp(Direction,'dip')
    okFault = cell(N,1);
    Nodes = cell(N,1);
    PatchX = cell(1,N);
    PatchY = cell(1,N);
    PatchZ = cell(1,N);
    PatchCountDip = 0;
    ASPatchSize = cell(N,1);
    ADPatchSize = cell(N,1);
    TotalPatchCount = cell(N,1);
    for i = 1:N
        % Retrieve data
        okFault{i} = FaultModels{i}.okFault;
        Nodes{i} = FaultModels{i}.Nodes;
        PatchX{i} = FaultModels{i}.PatchX;
        PatchY{i} = FaultModels{i}.PatchY;
        PatchZ{i} = FaultModels{i}.PatchZ;
        ASPatchSize{i} = FaultModels{i}.AlongStrikePatchSize;
        ADPatchSize{i} = FaultModels{i}.AlongDipPatchSize;
        PatchCountDip = PatchCountDip + FaultModels{i}.PatchCountDip;
        PatchCountStrike = FaultModels{i}.PatchCountStrike;
        TotalPatchCount{i} = FaultModels{i}.TotalPatchCount;

        % Populate the input fault model to a field
        FaultModelName = strcat('FaultModel',num2str(i));
        CombinedFaultModels.(FaultModelName) = FaultModels{i};
    end
    % Output the rest
    CombinedFaultModels.okFault = cell2mat(okFault);
    CombinedFaultModels.okFault(:,end) = 1:size(cell2mat(okFault),1);
    CombinedFaultModels.Nodes = cell2mat(Nodes);
    CombinedFaultModels.PatchX = cell2mat(PatchX);
    CombinedFaultModels.PatchY = cell2mat(PatchY);
    CombinedFaultModels.PatchZ = cell2mat(PatchZ);
    CombinedFaultModels.AlongStrikePatchSize = cell2mat(ASPatchSize);
    CombinedFaultModels.AlongDipPatchSize = cell2mat(ADPatchSize);
    CombinedFaultModels.PatchCountStrike = PatchCountStrike;
    CombinedFaultModels.PatchCountDip = PatchCountDip;
    CombinedFaultModels.FaultModelCount = N;
    CombinedFaultModels.TotalPatchCount = cell2mat(TotalPatchCount);

elseif strcmp(Direction,'separate')
    okFault = cell(N,1);
    Nodes = cell(N,1);
    PatchX = cell(1,N);
    PatchY = cell(1,N);
    PatchZ = cell(1,N);
    PatchCount = zeros(1,N);
    ASPatchSize = cell(N,1);
    ADPatchSize = cell(N,1);
    TotalPatchCount = cell(N,1);
    for i = 1:N
        % Retrieve data
        okFault{i} = FaultModels{i}.okFault;
        Nodes{i} = FaultModels{i}.Nodes;
        PatchX{i} = FaultModels{i}.PatchX;
        PatchY{i} = FaultModels{i}.PatchY;
        PatchZ{i} = FaultModels{i}.PatchZ;
        ASPatchSize{i} = FaultModels{i}.AlongStrikePatchSize;
        ADPatchSize{i} = FaultModels{i}.AlongDipPatchSize;
        PatchCount(i) = size(FaultModels{i}.SmoothMat,1);
        TotalPatchCount{i} = FaultModels{i}.TotalPatchCount;

        % Populate the input fault model to a field
        FaultModelName = strcat('FaultModel',num2str(i));
        CombinedFaultModels.(FaultModelName) = FaultModels{i};
    end

    % Arrange smoothing matrix

    SmoothMat = zeros(sum(PatchCount));
    BeginI = 0;
    EndI = cumsum(PatchCount);
    for i = 1:N
        Ind = (BeginI+1):EndI(i);
        SmoothMat(Ind,Ind) = FaultModels{i}.SmoothMat;
        BeginI = Ind(end);
    end
    % Output the rest
    CombinedFaultModels.okFault = cell2mat(okFault);
    CombinedFaultModels.okFault(:,end) = 1:size(cell2mat(okFault),1);
    CombinedFaultModels.Nodes = cell2mat(Nodes);
    CombinedFaultModels.PatchX = cell2mat(PatchX);
    CombinedFaultModels.PatchY = cell2mat(PatchY);
    CombinedFaultModels.PatchZ = cell2mat(PatchZ);
    CombinedFaultModels.AlongStrikePatchSize = cell2mat(ASPatchSize);
    CombinedFaultModels.AlongDipPatchSize = cell2mat(ADPatchSize);
    CombinedFaultModels.FaultModelCount = N;
    CombinedFaultModels.SmoothMat = SmoothMat;
    CombinedFaultModels.TotalPatchCount = cell2mat(TotalPatchCount);

else
    error_message = ['Direction can either be:\n' ...
        '"strike" or \n' ...
        '"dip"'];
    error('u:stuffed:it',error_message)
end

end
