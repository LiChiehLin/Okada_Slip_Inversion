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
% (Update: 2025.09.12)                                                    %
%   1. Fixed some syntax and made it more concise (maybe?)                %
%   2. Support 'along-strike' combination                                 %
%   3. Support 'along-dip' combination                                    %
%   4. Allows different Strike or Dip patch counts                        %
%                                                                         %
% Combine connected fault segments or make multiple fault models into 1   %
% fault data structure. Later processing only requires 1 input of fault   %
% data strucuture, so use this routine to comebine them all               %
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
N = length(FaultModels);
okFault = cell(N,1);
Nodes = cell(N,1);
PatchX = cell(1,N);
PatchY = cell(1,N);
PatchZ = cell(1,N);
StartXYZ = cell(N,1);
Length = cell(N,1);
Width = cell(N,1);
Strike = cell(N,1);
Dip = cell(N,1);
Depth = cell(N,1);
ASPatchSize = cell(N,1);
ADPatchSize = cell(N,1);
PatchCountDip = cell(N,1);
PatchCountStrike = cell(N,1);
TotalPatchCount = cell(N,1);
RowcountMat = cell(N,1);
ColcountMat = cell(N,1);
if strcmp(Direction,'strike')
    PatchCount = 0;
    for i = 1:N
        PatchCount = PatchCount + size(FaultModels{i}.okFault,1);
        StartXYZ{i} = FaultModels{i}.StartXYZ;
        Length{i} = FaultModels{i}.Length;
        Width{i} = FaultModels{i}.Width;
        Strike{i} = FaultModels{i}.Strike;
        Dip{i} = FaultModels{i}.Dip;
        Depth{i} = FaultModels{i}.Depth;
        ASPatchSize{i} = FaultModels{i}.AlongStrikePatchSize;
        ADPatchSize{i} = FaultModels{i}.AlongDipPatchSize;
        PatchCountDip{i} = FaultModels{i}.PatchCountDip;
        PatchCountStrike{i} = FaultModels{i}.PatchCountStrike;
        TotalPatchCount{i} = FaultModels{i}.TotalPatchCount;
        if isfield(FaultModels{i},'RowcountMat')
            RowcountMat{i} = FaultModels{i}.RowcountMat;
        end
        if isfield(FaultModels{i},'ColcountMat')
            ColcountMat{i} = FaultModels{i}.ColcountMat;
        end

        % Populate the input fault model to a field
        FaultModelName = strcat('FaultModel',num2str(i));
        CombinedFaultModels.(FaultModelName) = FaultModels{i};
    end

    % Arrange them together
    okFault = zeros(PatchCount,8);
    Nodes = zeros(PatchCount,4);
    PatchX = zeros(4,PatchCount);
    PatchY = zeros(4,PatchCount);
    PatchZ = zeros(4,PatchCount);

    PCDip1 = PatchCountDip{1};
    PCDip2 = PatchCountDip{2};
    PCStrike1 = PatchCountStrike{1};
    PCStrike2 = PatchCountStrike{2};
    if isempty(RowcountMat{1}) % Combine 2 models with no prior combination
        if PCDip2 > PCDip1 % The second model has more dip patches than the first
            Iter = PCDip1;
            Rowcountstore = zeros(Iter,2);
            Rowmax = ceil(PCDip2/PCDip1);
            Rowlimit = floor(PCDip2/PCDip1) + 0.8;
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind2last = 0;
                    Rowcountcum = PCDip2/PCDip1;
                else
                    Rowcountcum = Rowcountcum + (PCDip2/PCDip1 - floor(PCDip2/PCDip1));
                end
    
                if mod(PCDip2,PCDip1) ~= 0
                    if Rowcountcum >= Rowlimit
                        Rowcount = Rowmax;
                        Rowlimit = Rowlimit + 1;
                    else
                        Rowcount = Rowmax - 1;
                    end
                else
                    Rowcount = PCDip2/PCDip1;
                end
    
                Ind1 = (i-1)*PCStrike1+1:i*PCStrike1; % Indexing first fault 
                Ind2 = Ind2last+1:(Ind2last+1 + PCStrike2*Rowcount - 1); % Indexing second fault
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2))); % Indices to put in
                CIndlast = CInd(end); % Populate to the next fault model indexing
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Ind2last = Ind2(end);
                Rowcountstore(i,:) = [1,Rowcount];
            end

        elseif PCDip2 < PCDip1 % The first model has more dip patches than the second
            Iter = PCDip2;
            Rowcountstore = zeros(Iter,2);
            Rowmax = ceil(PCDip1/PCDip2);
            Rowlimit = floor(PCDip1/PCDip2) + 0.8;
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind1last = 0;
                    Rowcountcum = PCDip1/PCDip2;
                else
                    Rowcountcum = Rowcountcum + (PCDip1/PCDip2 - floor(PCDip1/PCDip2));
                end
    
                if mod(PCDip1,PCDip2) ~= 0
                    if Rowcountcum >= Rowlimit
                        Rowcount = Rowmax;
                        Rowlimit = Rowlimit + 1;
                    else
                        Rowcount = Rowmax - 1;
                    end
                else
                    Rowcount = PCDip1/PCDip2;
                end
                
                Ind1 = Ind1last + (1:Rowcount*PCStrike1); % Indexing first fault
                Ind2 = (i-1)*PCStrike2+1:i*PCStrike2; % Indexing second fault 
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2))); % Indices to put in
                CIndlast = CInd(end); % Populate to the next fault model indexing
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Ind1last = Ind1(end);
                Rowcountstore(i,:) = [Rowcount,1];
            end

        else % Both have the same dip patches
            Iter = PCDip1;
            Rowcountstore = zeros(Iter,2);
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind1last = 0;
                    Ind2last = 0;
                end

                Ind1 = Ind1last + (1:PatchCountStrike{1});
                Ind2 = Ind2last + (1:PatchCountStrike{2});
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2)));
                CIndlast = CInd(end);
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Ind1last = Ind1(end);
                Ind2last = Ind2(end);
                Rowcountstore(i,:) = [1,1];
            end
        end

    else % Combine 2 models with prior combination (Fault1 is prior combined and Fault2 is not)
        RCDip1 = size(RowcountMat{1},1);
        if PCDip2 > RCDip1 % The second model has more dip patches than the combined
            Iter = RCDip1;
            Rowcountstore = zeros(RCDip1,1);
            Rowmax = ceil(PCDip2/RCDip1);
            Rowlimit = floor(PCDip2/RCDip1) + 0.8;
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind1last = 0;
                    Ind2last = 0;
                    Rowcountcum = PCDip2/RCDip1;
                else
                    Rowcountcum = Rowcountcum + (PCDip2/RCDip1 - floor(PCDip2/RCDip1));
                end

                if mod(PCDip2,RCDip1) ~= 0
                    if Rowcountcum >= Rowlimit
                        Row2count = Rowmax;
                        Rowlimit = Rowlimit + 1;
                    else
                        Row2count = Rowmax - 1;
                    end
                else
                    Row2count = PCDip2/RCDip1;
                end

                Row1num = i;
                Ind1 = Ind1last + (1:sum(RowcountMat{1}(Row1num,:)*PatchCountStrike{1})); % Indexing first fault 
                Ind2 = Ind2last + (1:Row2count*PatchCountStrike{2}); % Indexing second fault
                if i == Iter && ~any(Ind2 == size(FaultModels{2}.okFault,1))
                    % This means by the end of the last iteration, it still
                    % hasn't exhausted the whole list, thus Row2count + 1
                    Row2count = Row2count + 1;
                    Ind2 = Ind2last + (1:Row2count*PatchCountStrike{2});
                end
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2)));
                CIndlast = CInd(end);
                Ind1last = Ind1(end);
                Ind2last = Ind2(end);
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Rowcountstore(i) = Row2count;
            end
            Rowcountstore = [RowcountMat{1},Rowcountstore];

        elseif PCDip2 < RCDip1 % The combined model has more dip patches than the second
            Iter = PCDip2;
            Rowcountstore = zeros(RCDip1,1);
            Rowmax = ceil(RCDip1/PCDip2);
            Rowlimit = floor(RCDip1/PCDip2) + 0.8;
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind1last = 0;
                    Ind2last = 0;
                    Row1numlast = 0;
                    Rowcountcum = RCDip1/PCDip2;
                else
                    Rowcountcum = Rowcountcum + (RCDip1/PCDip2 - floor(RCDip1/PCDip2));
                end

                if mod(RCDip1,PCDip2) ~= 0
                    if Rowcountcum >= Rowlimit
                        Row1count = Rowmax;
                        Rowlimit = Rowlimit + 1;
                    else
                        Row1count = Rowmax - 1;
                    end
                else
                    Row1count = RCDip1/PCDip2;
                end

                Row1num = (Row1numlast+1):((Row1numlast+1)+Row1count-1);
                Ind1 = Ind1last + (1:sum(RowcountMat{1}(Row1num,:)*PatchCountStrike{1})); % Indexing first fault 
                Ind2 = Ind2last + (1:PatchCountStrike{2}); % Indexing second fault
                if i == Iter && ~any(Ind1 == size(FaultModels{1}.okFault,1))
                    % This means by the end of the last iteration, it still
                    % hasn't exhausted the whole list, thus Row2count + 1
                    Row1count = Row1count + 1;
                    Row1num = i:(i+Row1count-1);
                    Ind1 = Ind1last + (1:sum(RowcountMat{1}(Row1num,:)*PatchCountStrike{1}));
                end
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2)));
                CIndlast = CInd(end);
                Ind1last = Ind1(end);
                Ind2last = Ind2(end);
                Row1numlast = Row1num(end);
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Rowcountstore(i) = Row1count;
            end
            Rowcountstore = [RowcountMat{1},Rowcountstore];

        else % Both have the same dip patches
            Iter = RCDip1;
            Rowcountstore = zeros(RCDip1,1);
            for i = 1:Iter
                if i == 1
                    CIndlast = 0;
                    Ind1last = 0;
                    Ind2last = 0;
                end

                Ind1 = Ind1last + (1:sum(RowcountMat{1}(i,:)*PatchCountStrike{1}));
                Ind2 = Ind2last + (1:PatchCountStrike{2});
                CInd = CIndlast + (1:(length(Ind1)+length(Ind2)));
                CIndlast = CInd(end);
                Ind1last = Ind1(end);
                Ind2last = Ind2(end);
                okFault(CInd,:) = [FaultModels{1}.okFault(Ind1,:);FaultModels{2}.okFault(Ind2,:)];
                Nodes(CInd,:) = [FaultModels{1}.Nodes(Ind1,:);FaultModels{2}.Nodes(Ind2,:)];
                PatchX(:,CInd) = [FaultModels{1}.PatchX(:,Ind1),FaultModels{2}.PatchX(:,Ind2)];
                PatchY(:,CInd) = [FaultModels{1}.PatchY(:,Ind1),FaultModels{2}.PatchY(:,Ind2)];
                PatchZ(:,CInd) = [FaultModels{1}.PatchZ(:,Ind1),FaultModels{2}.PatchZ(:,Ind2)];
                Rowcountstore(i) = 1;
            end
            Rowcountstore = [RowcountMat{1},Rowcountstore];
        end
    end
    Nodes(:,end) = 1:size(Nodes,1);
    okFault(:,end) = 1:size(okFault,1);
    % Put them back to cell so I don't have to change them in the end
    okFault = mat2cell(okFault,PatchCount);
    Nodes = mat2cell(Nodes,PatchCount);
    PatchX = mat2cell(PatchX,4);
    PatchY = mat2cell(PatchY,4);
    PatchZ = mat2cell(PatchZ,4);


elseif strcmp(Direction,'dip')
    PatchCount = 0;
    for i = 1:N
        PatchCount = PatchCount + size(FaultModels{i}.okFault,1);
        okFault{i} = FaultModels{i}.okFault;
        Nodes{i} = FaultModels{i}.Nodes;
        PatchX{i} = FaultModels{i}.PatchX;
        PatchY{i} = FaultModels{i}.PatchY;
        PatchZ{i} = FaultModels{i}.PatchZ;
        StartXYZ{i} = FaultModels{i}.StartXYZ;
        Length{i} = FaultModels{i}.Length;
        Width{i} = FaultModels{i}.Width;
        Strike{i} = FaultModels{i}.Strike;
        Dip{i} = FaultModels{i}.Dip;
        Depth{i} = FaultModels{i}.Depth;
        ASPatchSize{i} = FaultModels{i}.AlongStrikePatchSize;
        ADPatchSize{i} = FaultModels{i}.AlongDipPatchSize;
        PatchCountDip{i} = FaultModels{i}.PatchCountDip;
        PatchCountStrike{i} = FaultModels{i}.PatchCountStrike;
        TotalPatchCount{i} = FaultModels{i}.TotalPatchCount;
        if isfield(FaultModels{i},'ColcountMat')
            ColcountMat{i} = FaultModels{i}.ColcountMat;
        end

        % Populate the input fault model to a field
        FaultModelName = strcat('FaultModel',num2str(i));
        CombinedFaultModels.(FaultModelName) = FaultModels{i};
    end

else
    error_message = ['Direction can either be:\n' ...
        '"strike" or \n' ...
        '"dip" or \n' ...
        '"separate"'];
    error('u:stuffed:it',error_message)
end
% Output the combined results
CombinedFaultModels.okFault = cell2mat(okFault);
CombinedFaultModels.okFault(:,end) = 1:size(cell2mat(okFault),1);
CombinedFaultModels.Nodes = cell2mat(Nodes);
CombinedFaultModels.Nodes(:,end) = 1:size(cell2mat(Nodes),1);
CombinedFaultModels.PatchX = cell2mat(PatchX);
CombinedFaultModels.PatchY = cell2mat(PatchY);
CombinedFaultModels.PatchZ = cell2mat(PatchZ);
CombinedFaultModels.StartXYZ = cell2mat(StartXYZ);
CombinedFaultModels.Length = cell2mat(Length);
CombinedFaultModels.Width = cell2mat(Width);
CombinedFaultModels.Strike = cell2mat(Strike);
CombinedFaultModels.Dip = cell2mat(Dip);
CombinedFaultModels.Depth = cell2mat(Depth);
CombinedFaultModels.AlongStrikePatchSize = cell2mat(ASPatchSize);
CombinedFaultModels.AlongDipPatchSize = cell2mat(ADPatchSize);
CombinedFaultModels.PatchCountStrike = cell2mat(PatchCountStrike);
CombinedFaultModels.PatchCountDip = cell2mat(PatchCountDip);
CombinedFaultModels.FaultModelCount = N;
CombinedFaultModels.TotalPatchCount = cell2mat(TotalPatchCount);
if exist('Rowcountstore','var')
    CombinedFaultModels.RowcountMat = Rowcountstore;
end
if exist('Colcountstore','var')
    CombinedFaultModels.ColcountMat = Colcountstore;
end
end
