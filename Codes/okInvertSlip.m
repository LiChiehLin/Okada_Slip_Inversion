%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.09                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***            okInvertSlip.m               ***             %
%             ***********************************************             %
%                                                                         %
% (Update: 2025.03.26)                                                    %
%   Allowing putting weights on data. Users can use the covariance matrix %
%   made in okInSARCov.m and okMakeInSARCovMat.m                          %
%                                                                         %
% Invert fault slip based on the input Displacement, Green's function and %
% smoothing matrix. This function also searches the best smoothing param. %
%                                                                         %
% There should only be 1 smoothing matrix in 1 fault model. May be updated%
% to include more sophisticated smoothing matrices per user's need        %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure or cell of structure. Containing the displ.    %
%    matrix and its attributes                                            %
% 2. Dataset: Character or cell of character. The Field that is going to  %
%    be used in DataStruct. Be paired with DataStruct                     %
% 3. FaultModel: Strucuture. The fault geomtery data strucuture           %
% 4. GreenFuncDataset: Character or cell of character. Green's function   %
%    field names. Be paired with DataStruct                               %
% 5. GreenFuncPosition: Matrix or vector. Specify how the Green's         %
%    function should be arranged in the design matrix (matrix indices)    %
%    This is used in 2 scenarios: (See below and Examples)                %
%    5.1. Solve for fault slip with an uniform rake angle                 %
%    5.2. Solve for strike-slip and dip-slip simultaneously               % 
% 6. SmoothMatDataset: Character or cell of character. Smoothing matrix   %
%    field names                                                          %
% 7. varargin (Name-value pair)                                           %
%    7.1. 'solver': 'lsq' or 'nnlsq' (default: lsq)                       %
%         'lsq': Least-squares                                            %
%         'nnlsq': Non-negative least squares                             %
%    7.2. 'smoothsearch': Vector for the smoothing constant search        %
%         (default: 1e-10~1e4)                                            %
%    7.3. 'weight': Characters or matrix (default: uniform weighting)     %
%                                                                         %
% Example: (1 diplacement, 1 fault, uniform rake, Least squares)          %
% ModelSlip = okInvertSlip(DataStruct,'Dsample', ...                      %
%             FaultModel,'GreenLOS',[11],'SmoothMat', ...                 %
%             'solver','lsq')                                             %
%                                                                         %
% Example: (1 diplacement, 1 fault, strike-slip and dip-slip, Non-neg     %
%          Least Squares)                                                 %
% ModelSlip = okInvertSlip(DataStruct,'Dsample', ...                      %
%             FaultModel,{'GreenDS','GreenSS'},[11,12],'SmoothMat', ...   %
%             'solver','nnlsq')                                           %
%                                                                         %
% Example: (1 diplacement, 1 fault, strike-slip and dip-slip, Non-neg     %
%          Least Squares, data weighted by covariance)                    %
% ModelSlip = okInvertSlip(DataStruct,'Dsample', ...                      %
%             FaultModel,{'GreenDS','GreenSS'},[11,12],'SmoothMat', ...   %
%             'solver','nnlsq',...                                        %
%             'weight','Covariance')                                      %
%                                                                         %
% Example: (2 diplacement, 1 fault, strike-slip and dip-slip, Non-neg     %
%          Least Squares)                                                 %
% ModelSlip = okInvertSlip({Data1,Data2},{'Dsample','Dsample'}, ...       %
%             FaultModel,{'GreenD1DS','GreenD1SS','GreenD2DS','GreenD2SS'}%
%             ,[11,12,21,22],'SmoothMat','solver','nnlsq')                %
%                                                                         %
% Example: (2 diplacement, 1 fault, strike-slip and dip-slip, Non-neg     %
%          Least Squares, customized smoothing parameter search)          %
% ModelSlip = okInvertSlip({Data1,Data2},{'Dsample','Dsample'}, ...       %
%             FaultModel,{'GreenD1DS','GreenD1SS','GreenD2DS','GreenD2SS'}%
%             ,[11,12,21,22],'SmoothMat','solver','nnlsq', ...            %
%             'smoothsearch',[1e-5,1e-4,1e-3,1e-2,1e-1,0,1,10])           %
%                                                                         %
% Note that: To solve for strike-slip and dip-slip, dip-slip Green's      %
%   function should be put in the first column and strike-slip second.    %
%   Otherwise, the rake angle will be wrong.                              %
%                                                                         %
% Output:                                                                 %
% 1. ModelSlip: Structure. Slip inversion output                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ModelSlip = okInvertSlip(DataStruct,Dataset,FaultModel,GreenFuncDataset,GreenFuncPosition,SmoothMatDataset,varargin)
p = inputParser;
% Parse options:
% For solver, choose:
% 1. Least-Squares (lsq, default)
% 2. Non-negative least squares (nnlsq)
default_solver = 'lsq';
% For smoothing parameter search
n = -10:0.1:4;
default_smoothsearch = transpose(10.^n);
% For weighting
default_weight = [];
addParameter(p, 'solver', default_solver, @(x) ischar(x) && strcmp(x,'lsq') || strcmp(x,'nnlsq'));
addParameter(p, 'smoothsearch', default_smoothsearch, @(x) isnumeric(x) && isvector(x));
addParameter(p, 'weight', default_weight, @(x) ischar(x) || iscell(x) || ismatrix(x));

parse(p, varargin{:});
solver = p.Results.solver;
SmoothParam = p.Results.smoothsearch;
WeightParse = p.Results.weight;


% Take care of the smoothing constant search
if size(SmoothParam,1) < size(SmoothParam,2)
    SmoothParam = transpose(SmoothParam);
end


% Get the arrangement of the design matrix
% Arrange displacements and corresponding Green's functions
tmp = max(GreenFuncPosition);
GR = floor(tmp/10);
GC = mod(tmp,10);
N = length(DataStruct); 
Displ = cell(N,1);
LocalX = cell(N,1);
LocalY = cell(N,1);
Weight = cell(N,1);
G = cell(GR,GC);
DisplN = zeros(N,1);
% If the input Green's function datasets only occupy 1 column, then it will
% assume this is only solving for 1 slip with an uniform rake angle. The
% slip model output will not include the rake and total slip fields
if GC == 2
    flagslip2 = 1;
else
    flagslip2 = 0;
end

for i = 1:N
    if N == 1
        tmpStruct = DataStruct;
        tmpDisplDataset = Dataset;
    else
        tmpStruct = DataStruct{i};
        tmpDisplDataset = Dataset{i};
    end

    Displtmp = tmpStruct.(tmpDisplDataset).Displ;
    LocalXtmp = tmpStruct.(tmpDisplDataset).LocalX;
    LocalYtmp = tmpStruct.(tmpDisplDataset).LocalY;
    LonOrigin = tmpStruct.(tmpDisplDataset).LongitudeOrigin;
    

    % Remove NaN values
    rmNaN = find(~isnan(Displtmp));
    Displ{i} = reshape(Displtmp(rmNaN),[],1);
    LocalX{i} = reshape(LocalXtmp(rmNaN),[],1);
    LocalY{i} = reshape(LocalYtmp(rmNaN),[],1);
    

    % Link displacement with its associated Green's function
    LinkInd = find(i == floor(GreenFuncPosition./10));
    for j = LinkInd
        GInd = GreenFuncPosition(j);
        Grow = floor(GInd/10);
        Gcol = mod(GInd,10);
        tmpGreenDataset = GreenFuncDataset{j};
        Gtmp = FaultModel.(tmpGreenDataset);  
        Gtmp = Gtmp(rmNaN,:);

        % Store in pre-allocated cell
        G{Grow,Gcol} = Gtmp;
    end

    % Deal with the weighting 
    if ~isempty(WeightParse)
        if N == 1 && ischar(WeightParse) && ~isempty(WeightParse)
            tmpWeightDataset = WeightParse;
        elseif N ~= 1 && iscell(WeightParse) && ~isempty(WeightParse)
            tmpWeightDataset = WeightParse{i};
        end

    % Retrieve data
    Weighttmp = tmpStruct.(tmpDisplDataset).(tmpWeightDataset);

    % Remove NaN values
    Weight{i} = Weighttmp(rmNaN,:);

    % Store data count
    DisplN(i) = size(Weighttmp,1);
    end

    
end
% Unpack the arranged stuff
Displ = cell2mat(Displ);
LocalX = cell2mat(LocalX);
LocalY = cell2mat(LocalY);
G = cell2mat(G);


% Unpack weighting matrix
if ~isempty(WeightParse)
    W = zeros(numel(Displ),numel(Displ));
    BeginI = 0;
    for i = 1:length(Weight)
        Wtmp = Weight{i};
        EndI = sum(DisplN(1:i));
        Ind = (BeginI+1):EndI;

        W(Ind,Ind) = Wtmp;
        BeginI = Ind(end);
    end
else
    W = eye(numel(Displ));

end


% Arrange smoothing matrix
% Assume solving for strike-slip and dip-slip
Stmp = FaultModel.(SmoothMatDataset);
SmoothSize = size(Stmp,1);
if SmoothSize ~= size(G,2)
    S = zeros(SmoothSize*2);
    BeginI = 0;
    EndI = SmoothSize*(1:2);
    for i = 1:size(G,2)/size(Stmp,2)
        Ind = (BeginI+1):EndI(i);
        S(Ind,Ind) = Stmp;
        BeginI = Ind(end);
    end
else
    S = Stmp;
end

% Get the patch sizes
% Assuming solving for strike-slip and dip-slip
ASPatchSize = FaultModel.AlongStrikePatchSize;
ADPatchSize = FaultModel.AlongDipPatchSize;
PatchSizetmp = ASPatchSize.*ADPatchSize;
PatchCount = FaultModel.TotalPatchCount;
TotalPatchCount = sum(PatchCount.*2);
M = length(PatchCount);
PatchSize = cell(M,1);
for i = 1:M
    PatchSize{i} = repmat(PatchSizetmp(i),PatchCount(i)*2,1);
end
PatchSize = cell2mat(PatchSize);



% Start inversion
obsN = size(G,1);
patM = size(S,1);
Lcurve = zeros(length(SmoothParam),2);
FaultSlip = zeros(patM+1,length(SmoothParam));
ModelPrediction = zeros(obsN,length(SmoothParam));
Moment = zeros(length(SmoothParam),1);
MomentMagnitude = zeros(length(SmoothParam),1);
ResolutionMat = cell(length(SmoothParam),1);
RMSE = zeros(length(SmoothParam),1);
TotalSlip = zeros(SmoothSize,length(SmoothParam));
Rake = zeros(SmoothSize,length(SmoothParam));
for i = 1:length(SmoothParam)
    % Inversion
    alpha = SmoothParam(i);
    Gtmp = inv(W)*G; % Put the weighting here
    A = [Gtmp,ones(obsN,1); alpha.*S,zeros(patM,1)];

    dtmp = inv(W)*Displ(:); % Put the weighting here 
    d = [dtmp;zeros(patM,1)];
    if strcmp(solver,'lsq')
        m = A\d;
    elseif strcmp(solver,'nnlsq')
        m = lsqnonneg(A,d);
    end
    % Calculate resolution matrix
    Resol = (A'*A)\(A'*A);

    % Calculate seismic moment and moment magnitude
    M0 = sum(m(1:sum(TotalPatchCount)).*PatchSize.*30*10^9);
    Mw = (2/3)*(log10(M0)-9.1);

    % Calculate the residual and prepare for L-curve
    tmp = [G,ones(obsN,1)]*m;
    ResNorm = sqrt(sum((inv(W)*(Displ(:) - tmp(1:obsN))).^2));
    SolNorm = sqrt(sum((S*m(1:sum(TotalPatchCount))).^2));
    Pred = tmp(1:obsN);

    % Report inversion result
    avgSlip = mean(m(1:sum(TotalPatchCount)));
    disp(strcat('* Average slip:',32,num2str(avgSlip),32,'| Smoothing param:',32,num2str(alpha), ...
        32,'| Solution res:',32,num2str(SolNorm),32,'| Model res.:',32,num2str(ResNorm)))

    % Output inversion result
    FaultSlip(:,i) = m;
    Lcurve(i,1) = SolNorm;
    Lcurve(i,2) = ResNorm;
    ModelPrediction(:,i) = Pred;
    RMSE(i) = sqrt(mean((Displ(:) - Pred).^2));
    Moment(i) = M0;
    MomentMagnitude(i) = Mw;
    ResolutionMat{i} = Resol;

    % Output total slip and rake angle if input 2 Green's function
    if flagslip2 == 1
        tmp = reshape(m(1:sum(TotalPatchCount)),SmoothSize,2);
        TotalSlip(:,i) = sqrt(tmp(:,1).^2 + tmp(:,2).^2);
        Rake(:,i) = atand(tmp(:,1)./tmp(:,2));
    end
end

% Output to a structure to make it more organized
% Propagate input Displacement structure to the output
ModelSlip.Displ = Displ;
ModelSlip.LocalX = LocalX;
ModelSlip.LocalY = LocalY;
ModelSlip.FaultSlip = FaultSlip;
ModelSlip.TotalSlip = TotalSlip;
ModelSlip.Rake = Rake;
ModelSlip.Lcurve = Lcurve;
ModelSlip.GreenFunc = G;
ModelSlip.SmoothMat = S;
ModelSlip.SmoothParam = SmoothParam;
ModelSlip.SourceDispl = Dataset;
ModelSlip.FaultModel = FaultModel;
ModelSlip.ModelPrediction = ModelPrediction;
ModelSlip.RMSE = RMSE;
ModelSlip.Moment = Moment;
ModelSlip.MomentMagnitude = MomentMagnitude;
ModelSlip.Solver = solver;
ModelSlip.LongitudeOrigin = LonOrigin;
ModelSlip.ResolutionMat = ResolutionMat;
ModelSlip.WeightMat = W;
end
