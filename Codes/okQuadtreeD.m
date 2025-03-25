%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.06                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okQuadtreeD.m              ***             %
%             ***********************************************             %
%                                                                         %
% Perform quadtree downsampling based on the sampling criterion set in the%
% input parameter file made from okMakeDsampleParam.m                     %
%                                                                         %
% *Local mean removed ('variance'):                                       %
% Jonsson, S. (2002). Modeling volcano and earthquake deformation from    %
% satellite radar interferometric observations. Stanford University.      %
%                                                                         %
% *Bilinear ramp removed ('curvature'):                                   %
% Simons, M., Fialko, Y., & Rivera, L. (2002). Coseismic deformation from %
% the 1999 M w 7.1 Hector Mine, California, earthquake as inferred from   %
% InSAR and GPS observations. Bulletin of the Seismological Society of    %
% America, 92(4), 1390-1402.                                              %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. Data: Structure. Containing the displ. matrix and its attributes     %
% [2-10]: Parsed from okDsample.m                                         %
% 2. rstart: Numeric. Starting row                                        %
% 3. rend: Numeric. Ending row                                            %
% 4. cstart: Numeric. Starting column                                     %
% 5. cend: Numeric. Ending column                                         %
% 6. Criterion: String. Downsampling criterion. Choose from:              %
%   'variance': Method used in (Jonsson, 2002; see above)                 %
%   'curvature': Method used in (Simons et al., 2002; see above)          %
% 7. Tolerance: Numeric. Downsampling criterion threshold                 %
% 8. LeafMinSize: Numeric. Minimum allowed leaf size                      %
% 9. LeafMaxSize: Numeric. Maximum allowed leaf size                      %
% 10. NanPixelAllow: Numeric. Percentage of allowed NaN pixels in a leaf  %
%                                                                         %
% Output:                                                                 %
% 1. LeafCoord: Matrix. Each leaf's starting/ending rows and columns      %
%    col1: starting row                                                   %
%    col2: ending row                                                     %
%    col3: starting column                                                %
%    col4: ending column                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LeafCoord = okQuadtreeD(Data,rstart,rend,cstart,cend,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow)

Tol = Tolerance;
% Calculate the threshold of interest 
if strcmp(Criterion,'variance')
    % From (Jonsson, 2002)
    % Remove the local mean before calculating the variance
    LocalMean = mean(Data(rstart:rend,cstart:cend),'all','omitnan');
    ThreshVal = var(Data(rstart:rend,cstart:cend)-LocalMean,0,'all','omitnan');
elseif strcmp(Criterion,'curvature')
    % From (Simons et al., 2002)
    % Remove the local ramp before calculating the variance
    rMat = repmat(flipud(transpose(rstart:rend)),1,length(cstart:cend));
    cMat = repmat(cstart:cend,length(rstart:rend),1);
    % Model the ramp
    Gramp = [cMat(:),rMat(:),ones(size(cMat,1)*size(cMat,2),1)];
    tmp = Data(rstart:rend,cstart:cend);
    dramp = tmp(:);
    rmNaN = find(~isnan(dramp));
    % Remove NaN
    dramp = dramp(rmNaN);
    Gramp = Gramp(rmNaN,:);
    mramp = Gramp\dramp;
    RampVal = Gramp*mramp;

    ThreshVal = var(dramp-RampVal,0,'all','omitnan');
end

% Calculate the current leaf size
LeafSize = (rend-rstart+1)*(cend-cstart+1);
% Calculate the NaN pixel counts
NaNPixelRatio = sum(sum(isnan(Data(rstart:rend,cstart:cend)))) / (size(Data(rstart:rend,cstart:cend),1)*size(Data(rstart:rend,cstart:cend),2));

% Base case: 3 scenarios to stop except for 1
% 1. Criterion < tolerance
% 2. Leaf size < minimum leaf size allowed
% 3. NaN pixel ratio > NaN pixel ratio allowed
% when current leaf size is still larger than maximum leaf size allowed,
% then keeps dividing
if ThreshVal < Tol || LeafSize < LeafMinSize || NaNPixelRatio > NaNPixelAllow
    if LeafSize < LeafMaxSize
        LeafCoord = [rstart,rend,cstart,cend];
        return;
    end
end

% Divide the current matrix into four quadrants
rmid = floor((rstart + rend)/2);
cmid = floor((cstart + cend)/2);

% Call the recursive function if criterion is not met
% Q1
Q1 = okQuadtreeD(Data,rstart,rmid,cmid,cend,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow);
% Q2 
Q2 = okQuadtreeD(Data,rstart,rmid,cstart,cmid,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow);
% Q3
Q3 = okQuadtreeD(Data,rmid,rend,cstart,cmid,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow);
% Q4
Q4 = okQuadtreeD(Data,rmid,rend,cmid,cend,Criterion,Tolerance,LeafMinSize,LeafMaxSize,NaNPixelAllow);

% Output the coordinates
LeafCoord = [Q1;Q2;Q3;Q4];
end











