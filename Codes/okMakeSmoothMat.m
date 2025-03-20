%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.03                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***            okMakeSmoothMat.m            ***             %
%             ***********************************************             %
%                                                                         %
% Make smoothing matrix from the fault model generated from               %
% okMakeFaultModel.m leveraging okada85.m.                                %
% Currently only supports simple 2D Laplacian smoother                    %
%                                                                         %
% Each i and j is the linear indexing of the fault patches starting from  %
% the top to the bottm                                                    %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. FaultModel: Fault geometry from okMakeFaultModel.m                   %
%                                                                         %
% Example:                                                                %
% FaultModel = okMakeSmoothMat(FaultModel)                                %
%                                                                         %
% Output:                                                                 %
% 1. S: Smoothing matrix with dimension MxM                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SmoothModel = okMakeSmoothMat(FaultModel)
% Retrieve data
okFault = FaultModel.okFault;

% Smoothing function
smFunc = 'Laplacian';

% Find the along-strike and along-dip patch counts
ASpatch = FaultModel.PatchCountStrike;
ADpatch = FaultModel.PatchCountDip;


disp(' ')
disp("******* Constructing smoothing matrix okMakeSmoothMat.m *******")
disp(strcat('*** Along-strike patch count:',32,num2str(ASpatch)))
disp(strcat('*** Along-dip patch count:',32,num2str(ADpatch)))
disp(strcat('*** Smoothing function:',32,smFunc))


Dim = okFault(end,8);
S = zeros(Dim);
for i = 1:Dim
    % From Top-left parch to bottom-right parch
    ui = i-ASpatch;
    di = i+ASpatch;
    li = i-1;
    ri = i+1;

    S(i,i) = -4;
    % Detect edge cases, corners first
    if i == 1
    % 1. top-left corner
        S(i,ri) = 2;
        S(i,di) = 2;
    elseif i == ASpatch
    % 2. top-right corner
        S(i,li) = 2;
        S(i,di) = 2;
    elseif i == ASpatch*(ADpatch-1)
    % 3. bottom-left corner
        S(i,ui) = 2;
        S(i,ri) = 2;
    elseif i == ASpatch*ADpatch
    % 4. bottom-right corner
        S(i,ui) = 2;
        S(i,li) = 2;
    elseif i < ASpatch
    % 5. top row
        S(i,li) = 1;
        S(i,ri) = 1;
        S(i,di) = 2;
    elseif i > ASpatch*(ADpatch-1)
    % 6. bottom row
        S(i,li) = 1;
        S(i,ri) = 1;
        S(i,ui) = 2;
    elseif mod(i,ASpatch) == 1
    % 7. left column
        S(i,ui) = 1;
        S(i,di) = 1;
        S(i,ri) = 2;
    elseif mod(i,ASpatch) == 0
    % 8. right column
        S(i,ui) = 1;
        S(i,di) = 1;
        S(i,li) = 2;
    else
    % 9. inner patches
        S(i,li) = 1;
        S(i,ri) = 1;
        S(i,ui) = 1;
        S(i,di) = 1;
    end
    
end

% Output
SmoothModel = FaultModel;
SmoothModel.SmoothMat = S;

end



