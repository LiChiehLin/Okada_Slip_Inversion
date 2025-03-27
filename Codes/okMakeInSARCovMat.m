%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.26                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***           okMakeInSARCovMat.m           ***             %
%             ***********************************************             %
%                                                                         %
% Make the covariance matrix based on the autocorrelation result given the%
% function fit                                                            %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement, and   %
%    the local coordinates                                                %
% 3. func: Character. The fitted function when making auto-correlation    %
%                                                                         %
% Example:                                                                %
% DataStruct = okMakeInSARCovMat(DataStruct,'Dsample',...                 %
%       'exp1')                                                           %
%                                                                         %
% Output:                                                                 %
% Covariance: Structure. Contains the covariance matrix                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Covariance = okMakeInSARCovMat(DataStruct,Dataset,func)
Covariance = DataStruct;
% Retrieve data
LocalX = DataStruct.(Dataset).LocalX;
LocalY = DataStruct.(Dataset).LocalY;
Function = DataStruct.(Dataset).(func);
N = numel(LocalX);

disp(strcat('*** Making covariance matrix based on function:',32,func))

WeightMat = zeros(N,N);
for i = 1:N
    x = LocalX(i);
    y = LocalY(i);
    Dist = sqrt((LocalX - x).^2 + (LocalY - y).^2);
    Corr = Function(Dist);

    WeightMat(i,:) = Corr;
end

% Output
Covariance.(Dataset).Covariance = WeightMat;

end




