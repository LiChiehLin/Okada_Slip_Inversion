%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.25                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okInSARCov.m               ***             %
%             ***********************************************             %
%                                                                         %
% Perform auto-correlation to calculate the covariance from pixel-to-pixel%
% This is further used to produce covariance matrix for weighting InSAR   %
% measurements during inversion                                           %
%                                                                         %
% Next steps:                                                             %
% 1. Downsampling (okDsample.m)                                           %
% 2. Make covaraince matrix (okMakeInSARCovMat.m)                         %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. Dataset: Character. The field which contains the displacement, and   %
%    the local coordinates                                                %
%                                                                         %
% Name-Value pairs:                                                       %
% 'deramp': Deramp or not before auto-correlation [0 or 1, default=1]     %
% 'func': Functions to fit the auto-correlation                           %
%    [default='exp1','poly3','poly5','poly7','poly9']                     %
%                                                                         %
% Example:                                                                %
% DataStruct = okInSARCov(DataStruct,'Subset', ...                        %
%       'deramp',1, ...                                                   %
%       'func',{'exp1','poly5'})                                          %
%                                                                         %
% Note that:                                                              %
% This routine uses the built-in function fit to perform the curve-fitting%
%                                                                         %
% Output:                                                                 %
% AutoCorr: Structure. Contain the auto-correlation results               %
%   AutoCorrelation: Auto-correlation coefficients                        %
%   RadialDistance: The radial distance from the center to other pixels   %
%   FittedFunction: The functions used to fit the pattern                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function AutoCorr = okInSARCov(DataStruct,Dataset,varargin)
warning('off','all')
% Parse function parameters
p = inputParser;
default_deramp = 1;
default_func = {'exp1','poly3','poly5','poly7','poly9'};
addParameter(p, 'deramp', default_deramp, @(x) isnumeric(x) && x == 1 || x == 0);
addParameter(p, 'func', default_func, @(x) iscell(x) || ischar(x));

parse(p, varargin{:});
deramp = p.Results.deramp;
func = p.Results.func;

if ~iscell(func)
    func = cellstr(func);
end


AutoCorr = DataStruct;
% Retrieve the data
Displ = DataStruct.(Dataset).Displ;
LocalX = DataStruct.(Dataset).LocalX;
LocalY = DataStruct.(Dataset).LocalY;


% Get the size of the displacement
N = numel(Displ);
[Row,Col] = size(Displ);
RowMat = repmat(transpose(fliplr(1:Row)),1,Col);
ColMat = repmat(1:Col,Row,1);

% Calculate average row and column pixel sizes
Rowtmp = zeros(Row-1,1);
Coltmp = zeros(Col-1,1);
for r = 1:length(Rowtmp)
    Rowtmp(r) = abs(LocalY(r+1,1) - LocalY(r,1)); 
end
for c = 1:length(Coltmp)
    Coltmp(c) = abs(LocalX(1,c) - LocalX(1,c+1));
end
PixelSizeR = mean(Rowtmp);
PixelSizeC = mean(Coltmp);

disp(strcat('Row pixel size:',32,num2str(PixelSizeR)))
disp(strcat('Col pixel size:',32,num2str(PixelSizeC)))
disp(' ')
%% Deramp if it is asked to do so
if deramp == 1
    disp('*** Deramping')
    rmNaN = find(~isnan(Displ));
    G = [RowMat(:),ColMat(:),ones(N,1)];
    GrmNaN = G(rmNaN,:);
    drmNaN = Displ(rmNaN);

    mramp = GrmNaN\drmNaN;
    Displ = Displ - reshape(G*mramp,Row,Col);
    
end


%% Find the autocorrelation
% Convert all NaN to 0
disp('*** Calculating auto-correlation')
Displ(isnan(Displ)) = 0;

% FFT
DisplFFT = fft2(Displ);

% Get the power spectrum of it
PowerSpec = DisplFFT.*conj(DisplFFT);


% Inverse FFT to get the autocorrelation
% Based on Wiener-Khinchin theorem
C = ifft2(PowerSpec);


% Normalize to the maximum is 1
Cnorm = C./N;
CenterVal = max(Cnorm(:));
Cnorm = Cnorm./CenterVal;

% Shift the zero-frequency to the middle
AutoCorrelation = fftshift(Cnorm);


%% Calculate the radial distance
disp('*** Calculating radian distance between the center and other pixels')
disp(' ')
[CenterR,CenterC] = ind2sub([Row,Col],find(AutoCorrelation==1));
Dist = sqrt(((RowMat-CenterR).*PixelSizeR).^2 + ((ColMat - CenterC).*PixelSizeC).^2);

%% Function fitting 
x = Dist(:);
y = double(AutoCorrelation(:));
disp('*** Fitting functions to radial distance v.s. correlation coefficient')
disp('*** This would take a little while...')
for i = 1:length(func)
    disp(strcat('Function:',32,func{i}))

    % Fit function
    % Take extra care when function is exp1
    if strcmp(func{i},'exp1')
        lb = [0, -1e-2];
        ub = [1, -1e-6];
        m = fit(x,y,func{i},'Lower',lb,'Upper',ub);
    else
        m = fit(x,y,func{i});
    end
    
    % Output
    AutoCorr.(Dataset).(func{i}) = m;
end


%% Output 
AutoCorr.(Dataset).AutoCorrelation = AutoCorrelation;
AutoCorr.(Dataset).RadialDistance = Dist;
AutoCorr.(Dataset).FittedFunction = func;


end
