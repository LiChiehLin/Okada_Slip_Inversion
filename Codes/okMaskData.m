%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.19                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okMaskData.m               ***             %
%             ***********************************************             %
%                                                                         %
% Mask displacement given the mask criterion                              %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Data structure made by okLoadData.m           %
% 2. TobeMasked: Character. The matrix to be masked                       %
% 3. MaskCriterion: Matrix. The criterion used for masking                %
% 4. Threshold: Numeric. The threshold for masking (optional)             %
%                                                                         %
% If threshold is input, then it will treat as a threshold for            %
% MaskCriterion. If threshold is left un-input, then it assumes the       %
% MaskCriterion is already a mask                                         %
%                                                                         %
% Example:                                                                %
% MaskStruct = okMaskData(DataStruct,'Original','Displ', ...              %
%       DataStruct.Original.Coherence,'Threshold',0.35)                   %
% MaskStruct = okMaskData(DataStruct,'Original','Displ', ...              %
%       Mask)                                                             %
%                                                                         %
% Output:                                                                 %
% 1. MaskStruct: Structure. Contains the masked stuff                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MaskStruct = okMaskData(DataStruct,Dataset,TobeMasked,MaskCriterion,varargin)
% Parse the parameter
p = inputParser;
default_threshold = 'auto';
addParameter(p, 'Threshold', default_threshold, @(x) isnumeric(x) && length(x)==1);

parse(p, varargin{:});
Threshold = p.Results.Threshold;

disp(' ')
disp(strcat('*** Threshold:',32,num2str(Threshold)))

% Prepare the mask
OutFieldname = strcat(TobeMasked,'Orig');
if ischar(Threshold)
    Mask = MaskCriterion;
elseif ~ischar(Threshold)
    Mask = MaskCriterion;
    Mask(Mask < Threshold) = nan;
end

% Start masking
disp(strcat('*** Masking dataset',TobeMasked))
Data = DataStruct.(Dataset).(TobeMasked);
DataMasked = Data.*Mask;

% Output
MaskStruct = DataStruct;
MaskStruct.(Dataset).(OutFieldname) = Data;
MaskStruct.(Dataset).(TobeMasked) = DataMasked;
MaskStruct.(Dataset).MaskThreshold = Threshold;
end