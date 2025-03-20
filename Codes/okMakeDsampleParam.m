%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.07                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***           okMakeDsampleParam.m          ***             %
%             ***********************************************             %
%                                                                         %
% Make the parameter file for okDsample.m                                 %
% Leave blank to use the hard-coded default values, not recommended       %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. FuncType: Downsampling function type. Each requires different input  %
%    parameters. See below the supprted function and parameters           %
%    1.1 Quadtree:                                                        %
%        'Criterion': ['variance' or 'curvature'; default: variance]      %
%        See okQuadtreeD.m for more detail about the criterion            %
%        'Tolerance': [x > 0; default: 0.05]                              %
%        'LeafMinSize': [x > 0; default: 100]                             %
%        'LeafMaxSize': [x > 0; default: 4000]                            %
%        'NaNPixelAllow': [x > 0 && x < 1; default: 0.5]                  %
%    1.2 Uniform:                                                         %
%    1.3 R-based:                                                         %
%                                                                         %
% Output:                                                                 %
% 1. FuncParam: Structure. Parameters ready for okDsample.m               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FuncParam = okMakeDsampleParam(FuncType,varargin)
p = inputParser;
% Make different parameter sets based on the function type
if strcmp(FuncType,'quadtree') || strcmp(FuncType,'Quadtree')
    default_criterion = 'variance';
    default_tolerance = 0.05;
    default_minsize = 100;
    default_maxsize = 4000;
    default_nanallow = 0.5;
    addParameter(p, 'Criterion', default_criterion, @(x) validCriterion(x));
    addParameter(p, 'Tolerance', default_tolerance, @(x) validTolerance(x));
    addParameter(p, 'LeafMinSize', default_minsize, @(x) validLeafSize(x));
    addParameter(p, 'LeafMaxSize', default_maxsize, @(x) validLeafSize(x));
    addParameter(p, 'NaNPixelAllow', default_nanallow, @(x) validNaNAllow(x));

    % parse the input to p
    parse(p, varargin{:});

    % Store all the parameters
    FuncParam.FunctionType = FuncType;
    FuncParam.Criterion = p.Results.Criterion;
    FuncParam.Tolerance = p.Results.Tolerance;
    FuncParam.LeafMinSize = p.Results.LeafMinSize;
    FuncParam.LeafMaxSize = p.Results.LeafMaxSize;
    FuncParam.NaNPixelAllow = p.Results.NaNPixelAllow;
    disp(FuncParam);


elseif strcmp(FuncType,'uniform') || strcmp(FuncType,'Uniform')
    % To-be updated

elseif strcmp(FuncType,'rosampling') || strcmp(FuncType,'Rosampling') || strcmp(FuncType,'R-based')
    % To-be updated

else
    error_message = ['Cannot match the supporting function type. Now supports: \n', ...
        ' 1. Quadtree downsampling (put [quadtree] or [Quadtree]) (Jonsson 2002 or Simons et al., 2002) \n', ...
        ' 2. Uniform downsampling (put [uniform] or [Uniform])\n', ...
        ' 3. R-based downsampling (put [rosampling] or [Rosampling] or [R-based]) (Lohman and Simons, 2005)\n'];
    error('u:stuffed:it',error_message)

end

end

%% Code blocks for checking input argumments (case sensitive)
% Quadtree
function check = validCriterion(x)
if strcmp(x,'variance') || strcmp(x,'curvature')
    check = true;
else
    error_message = ['Quadtree criterion: either [variance] or [curvature] \n', ...
        '  variance: (Jonsson, 2002) \n', ...
        '  curvature: (Simons et al., 2002) \n'];
    error('u:stuffed:it',error_message)
end
end

function check = validTolerance(x)
if isnumeric(x) && isscalar(x) && x > 0
    check = true;
else
    error_message = 'Tolerance: \n Has to be a number larger than 0';
    error('u:stuffed:it',error_message)
end
end

function check = validLeafSize(x)
if isnumeric(x) && isscalar(x) && x > 0
    check = true;
else
    error_message = 'Leaf size: Has to be a number larger than 0';
    error('u:stuffed:it',error_message)
end
end

function check = validNaNAllow(x)
if isnumeric(x) && isscalar(x) && x > 0 && x < 1
    check = true;
else
    error_message = 'NaN pixel allowed: Has to be a number between 0 and 1';
    error('u:stuffed:it',error_message)
end
end