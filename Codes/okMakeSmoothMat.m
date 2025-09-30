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
% (Update: 2025.09.26)                                                    %
%  'dist-based': if cannot find patches within the distance, then smooth  %
%  it with its connecting patches                                         %
% (Update: 2025.09.14)                                                    %
%   Allows three options:                                                 %
%     1. 'equidist' (default):                                            %
%        Simple 2D Laplacian smoothing. Only works with the fault model   %
%        having the same patch sizes                                      %
%     2. 'dist-weighted':                                                 %
%        Finds the connected patches and weight them with the inverse of  %
%        their centroid distances                                         %
%     3. 'dist-based':                                                    %
%        Smooths the patches within the input distance threshold. Also    %
%        weights them with the centroid distances                         %
%                                                                         %
% Make smoothing matrix from the fault model generated from               %
% okMakeFaultModel.m leveraging okada85.m.                                %
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
function SmoothModel = okMakeSmoothMat(FaultModel,varargin)
p = inputParser;
% Parse options for the smoothing matrix
default_method = 'equidist';
default_dist = [];
default_tol = 1e-3;
addParameter(p,'method',default_method, @(x) ischar(x) && (strcmp(x,'equidist') || strcmp(x,'dist-weighted') || strcmp(x,'dist-based')));
addParameter(p,'dist',default_dist, @(x) isnumeric(x));
addParameter(p,'tol',default_tol, @(x) isnumeric(x));


parse(p, varargin{:});
method = p.Results.method;
dist = p.Results.dist;
tol = p.Results.tol;

% Check if "dist-based" option is put without the distance parameter
if strcmp(method,'dist-based') && isempty(dist)
    error('"dist-based" option requires option "dist"')
end

% Retrieve data
okFault = FaultModel.okFault;
PatchX = FaultModel.PatchX;
PatchY = FaultModel.PatchY;
PatchZ = FaultModel.PatchZ;
ASpatch = FaultModel.PatchCountStrike;
ADpatch = FaultModel.PatchCountDip;
ASpatchsize = FaultModel.AlongStrikePatchSize;
ADpatchsize = FaultModel.AlongDipPatchSize;

% If the fault model is combined with different patch sizes, then
% "equidist" will not work
if any(ASpatchsize ~= ADpatchsize) && strcmp(method,'equidist')
    error(['Combined fault model with different patch sizes does not allow "equidist". '  ...
        'try "dist-based" or "dist-weighted"'])
end

% Smoothing function
smFunc = method;

disp(' ')
disp("******* Constructing smoothing matrix okMakeSmoothMat.m *******")
disp(strcat('*** Smoothing function:',32,smFunc))
disp(strcat('*** Tolerance:',32,num2str(tol)))


Dim = okFault(end,8);
S = zeros(Dim);
SmoothModel = FaultModel;
if strcmp(method,'equidist')
    for i = 1:Dim
        % From Top-left patch to bottom-right patch
        ui = i-sum(ASpatch);
        di = i+sum(ASpatch);
        li = i-1;
        ri = i+1;
    
        S(i,i) = -4;
        % Detect edge cases, corners first
        if i == 1
        % 1. top-left corner
            S(i,ri) = 2;
            S(i,di) = 2;
        elseif i == sum(ASpatch)
        % 2. top-right corner
            S(i,li) = 2;
            S(i,di) = 2;
        elseif i == sum(ASpatch.*ADpatch) - sum(ASpatch)
        % 3. bottom-left corner
            S(i,ui) = 2;
            S(i,ri) = 2;
        elseif i == sum(ASpatch.*ADpatch)
        % 4. bottom-right corner
            S(i,ui) = 2;
            S(i,li) = 2;
        elseif i < sum(ASpatch)
        % 5. top row
            S(i,li) = 1;
            S(i,ri) = 1;
            S(i,di) = 2;
        elseif i > sum(ASpatch.*ADpatch) - sum(ASpatch)
        % 6. bottom row
            S(i,li) = 1;
            S(i,ri) = 1;
            S(i,ui) = 2;
        elseif mod(i,sum(ASpatch)) == 1
        % 7. left column
            S(i,ui) = 1;
            S(i,di) = 1;
            S(i,ri) = 2;
        elseif mod(i,sum(ASpatch)) == 0
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

elseif strcmp(method,'dist-weighted')
    % Find the connecting patches
    Neighbors = cell(size(PatchX,2),1);
    for i = 1:size(PatchX,2)
        % For each patch
        ntmp = i;
        for j = 1:4
            % For each side of the patch
            AInd = j;
            if j == 4
                BInd = 1;
            else
                BInd = j + 1;
            end
            % These are the coordinates of the one looking for neighbors
            Ax = PatchX(AInd,i); Bx = PatchX(BInd,i);
            Ay = PatchY(AInd,i); By = PatchY(BInd,i);
            Az = PatchZ(AInd,i); Bz = PatchZ(BInd,i);
            AB = [Bx-Ax,By-Ay,Bz-Az];
    
            for k = 1:4
                % For each side of all other patches
                CInd = k;
                if k == 4
                    DInd = 1;
                else
                    DInd = k + 1;
                end
                % These are the coordinates of the searched sides
                Cx = PatchX(CInd,:); Dx = PatchX(DInd,:);
                Cy = PatchY(CInd,:); Dy = PatchY(DInd,:);
                Cz = PatchZ(CInd,:); Dz = PatchZ(DInd,:);
    
                CA = [Cx-Ax;Cy-Ay;Cz-Az]';
                DA = [Dx-Ax;Dy-Ay;Dz-Az]';
    
                % 1. Cross product to see if they are colinear
                ABcrossCA = vecnorm(cross(repmat(AB,[size(PatchX,2),1]),CA),2,2);
                ABcrossDA = vecnorm(cross(repmat(AB,[size(PatchX,2),1]),DA),2,2);
    
                % 2. Project onto the line
                unitvec = AB/norm(AB);
        
                % 3. Set the bounds
                Abound = [0,dot(AB,unitvec)];
                Bbound = [dot(CA,repmat(unitvec,[size(PatchX,2),1]),2), dot(DA,repmat(unitvec,[size(PatchX,2),1]),2)];
    
                % 4. Determine if they overlapped
                tmp = zeros(size(Bbound,1),1);
                for l = 1:size(tmp,1)
                    if max([min(Abound),min(Bbound(l,:))]) + tol < min([max(Abound),max(Bbound(l,:))])
                        tmp(l) = true;
                    else
                        tmp(l) = false;
                    end
                end
                tmp = logical(tmp);
    
                % Determine if they are neighbors
                Neighborlogic = (ABcrossCA < tol) & (ABcrossDA < tol) & tmp;
                Neighborlogic(i) = false; % Maskout the current patch
                Ind = find(Neighborlogic);
            
                % Store the patch number
                ntmp = [ntmp,Ind'];
            end
        end
        ntmp = [ntmp(1),sort(unique(ntmp(2:end)),'ascend')]; %  clean up and store

        % Calculate the distance of the current patch to its connected patches
        N = length(ntmp);
        NInd = ntmp(2:end);
        Ccoord = okFault(i,1:3);
        Ndist = sqrt((okFault(:,1) - Ccoord(1)).^2 + (okFault(:,2) - Ccoord(2)).^2 + (okFault(:,3) - Ccoord(3)).^2);
        
        % Invert the distances so that larger distance will have lower
        % weights
        Mindist = min(Ndist(NInd));
        Maxdist = max(Ndist(NInd));
        Ndist = Maxdist + Mindist - Ndist;
        Normfactor = sum(Ndist(NInd));

        % Calculate the Laplacian numbers
        Laplanum = (Ndist./Normfactor)*(N-1);

        % Put the smoothing numbers in the cooresponding elements
        S(i,i) = -(N-1);
        S(i,NInd) = Laplanum(NInd);
        
        % Also output the connected patches to the fault model
        Neighbors{i} = ntmp;
    end
    SmoothModel.Neighbors = Neighbors;

elseif strcmp(method,'dist-based')
    Neighbors = cell(size(PatchX,2),1);
    % Construct the distance matrix to find patches to be smoothed together
    Distmat = zeros(size(okFault,1),size(okFault,1));
    for i = 1:size(okFault,1)
        Ccoord = okFault(i,1:3);
        Distmat(:,i) = sqrt((Ccoord(1) - okFault(:,1)).^2 + (Ccoord(2) - okFault(:,2)).^2 + (Ccoord(3) - okFault(:,3)).^2);
    end

    % Find the patches to be smoothed
    for i = 1:size(Distmat,1)
        ntmp = find(Distmat(i,:) <= dist);
        ntmp(ntmp == i) = [];
        ntmp = [i,ntmp]; % This is just to be consistent with 'dist-weighted'
        N = length(ntmp);
        NInd = ntmp(2:end);

        if ~isempty(NInd)
            % Invert the distance
            Mindist = min(Distmat(i,NInd));
            Maxdist = max(Distmat(i,NInd));
            Ndist = Distmat(i,:);
            Ndist = Maxdist + Mindist - Ndist;
            Normfactor = sum(Ndist(NInd));
        else
            % If no patches were found, then smooth with the one connected
            % to it (Copy from the above)
            for j = 1:4
                AInd = j;
                if j == 4
                    BInd = 1;
                else
                    BInd = j + 1;
                end
                % These are the coordinates of the one looking for neighbors
                Ax = PatchX(AInd,i); Bx = PatchX(BInd,i);
                Ay = PatchY(AInd,i); By = PatchY(BInd,i);
                Az = PatchZ(AInd,i); Bz = PatchZ(BInd,i);
                AB = [Bx-Ax,By-Ay,Bz-Az];
        
                for k = 1:4
                    % For each side of all other patches
                    CInd = k;
                    if k == 4
                        DInd = 1;
                    else
                        DInd = k + 1;
                    end
                    % These are the coordinates of the searched sides
                    Cx = PatchX(CInd,:); Dx = PatchX(DInd,:);
                    Cy = PatchY(CInd,:); Dy = PatchY(DInd,:);
                    Cz = PatchZ(CInd,:); Dz = PatchZ(DInd,:);
        
                    CA = [Cx-Ax;Cy-Ay;Cz-Az]';
                    DA = [Dx-Ax;Dy-Ay;Dz-Az]';
        
                    % 1. Cross product to see if they are colinear
                    % Might need to put a tolerance here to avoid float number
                    % errors
                    ABcrossCA = vecnorm(cross(repmat(AB,[size(PatchX,2),1]),CA),2,2);
                    ABcrossDA = vecnorm(cross(repmat(AB,[size(PatchX,2),1]),DA),2,2);
        
                    % 2. Project onto the line
                    unitvec = AB/norm(AB);
            
                    % 3. Set the bounds
                    Abound = [0,dot(AB,unitvec)];
                    Bbound = [dot(CA,repmat(unitvec,[size(PatchX,2),1]),2), dot(DA,repmat(unitvec,[size(PatchX,2),1]),2)];
        
                    % 4. Determine if they overlapped
                    tmp = zeros(size(Bbound,1),1);
                    for l = 1:size(tmp,1)
                        if max([min(Abound),min(Bbound(l,:))]) + tol < min([max(Abound),max(Bbound(l,:))])
                            tmp(l) = true;
                        else
                            tmp(l) = false;
                        end
                    end
                    tmp = logical(tmp);
        
                    % Determine if they are neighbors
                    Neighborlogic = (ABcrossCA < tol) & (ABcrossDA < tol) & tmp;
                    Neighborlogic(i) = false; % Maskout the current patch
                    Ind = find(Neighborlogic);
                
                    % Store the patch number
                    ntmp = [ntmp,Ind'];
                end
            end
            % Invert the distance
            N = length(ntmp);
            NInd = ntmp(2:end);
            Mindist = min(Distmat(i,NInd));
            Maxdist = max(Distmat(i,NInd));
            Ndist = Distmat(i,:);
            Ndist = Maxdist + Mindist - Ndist;
            Normfactor = sum(Ndist(NInd));
        end

        % Calculate the Laplacian numbers
        Laplanum = (Ndist./Normfactor)*(N-1);

        % Put the smoothing numbers in the cooresponding elements
        S(i,i) = -(N-1);
        S(i,NInd) = Laplanum(NInd);

        % Also output the connected patches to the fault model
        Neighbors{i} = ntmp;
    end
    SmoothModel.Neighbors = Neighbors;

end

% Output
SmoothModel.SmoothMat = S;

end



