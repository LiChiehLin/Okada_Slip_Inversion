%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.10                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***                okPlot.m                 ***             %
%             ***********************************************             %
%                                                                         %
% (Update 2025.05.13)                                                     %
%   Add 1 figure type:                                                    %
%   1. GNSS: Plot the GNSS measurements                                   %
%   Supports different 'Inversion' figure types. Now supports InSAR and   %
%   GNSS model fit plots                                                  %
% (Update 2025.03.26)                                                     %
%   Add 2 figure types:                                                   %
%   1. AutoCorr: Plot the auto-correlation result                         %
%   2. Covaraince: Plot the covaraince matrix of InSAR measurements       %
%                                                                         %
% A holistic plotting routine for make figures according to the FigType   %
% and the input data structure                                            %
% See below the supported figure types:                                   %
%                                                                         %
% 1. FigType='Displ'                                                      %
%    To visualize the read-in dispalcement to determine where the fault is%
%    and to make subsets according to the coordinates                     %
% 2. FigType='GNSS'                                                       %
%    To visualize the read-in GNSS measurements                           %
% 3. FigType='Azimuth'                                                    %
%    To visualize the Azimuth angle                                       %
% 4. FigType='Incidence'                                                  %
%    To visualize the Incidence angle                                     %
% 5. FigType='FaultModel'                                                 %
%    To visualize the fault model made from okMakeFaultModel.m            %
% 6. FigType='Dsample'                                                    %
%    To visualzie the downsampled result. Both data points and tree leaves%
% 7. FigType='GreenFunc'                                                  %
%    To visualize the forward okada Green's function                      %
% 8. FigType='Smooth'                                                     %
%    To visualize the smoothing matrix with matlab spy command            %
% 9. FigType='Inversion'                                                  %
%    To visualize the inversion results                                   %
% 10. FigType='AutoCorr'                                                  %
%    To visualize the InSAR auto-correlation result                       %
% 11. FigType='Covariance'                                                %
%    To visualize the covariance matrix                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% See below for the sub figure options and examples                       %
% FigType = 'Displ'                                                       %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
% okPlot(DataStruct.Subset,'Displ','clim',[-2 2])                         %
%                                                                         %
% FigType = 'GNSS'                                                        %
% 'scale': Multiplier for visualization                                   %
% okPlot(DataStruct.Subset,'GNSS','scale',1000)                           %
%                                                                         %
% FigType = 'Azimuth'                                                     %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
% okPlot(DataStruct.Subset,'Azimuth')                                     %
%                                                                         %
% FigType = 'Incidence'                                                   %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
% okPlot(DataStruct.Subset,'Incidence')                                   %
%                                                                         %
% FigType = 'FaultModel'                                                  %
% 'displ': Plot displacement field on top of the fault model              %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
% 'GNSS': Plot GNSS measurement on top of the fault model                 %
% 'scale': Multiplier for visualization                                   %
% okPlot(FaultModel,'FaultModel','displ','DataStruct.Subset', ...         %
%   'clim',[-2 2])                                                        %
% okPlot(FaultModel,'FaultModel','GNSS','DataStruct.Original', ...        %
%   'scale',1000)                                                         %
%                                                                         %
% FigType = 'Dsample'                                                     %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'MarkerSize': Point size                                                %
% 'title': Figure title                                                   %
% okPlot(DataStruct,'Dsample','clim',[-2 2],'MarkerSize',20)              %
%                                                                         %
% FigType = 'GreenFunc'                                                   %
% 'displ': The data structure used to make the Green's function           %
%    This is for the sake of the coordinates of the Green's function. So  %
%    input the one that you used to make the Green's function             %
% 'dset': Green's function field name                                     %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'MarkerSize': Point size                                                %
% 'title': Figure title                                                   %
% okPlot(FaultModel,'GreenFunc','dset','GreenAscDS', ...                  %
%   'displ',DataStruct.Dsample,'MarkerSize',20)                           %
%                                                                         %
% FigType = 'Smooth'                                                      %
% 'title': Figure title                                                   %
% okPlot(FaultModel,'Smooth')                                             %
%                                                                         %
% FigType = 'Inversion'                                                   %
% 'residual': Plot the residual                                           %
%   This has to be used with 'n' to specify which result to plot and 'n'  %
%   can only be 1 number                                                  %
% 'dset': Plot which dataset of the model prediction                      %
% 'lcurve': Plot the L-curve                                              %
% 'slip1': Plot the slip distribution (Normally is the Dip-slip)          %
% 'slip2': Plot the slip distribution (Normally is the Strike-slip)       %
% 'totalslip': Plot the total slip                                        %
% 'rake': Plot the rake angle                                             %
% 'n': Plot the specified smoothing parameter results                     %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
% okPlot(SlipModel,'Inversion','lcurve',1)                                %
% okPlot(SlipModel,'Inversion','lcurve',1,'n',95)                         %
% okPlot(SlipModel,'Inversion','residual',1,'n',95)                       %
% okPlot(SlipModel,'Inversion','residual',1,'n',95,'dset','GNSS')         %
% okPlot(SlipModel,'Inversion','residual',1,'n',95,'dset','InSAR')        %
% okPlot(SlipModel,'Inversion','slip1',1)                                 %
% okPlot(SlipModel,'Inversion','slip1',1,'slip2',1,'totalslip',1)         %
% okPlot(SlipModel,'Inversion','totalslip',1,'rake',1,'n',95)             %
%                                                                         %
% FigType = 'AutoCorr'                                                    %
% 'func': Plot the fitted function                                        %
% 'one': Plot all in 1 figure                                             %
% 'title': Figure title                                                   %
%                                                                         %
% FigType = 'Covariance'                                                  %
% 'clim': Colorbar limits                                                 %
% 'cmap': Colormap                                                        %
% 'title': Figure title                                                   %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. DataStruct: Structure. Containing the displ. matrix and its          %
%    attributes                                                           %
% 2. FigType: Character. Specify the figure type                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function okPlot(DataStruct,FigType,varargin)
p = inputParser;
default_title = [];
default_clim = 0;
default_cmap = 'default';
default_displ = 0;
default_markersize = 5;
default_dset = 0;
default_n = 0;
default_residual = 0;
default_lcurve = 0;
default_slip1 = 0;
default_slip2 = 0;
default_totalslip = 0;
default_rake = 0;
default_func = [];
default_one = 0;
default_scale = 1;
default_linewidth = 0.5;

addParameter(p, 'title', default_title, @(x) ischar(x));
addParameter(p, 'clim', default_clim, @(x) isnumeric(x) && length(x)==2 && x(1)<x(2));
addParameter(p, 'cmap', default_cmap, @(x) ischar(x));
addParameter(p, 'displ', default_displ, @(x) isstruct(x) || ismatrix(x));
addParameter(p, 'MarkerSize', default_markersize, @(x) isnumeric(x));
addParameter(p, 'dset', default_dset, @(x) ischar(x));
addParameter(p, 'n', default_n, @(x) isnumeric(x));
addParameter(p, 'residual', default_residual, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'lcurve', default_lcurve, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'slip1', default_slip1, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'slip2', default_slip2, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'totalslip', default_totalslip, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'rake', default_rake, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'func', default_func, @(x) ischar(x) || iscell(x));
addParameter(p, 'one', default_one, @(x) isnumeric(x) && x==0 || x==1);
addParameter(p, 'scale', default_scale, @(x) isnumeric(x));
addParameter(p, 'GNSS', default_displ, @(x) isstruct(x));
addParameter(p, 'LineWidth', default_linewidth, @(x) isnumeric(x) && x > 0);

parse(p, varargin{:});

% Plot displacement
if strcmp(FigType,'Displ') || strcmp(FigType,'DisplOrig')
    disp('*** Plotting displacement field')
    if strcmp(FigType,'Displ')
        Displ = DataStruct.Displ;
    else
        Displ = DataStruct.DisplOrig;
    end
    
    [Row,Col] = size(Displ);
    if any(strcmp(fieldnames(DataStruct),'LocalX')) && any(strcmp(fieldnames(DataStruct),'LocalY'))
        disp('* Plot with the local coordinate system')
        X = DataStruct.LocalX;
        Y = DataStruct.LocalY;
        LocalFlag = 1;
    else
        disp('* Plot with row and column indices')
        X = repmat(1:Col,Row,1);
        Y = repmat(flipud(transpose(1:Row)),1,Col);
        LocalFlag = 0;
    end

    % Parse plotting parameters
    c = p.Results.clim;
    if c == 0
        c = [min(Displ(:)),max(Displ(:))];
    end
    cmap = p.Results.cmap;
    t = p.Results.title;

    % Report some attributes
    disp(strcat('clim:',32,num2str(c)))
    disp(strcat('colormap:',32,cmap))

    % Start plotting
    figure()
    if LocalFlag == 1
        imagesc([X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Displ);
        clim(c);colormap(cmap)
        title(t)
        set(gca, 'YDir','normal');
        xlabel('Local X')
        ylabel('Local Y')
        axis equal
    elseif LocalFlag == 0
        imagesc([X(1,1),X(1,end)],[Y(end,1),Y(1,1)],Displ);
        clim(c);colormap(cmap)
        set(gca, 'YDir','reverse');
        title(t)
        xlabel('Column index')
        ylabel('Row index')
        axis equal
    end
    colorbar()

elseif strcmp(FigType,'GNSS')
    disp('*** Plotting GNSS displacement/velocity vectors')
    X = DataStruct.LocalX;
    Y = DataStruct.LocalY;
    EWdispl = DataStruct.DisplEW;
    NSdispl = DataStruct.DisplNS;
    UDdispl = DataStruct.DisplUD;
    scale = p.Results.scale;
    LW = p.Results.LineWidth;

    %  Plot horizontal and vertical displacements
    figure()
    subplot(1,2,1)
    quiver(X,Y,EWdispl*scale,NSdispl*scale,'off','LineWidth',LW);title('Horizontal motion')
    xlabel('Local X (m)');ylabel('Local Y (m)');axis equal
    subplot(1,2,2)
    quiver(X,Y,zeros(length(X),1),UDdispl*scale,'off','LineWidth',LW);title('Vertical motion')
    xlabel('Local X (m)');ylabel('Local Y (m)');axis equal
    


% Plot azimuth angle
elseif strcmp(FigType,'Azimuth')
    disp('*** Plotting Azimuth angle')
    Azi = DataStruct.Azimuth;
    [Row,Col] = size(Azi);
    if any(strcmp(fieldnames(DataStruct),'LocalX')) && any(strcmp(fieldnames(DataStruct),'LocalY'))
        disp('* Plot with the local coordinate system')
        X = DataStruct.LocalX;
        Y = DataStruct.LocalY;
        LocalFlag = 1;
    else
        disp('* Plot with row and column indices')
        X = repmat(1:Col,Row,1);
        Y = repmat(flipud(transpose(1:Row)),1,Col);
        LocalFlag = 0;
    end

    % Parse plotting parameters
    c = p.Results.clim;
    if c == 0
        c = [min(Azi(:)),max(Azi(:))];
    end
    cmap = p.Results.cmap;
    t = p.Results.title;

    % Report some attributes
    disp(strcat('clim:',32,num2str(c)))
    disp(strcat('colormap:',32,cmap))

    % Start plotting
    figure()
    if LocalFlag == 1
        imagesc([X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Azi);
        clim(c);colormap(cmap)
        title(t)
        set(gca, 'YDir','normal');
        xlabel('Local X')
        ylabel('Local Y')
        axis equal
    elseif LocalFlag == 0
        imagesc([X(1,1),X(1,end)],[Y(end,1),Y(1,1)],Azi);
        clim(c);colormap(cmap)
        title(t)
        set(gca, 'YDir','reverse');
        xlabel('Column index')
        ylabel('Row index')
        axis equal
    end
    colorbar()

% Plot Incidence angle
elseif strcmp(FigType,'Incidence')
    disp('*** Plotting Incidence angle')
    Inc = DataStruct.Incidence;
    [Row,Col] = size(Inc);
    if any(strcmp(fieldnames(DataStruct),'LocalX')) && any(strcmp(fieldnames(DataStruct),'LocalY'))
        disp('* Plot with the local coordinate system')
        X = DataStruct.LocalX;
        Y = DataStruct.LocalY;
        LocalFlag = 1;
    else
        disp('* Plot with row and column indices')
        X = repmat(1:Col,Row,1);
        Y = repmat(flipud(transpose(1:Row)),1,Col);
        LocalFlag = 0;
    end

    % Parse plotting parameters
    c = p.Results.clim;
    if c == 0
        c = [min(Inc(:)),max(Inc(:))];
    end
    cmap = p.Results.cmap;
    t = p.Results.title;

    % Report some attributes
    disp(strcat('clim:',32,num2str(c)))
    disp(strcat('colormap:',32,cmap))

    % Start plotting
    figure()
    if LocalFlag == 1
        imagesc([X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Inc);
        clim(c);colormap(cmap)
        title(t)
        set(gca, 'YDir','normal');
        xlabel('Local X')
        ylabel('Local Y')
        axis equal
    elseif LocalFlag == 0
        imagesc([X(1,1),X(1,end)],[Y(end,1),Y(1,1)],Inc);
        clim(c);colormap(cmap)
        title(t)
        set(gca, 'YDir','reverse');
        xlabel('Column index')
        ylabel('Row index')
        axis equal
    end
    colorbar()

elseif strcmp(FigType,'Coherence')
    disp('*** Plotting Coherence angle')
    Coh = DataStruct.Coherence;
    [Row,Col] = size(Coh);
    if any(strcmp(fieldnames(DataStruct),'LocalX')) && any(strcmp(fieldnames(DataStruct),'LocalY'))
        disp('* Plot with the local coordinate system')
        X = DataStruct.LocalX;
        Y = DataStruct.LocalY;
        LocalFlag = 1;
    else
        disp('* Plot with row and column indices')
        X = repmat(1:Col,Row,1);
        Y = repmat(flipud(transpose(1:Row)),1,Col);
        LocalFlag = 0;
    end

    % Parse plotting parameters
    c = p.Results.clim;
    if c == 0
        c = [min(Coh(:)),max(Coh(:))];
    end
    cmap = p.Results.cmap;
    t = p.Results.title;

    % Report some attributes
    disp(strcat('clim:',32,num2str(c)))
    disp(strcat('colormap:',32,cmap))

    % Start plotting
    figure();title(t)
    if LocalFlag == 1
        imagesc([X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Coh);
        clim(c);colormap(cmap)
        set(gca, 'YDir','normal');
        xlabel('Local X')
        ylabel('Local Y')
        axis equal
    elseif LocalFlag == 0
        imagesc([X(1,1),X(1,end)],[Y(end,1),Y(1,1)],Coh);
        clim(c);colormap(cmap)
        set(gca, 'YDir','reverse');
        xlabel('Column index')
        ylabel('Row index')
        axis equal
    end
    colorbar()

% Plot fault model
elseif strcmp(FigType,'FaultModel')
    okFault = DataStruct.okFault;
    PatchX = DataStruct.PatchX;
    PatchY = DataStruct.PatchY;
    PatchZ = DataStruct.PatchZ;

    % Parse plotting parameters
    DisplParse = p.Results.displ;
    c = p.Results.clim;
    cmap = p.Results.cmap;
    t = p.Results.title;
    GNSSParse = p.Results.GNSS;
    scale = p.Results.scale;

    % Determine whether plot InSAR field or GNSS
    if isstruct(DisplParse)
        Displ = DisplParse.Displ;
        X = DisplParse.LocalX;
        Y = DisplParse.LocalY;
        DisplFlag = 1;
    elseif ismatrix(DisplParse) && DisplParse ~= 0
        Displ = DisplParse;
        [Row,Col] = size(Displ);
        X = repmat(1:Col,Row,1);
        Y = repmat(flipud(transpose(1:Row)),1,Col);
        DisplFlag = 1;
    elseif DisplParse == 0
        DisplFlag = 0;
    end
    if isstruct(GNSSParse)
        X = GNSSParse.LocalX;
        Y = GNSSParse.LocalY;
        EWdispl = GNSSParse.DisplEW;
        NSdispl = GNSSParse.DisplNS; 
        UDdispl = GNSSParse.DisplUD;
        GNSSFlag = 1;
    elseif ~isstruct(GNSSParse) && GNSSParse == 0
        GNSSFlag = 0;
    end

    % Make color limits
    if c==0 && DisplFlag == 1
        c = [min(Displ(:)),max(Displ(:))];
    end


    % Report some attributes
    disp(strcat('clim:',32,num2str(c)))
    disp(strcat('colormap:',32,cmap))

    % Start plotting
    figure();hold on;axis equal;title(t)
    if DisplFlag == 1
        imagesc([X(1,1),X(1,end)],[Y(1,1),Y(end,1)],Displ);
        set(gca, 'YDir','normal');
        clim(c);colormap(cmap);
    end
    if GNSSFlag == 1
        quiver3(X,Y,zeros(length(X),1),EWdispl*scale,NSdispl*scale,zeros(length(X),1),'b')
        quiver3(X,Y,zeros(length(X),1),zeros(length(X),1),zeros(length(X),1),UDdispl*scale,'r')
    end
    fill3(PatchX,PatchY,PatchZ,'black','FaceAlpha',0.3,'EdgeColor','k')
    scatter3(okFault(:,1),okFault(:,2),okFault(:,3),'r.')
    xlabel('X');ylabel('Y');zlabel('Depth (m)')
    colorbar();
    hold off

% Plot downsampled data
elseif strcmp(FigType,'Dsample')
    % Retrieve data
    % Undownsampled 
    SourceDataset = DataStruct.Dsample.SourceDataset;
    Displ = DataStruct.(SourceDataset).Displ;
    [Row,Col] = size(Displ);
    X = 1:Col;
    Y = fliplr(1:Row);
    % Downsampled
    LeafCoord = DataStruct.Dsample.LeafCoord;
    DisplDsample = DataStruct.Dsample.Displ;
    XDsample = DataStruct.Dsample.LocalX;
    YDsample = DataStruct.Dsample.LocalY;

    % Parse parameters
    c = p.Results.clim;
    if c == 0
        c = [min(Displ(:)),max(Displ(:))];
    end
    cmap = p.Results.cmap;
    MarkerSize = p.Results.MarkerSize;
    t = p.Results.title;

    % Plot the tree leaves
    figure();hold on;title(t)
    imagesc([X(1),X(end)],[Y(end),Y(1)],Displ);
    set(gca, 'YDir','reverse');
    xlabel('X');ylabel('Y')
    clim(c);colormap(cmap);colorbar()
    fill(transpose([LeafCoord(:,3),LeafCoord(:,4),LeafCoord(:,4),LeafCoord(:,3)]),transpose([LeafCoord(:,1),LeafCoord(:,1),LeafCoord(:,2),LeafCoord(:,2)]),'black','FaceAlpha',0)
    hold off

    % Plot the downsampled points
    figure();title(t)
    scatter(XDsample,YDsample,MarkerSize,DisplDsample,'filled')
    xlabel('X');ylabel('Y')
    clim(c);colormap(cmap);colorbar()

elseif strcmp(FigType,'GreenFunc')
    DisplParse = p.Results.displ;
    dset = p.Results.dset;
    c = p.Results.clim;
    cmap = p.Results.cmap;
    MarkerSize = p.Results.MarkerSize;
    t = p.Results.title;
    

    % Check if the input field name exist
    if any(strcmp(fieldnames(DataStruct),dset))
        LocalX = DisplParse.LocalX;
        LocalY = DisplParse.LocalY;
        GreenFunctmp = DataStruct.(dset);
        if isstruct(GreenFunctmp)
            fields = fieldnames(GreenFunctmp);
            N = numel(fields);
            GreenFunc = cell(N,1);
            for i = 1:N
                GreenFunc{i} = sum(GreenFunctmp.(fields{i}),2);
            end
        else
            GreenFunc = sum(GreenFunctmp,2);
            N = 1;
        end


        

        % Start plotting
        if N == 1
            % This is basically the LOS Green's
            if c == 0
                c = [min(GreenFunc(:)),max(GreenFunc(:))];
            end
            % Just in case both clim are equal
            if c(1) == c(2)
                c(1) = c(1)-0.0001;
                c(2) = c(2)+0.0001;
            end

            figure()
            scatter(LocalX,LocalY,MarkerSize,GreenFunc,'filled')
            clim(c);colormap(cmap);colorbar();title(t)
            xlabel('Local X');ylabel('Local Y')
        else
            % This is 3D Green's
            % Left to right: EW, NS, UD
            Direc = {'EW','NS','UD'};
            figure();hold on
            for i = 1:N
                if c == 0
                    tmp = GreenFunc{i};
                    cplot = [min(tmp(:)),max(tmp(:))];
                else
                    cplot = c;
                end
                % Just in case both clim are equal
                if cplot(1) == cplot(2)
                    cplot(1) = cplot(1)-0.0001;
                    cplot(2) = cplot(2)+0.0001;
                end
                subplot(1,N,i)
                scatter(LocalX,LocalY,MarkerSize,GreenFunc{i},'filled')
                clim(cplot);colormap(cmap);colorbar();title(strcat(t,32,Direc{i}))
                xlabel('Local X');ylabel('Local Y')
            end

        end
    else
        erorr('Cannot find input Green function field name')
    end

elseif strcmp(FigType,'Smooth')
    S = DataStruct.SmoothMat;
    t = p.Results.title;

    figure();title(t)
    spy(S)
    xlabel('Patch index');ylabel('Patch index')

elseif strcmp(FigType,'Inversion')
    Lcurve = DataStruct.Lcurve;
    SmoothParam = DataStruct.SmoothParam;
    Displ = DataStruct.Displ;
    LocalX = DataStruct.LocalX;
    LocalY = DataStruct.LocalY;
    ModelPrediction = DataStruct.ModelPrediction;
    DataType = DataStruct.DataType;
    FaultSlip = DataStruct.FaultSlip;
    TotalSlip = DataStruct.TotalSlip;
    Rake = DataStruct.Rake;
    PatchX = DataStruct.FaultModel.PatchX;
    PatchY = DataStruct.FaultModel.PatchY;
    PatchZ = DataStruct.FaultModel.PatchZ;
    PatchCount = size(DataStruct.FaultModel.SmoothMat,1);
    N = length(SmoothParam);
    Inc = ceil(N/10);

    % Parse parameters
    c = p.Results.clim;
    cmap = p.Results.cmap;
    n = p.Results.n;
    residual = p.Results.residual;
    lcurve = p.Results.lcurve;
    slip1 = p.Results.slip1;
    slip2 = p.Results.slip2;
    totalslip = p.Results.totalslip;
    rake = p.Results.rake;
    MarkerSize = p.Results.MarkerSize;
    t = p.Results.title;
    dset = p.Results.dset;
    scale = p.Results.scale;
    LW = p.Results.LineWidth;
    

    % Plot the L-curve
    if lcurve == 1
        if n == 0
            figure()
            plot(Lcurve(:,1),Lcurve(:,2),'-o')
            xlabel('Solution residual ||Lm||^2 ');ylabel('Model residual ||Gm-d||^2')
            title(t)
        elseif n ~= 0 && length(n) == 1
            figure();hold on
            plot(Lcurve(:,1),Lcurve(:,2),'-o')
            plot(Lcurve(n,1),Lcurve(n,2),'r*')
            xlabel('Solution residual ||Lm||^2 ');ylabel('Model residual ||Gm-d||^2')
            title(t)
        end
    end

    if n ~= 0
        % Just make figures for the input numbers
        plotseq = n;
    else
        % Make subplots for a quick skim
        plotseq = 1:Inc:N;
    end
    
    % Plot slip1 (Dip-slip, usually)
    if slip1 == 1
        figure()
        for i = plotseq
            nexttile
            PatchColor = FaultSlip(1:PatchCount,i);
            fill3(PatchX,PatchY,PatchZ,PatchColor);
            if c == 0
                c2 = [min(PatchColor),max(PatchColor)];
                % Just in case both clims are equal
                if c2(1) == c2(2)
                    c2(1) = c2(1)-0.0001;
                    c2(2) = c2(2)+0.0001;
                end
                clim(c2)
            else
                clim(c)
            end
            

            colormap(cmap);colorbar()
            xlabel('X'),ylabel('y')
            title(t);
            subtitle(strcat('Smoothing:',32,num2str(SmoothParam(i)),32,'n:',32,num2str(i)))
        end
    end
    % Plot slip2 (Stirke-slip, usually)
    if slip2 == 1
        figure()
        for i = plotseq
            nexttile
            PatchColor = FaultSlip(PatchCount+1:PatchCount*2,i);
            fill3(PatchX,PatchY,PatchZ,PatchColor);
            if c == 0
                c2 = [min(PatchColor),max(PatchColor)];
                % Just in case both clims are equal
                if c2(1) == c2(2)
                    c2(1) = c2(1)-0.0001;
                    c2(2) = c2(2)+0.0001;
                end
                clim(c2)
            else
                clim(c)
            end

            colormap(cmap);colorbar()
            xlabel('X'),ylabel('y')
            title(t);
            subtitle(strcat('Smoothing:',32,num2str(SmoothParam(i)),32,'n:',32,num2str(i)))
        end
    end
    % Plot totalslip 
    if totalslip == 1
        figure()
        for i = plotseq
            nexttile
            PatchColor = TotalSlip(:,i);
            fill3(PatchX,PatchY,PatchZ,PatchColor);
            if c == 0
                c2 = [min(PatchColor),max(PatchColor)];
                % Just in case both clims are equal
                if c2(1) == c2(2)
                    c2(1) = c2(1)-0.0001;
                    c2(2) = c2(2)+0.0001;
                end
                clim(c2)
            else
                clim(c)
            end

            colormap(cmap);colorbar()
            xlabel('X'),ylabel('y')
            title(t);
            subtitle(strcat('Smoothing:',32,num2str(SmoothParam(i)),32,'n:',32,num2str(i)))
        end
    end
    % Plot rake
    if rake == 1
        figure()
        for i = plotseq
            nexttile
            PatchColor = Rake(:,i);
            fill3(PatchX,PatchY,PatchZ,PatchColor);
            if c == 0
                c2 = [min(PatchColor),max(PatchColor)];
                % Just in case both clims are equal
                if c2(1) == c2(2)
                    c2(1) = c2(1)-0.0001;
                    c2(2) = c2(2)+0.0001;
                end
                clim(c2)
            else
                clim(c)
            end

            colormap(cmap);colorbar()
            xlabel('X'),ylabel('y')
            title(t);
            subtitle(strcat('Smoothing:',32,num2str(SmoothParam(i)),32,'n:',32,num2str(i)))
        end
    end

    % Plot the residual
    if residual == 1 && n ~= 0 && length(n) == 1
        if dset == 0
            if c == 0
                c = [min(Displ(:)),max(Displ(:))];
            end
            rmNaN = find(~isnan(Displ(:)));
            Pred = ModelPrediction(:,n);
            Diff = Displ(rmNaN)-Pred(rmNaN);
            c2 =[min(Diff(:)),max(Diff(:))];
            RMSE = sqrt(sum(Diff.^2)/length(Displ(rmNaN)));
    
            disp('* Treat all displacement points as InSAR points')
            figure()
            subplot(1,3,1)
            scatter(LocalX,LocalY,MarkerSize,Displ,'filled')
            clim(c);colormap(cmap);colorbar()
            xlabel('X'),ylabel('Y'),title('Observation')
            subplot(1,3,2)
            scatter(LocalX,LocalY,MarkerSize,ModelPrediction(:,n),'filled')
            clim(c);colormap(cmap);colorbar()
            xlabel('X'),ylabel('Y'),title('Model prediction');subtitle(strcat('Smoothing:',32,num2str(SmoothParam(n))))
            subplot(1,3,3)
            scatter(LocalX,LocalY,MarkerSize,Displ-ModelPrediction(:,n),'filled')
            clim(c2);colormap(cmap)
            xlabel('X'),ylabel('Y'),title('Residual (Obs. - Pred.)');subtitle(strcat('RMSE:',32,num2str(RMSE)))
            colorbar()
        elseif strcmp(dset,'GNSS')
            ind = DataType{strcmp(DataType,'GNSS'),2};
            Xtmp = LocalX(ind);
            Comp1 = sum(Xtmp(1) == Xtmp);
            
            if Comp1 == 2
                % This implies only EW and NS were used in inversion
                Dlength = length(ind)/2;
                X = LocalX(ind(1):ind(Dlength));
                Y = LocalY(ind(1):ind(Dlength));
                EWdispl = Displ(ind(1):ind(Dlength));
                NSdispl = Displ(ind(Dlength+1):ind(Dlength*2));
                UDdispl = nan(Dlength,1);
                EWPred = ModelPrediction(ind(1):ind(Dlength),n);
                NSPred = ModelPrediction(ind(Dlength+1):ind(Dlength*2),n);
                UDPred = nan(Dlength,1);
            elseif Comp1 ~= 2
                % This implies all three EW, NS and UD were used in
                % inversion
                Dlength = length(ind)/3;
                X = LocalX(ind(1):ind(Dlength));
                Y = LocalY(ind(1):ind(Dlength));
                EWdispl = Displ(ind(1):ind(Dlength));
                NSdispl = Displ(ind(Dlength+1):ind(Dlength*2));
                UDdispl = Displ(ind(Dlength*2+1):ind(Dlength*3));
                EWPred = ModelPrediction(ind(1):ind(Dlength),n);
                NSPred = ModelPrediction(ind(Dlength+1):ind(Dlength*2),n);
                UDPred = ModelPrediction(ind(Dlength*2+1):ind(Dlength*3),n);
            end
            EWRMSE = sqrt(mean((EWdispl - EWPred).^2,'omitnan'));
            NSRMSE = sqrt(mean((NSdispl - NSPred).^2,'omitnan'));
            UDRMSE = sqrt(mean((UDdispl - UDPred).^2,'omitnan'));
            RMSE = mean([EWRMSE,NSRMSE,UDRMSE],'omitnan');
            
            disp('* Plot GNSS model fit')
            figure()
            subplot(1,2,1);axis equal;hold on
            quiver(X,Y,EWdispl*scale,NSdispl*scale,'off','b','LineWidth',LW)
            quiver(X,Y,EWPred*scale,NSPred*scale,'off','r','LineWidth',LW)
            legend({'Data','Modeled'})
            title(strcat(t,32,'Horizontal motion'))
            subtitle(strcat('RMSE:',32,num2str(RMSE)))
            subplot(1,2,2);axis equal;hold on
            quiver(X,Y,zeros(length(UDdispl),1),UDdispl*scale,'off','b','LineWidth',LW)
            quiver(X,Y,zeros(length(UDdispl),1),UDPred*scale,'off','r','LineWidth',LW)
            legend({'Data','Modeled'})
            title(strcat(t,32,'Vertical motion'))
            subtitle(strcat('RMSE:',32,num2str(RMSE)))

        elseif strcmp(dset,'InSAR')
            ind = [DataType{strcmp(DataType,'InSAR'),2}];

            X = LocalX(ind);
            Y = LocalY(ind);
            InSARdispl = Displ(ind);
            InSARpred = ModelPrediction(ind);
            if c == 0
                c = [min(InSARdispl(:)),max(InSARdispl(:))];
            end

            rmNaN = find(~isnan(InSARdispl(:)));
            Diff = InSARdispl(rmNaN)-InSARpred(rmNaN);
            c2 =[min(Diff(:)),max(Diff(:))];
            RMSE = sqrt(sum(Diff.^2)/length(InSARdispl(rmNaN)));

            disp('* Plot InSAR model fit')
            figure()
            subplot(1,3,1)
            scatter(X,Y,MarkerSize,InSARdispl,'filled')
            clim(c);colormap(cmap);colorbar()
            xlabel('X'),ylabel('Y'),title('Observation')
            subplot(1,3,2)
            scatter(X,Y,MarkerSize,InSARpred(:,n),'filled')
            clim(c);colormap(cmap);colorbar()
            xlabel('X'),ylabel('Y'),title('Model prediction');subtitle(strcat('Smoothing:',32,num2str(SmoothParam(n))))
            subplot(1,3,3)
            scatter(X,Y,MarkerSize,InSARdispl-InSARpred(:,n),'filled')
            clim(c2);colormap(cmap)
            xlabel('X'),ylabel('Y'),title('Residual (Obs. - Pred.)');subtitle(strcat('RMSE:',32,num2str(RMSE)))
            colorbar()
        end
    end

elseif strcmp(FigType,'AutoCorr')
    AutoCorrelation = DataStruct.AutoCorrelation;
    Dist = DataStruct.RadialDistance;
    x = Dist(:);
    y = AutoCorrelation(:);

    % Parse parameter
    func = p.Results.func;
    one = p.Results.one;
    t = p.Results.title;

    % Check title and func
    if isempty(t)
        t = 'Radial distance v.s. Auto-correlation';
    end
    if isempty(func)
        % User did not put any function input
        func = DataStruct.FittedFunction;
    elseif ~iscell(func)
        % User only put 1 function input
        func = cellstr(func);
    end



    disp('*** Plotting fitted function')
    if one == 1
        figure();hold on
        plot(x,y,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'DisplayName','Pixels')
        xlabel('Radial distance (m)');ylabel('Auto-correlation coefficient')
        for i = 1:length(func)
            disp(func{i})
            Function = DataStruct.(func{i});
            ypred = Function(x);

            plot(x,ypred,'DisplayName',func{i},'LineWidth',3)
        end
        legend;title(t)

    else
        for i = 1:length(func)
            disp(func{i})
            Function = DataStruct.(func{i});
            ypred = Function(x);

            figure();hold on
            plot(x,y,'.','MarkerFaceColor',[0.5 0.5 0.5],'MarkerEdgeColor',[0.5 0.5 0.5],'DisplayName','Pixels')
            plot(x,ypred,'DisplayName',func{i},'LineWidth',3)
            xlabel('Radial distance (m)');ylabel('Auto-correlation coefficient')
            legend;title(strcat(t,32,func{i}))
        end
    end
elseif strcmp(FigType,'Covariance')
    Covariance = DataStruct.Covariance;

    % Parse parameter
    c = p.Results.clim;
    cmap = p.Results.cmap;
    t = p.Results.title;

    if c == 0
        c = [min(Covariance(:)), max(Covariance(:))];

    end
    if isempty(t)
        t = 'Covariance matrix';
    end
    imagesc(Covariance);title(t);clim(c);colorbar();colormap(cmap)
    xlabel('Pixel index');ylabel('Pixel index')
            
end


end
