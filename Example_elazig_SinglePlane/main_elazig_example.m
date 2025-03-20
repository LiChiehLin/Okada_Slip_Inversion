clear
%% okLoadData.m
datadir = 'elazig';
Displname = strcat(datadir,'/','UnwrapPhase.tif');
Aziname = strcat(datadir,'/','Azimuth.tif');
Incname = strcat(datadir,'/','Incidence.tif');
Cohname = strcat(datadir,'/','Coherence.tif');

DataStruct = okLoadData('Displacement',Displname, ...
    'Azimuth',Aziname, ...
    'Incidence',Incname, ...
    'Coherence',Cohname, ...
    'Processor','ISCE', ...
    'LookDir','r');
% Convert radian to meters
wavel = 0.055;
DataStruct.Original.Displ = DataStruct.Original.Displ.*wavel./4./pi;

%% Mask data (okMaskData.m)
Threshold = 0.35;
DataStruct = okMaskData(DataStruct,'Original','Displ',DataStruct.Original.Coherence,'Threshold',Threshold);

%% Crop the data (okMakeDataSubset.m)
rmin = 3500;
rmax = 5932;
cmin = 1;
cmax = 4926;
DataStruct = okMakeDataSubset(DataStruct,'Original',rmin,rmax,cmin,cmax);

%% Convert to local coordinates (okLLtoLocalXY.m)
lon_origin = 38;
DataStruct = okLLtoLocalXY(DataStruct,'Subset',lon_origin);

%% Downsample (okDsample.m);
DsampleParam = okMakeDsampleParam('quadtree', ...
    'Criterion','curvature', ...
    'Tolerance',0.00007, ...
    'LeafMinSize',1000, ...
    'LeafMaxSize',200000, ...
    'NaNPixelAllow',0.5);
DataStruct = okDsample(DataStruct,'Subset',DsampleParam);

%% Make fault model (okMakeFaultModel.m)
StartXYZ = [359099, 4253530, 0];
Length = 60000;
Width = 20000;
Strike = 242;
Dip = 80;
PatchStrike = 24;
PatchDip = 8;
FaultModel = okMakeFaultModel(StartXYZ,Length,Width,Strike,Dip,PatchStrike,PatchDip);

%% Make Green's function (okMakeGreenFunc.m)
RakeDS = -90;
RakeSS = 0;
SlipDS = 1;
SlipSS = 1;
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeDS,SlipDS,0,'GreenDS');
FaultModel = okMakeGreenFunc(DataStruct,'Dsample','LOS',FaultModel,RakeSS,SlipSS,0,'GreenSS');

%% Make Smoothing matrix (okMakeSmoothMat.m)
FaultModel = okMakeSmoothMat(FaultModel);

%% Invert fault slip (okInvertSlip.m)
SmoothParam = [0.00001,0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10];
SlipModel = okInvertSlip(DataStruct,'Dsample',FaultModel,{'GreenDS','GreenSS'},[11,12],'SmoothMat', ...
    'solver','nnlsq', ...
    'smoothsearch',SmoothParam);

%okPlot(SlipModel,'Inversion','slip1',1,'n',5,'title','Dip-slip')
%okPlot(SlipModel,'Inversion','slip2',1,'n',5,'title','Strike-slip')
okPlot(SlipModel,'Inversion','totalslip',1,'n',5)
%okPlot(SlipModel,'Inversion','rake',1,'n',5)
%okPlot(SlipModel,'Inversion','residual',1,'n',5,'MarkerSize',30)
%okPlot(SlipModel,'Inversion','lcurve',1,'n',5,'title','L-curve')

%% Output to GMT (okOutputGMT.m)
okOutputGMT(SlipModel,'slip', ...
    'n',5);
okOutputGMT(SlipModel,'lcurve', ...
    'n',5);
okOutputGMT(SlipModel,'green', ...
    'n',5);
okOutputGMT(SlipModel,'displ', ...
    'n',5);

