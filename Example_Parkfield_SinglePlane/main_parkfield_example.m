clear 
%% Load data (okLoadData.m)
dataname = 'Houlie_et_al_2014.txt';
table = readtable(dataname);
DataStruct = okLoadData('Longitude',table2array(table(:,'Lon_deg'))-360, ...
    'Latitude',table2array(table(:,'Lat_deg')), ...
    'DisplacementEW',table2array(table(:,'De_cm')).*0.01, ...
    'DisplacementNS',table2array(table(:,'Dn_cm')).*0.01, ...
    'DisplacementUD',table2array(table(:,'Dup_cm')).*nan, ...
    'StationName',table2cell(table(:,'Site')));

%% Convert to local coordinates (okLLtoLocalXY.m)
lon_origin = -120;
DataStruct = okLLtoLocalXY(DataStruct,'Original',lon_origin);

okPlot(DataStruct.Original,'GNSS','scale',100000,'LineWidth',1)
%% Make fault model (okMakeFaultModel.m)
StartXYZ = [220014,3962760,0];
Length = 40000;
Width = 30000;
Strike = 318;
Dip = 90;
PatchStrike = 20;
PatchDip = 15;
FaultModel = okMakeFaultModel(StartXYZ,Length,Width,Strike,Dip,PatchStrike,PatchDip);

okPlot(FaultModel,'FaultModel', ...
    'GNSS',DataStruct.Original, ...
    'scale',1000);

%% Make Green's function (okMakeGreenFunc.m)
RakeDS = -90;
RakeSS = 180;
SlipDS = 0;
SlipSS = 1;
FaultModel = okMakeGreenFunc(DataStruct,'Original','3D',FaultModel,RakeDS,SlipDS,0,'GreenDS');
FaultModel = okMakeGreenFunc(DataStruct,'Original','3D',FaultModel,RakeSS,SlipSS,0,'GreenSS');

okPlot(FaultModel,'GreenFunc', ...
    'dset','GreenSS', ...
    'displ',DataStruct.Original, ...
    'MarkerSize',20, ...
    'title','SS')

okPlot(FaultModel,'GreenFunc', ...
    'dset','GreenDS', ...
    'displ',DataStruct.Original, ...
    'MarkerSize',20, ...
    'title','DS')

%% Make Smoothing matrix (okMakeSmoothMat.m)
FaultModel = okMakeSmoothMat(FaultModel);
okPlot(FaultModel,'Smooth')

%% Invert fault slip (okInvertSlip.m)
SmoothParam = [0.00001,0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,1,10];
GreenFuncPosition = [11,12];
SlipModel = okInvertSlip(DataStruct,'Original',FaultModel,{'GreenSS','GreenDS'},GreenFuncPosition,'SmoothMat', ...
    'solver','nnlsq', ...
    'smoothsearch',SmoothParam);

%% Plot the result
okPlot(SlipModel,'Inversion','slip1',1,'n',4,'title','Strike-slip')
okPlot(SlipModel,'Inversion','lcurve',1,'n',4)
okPlot(SlipModel,'Inversion','residual',1,'n',4,'dset','GNSS','scale',100000)
%% Output 
okOutputGMT(SlipModel,'slip','n',4,'o','Parkfield')
