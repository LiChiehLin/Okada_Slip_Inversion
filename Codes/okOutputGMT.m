%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.03.17                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***              okOutputGMT.m              ***             %
%             ***********************************************             %
%                                                                         %
% Output to GMT plottable format                                          %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. SlipModel: Structure. Slip inversion result from okInvertSlip.m      %
% 2. OutType: Character. Specify the output dataset. See below            %
%    2.1. 'slip': Fault slip distribution. 2 files will be produced:      %
%         "_Solution.txt": The slip inversion parameters                  %
%         "_GMT.txt": GMT compatible file (psxy -L)                       %
%    2.2. 'lcurve': The L-curve search result                             %
%         "_Lcurve.txt"                                                   %
%    2.3. 'displ': The displacement and model prediction                  %
%         "_Displ.txt"                                                    %
%    2.4. 'green': The Okada Green's function used in inversion           %
%         "_Green.txt"                                                    %
%                                                                         %
% Function options:                                                       %
% 'n': Specify which slip inversion (With the search of smoothing)        %
% 'o': Specify the prefix output filename, default: okSlipResult          %
%                                                                         %
% Example:                                                                %
% okOutputGMT(SlipModel,'slip','n',90,'o','Elazig')                       %
% okOutputGMT(SlipModel,'lcurve')                                         %
% okOutputGMT(SlipModel,'displ','n',90)                                   %
% okOutputGMT(SlipModel,'green')                                          %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function okOutputGMT(SlipModel,OutType,varargin)
% parse the options
p = inputParser;
default_n = 1;
default_o = 'okSlipResult';
addParameter(p, 'n', default_n, @(x) isnumeric(x) && length(x)==1);
addParameter(p, 'o', default_o, @(x) ischar(x));


parse(p, varargin{:});
n = p.Results.n;
OutPrefix = p.Results.o;
if strcmp(OutType,'slip')
    SolName = strcat(OutPrefix,'_Solution.txt');
    Slip1GMTName = strcat(OutPrefix,'_Slip1_GMT.txt');
    Slip2GMTName = strcat(OutPrefix,'_Slip2_GMT.txt');
    TotalSlipGMTName = strcat(OutPrefix,'_TotalSlip_GMT.txt');
    % Output:
    % 1. slip inversion report
    % 2. GMT plottable format (psxy -L)
    
    % Retrieve the attributes
    % From slip distribution
    SeisMoment = SlipModel.Moment(n);
    SeisMomentMag = SlipModel.MomentMagnitude(n);
    FaultSlipParse = SlipModel.FaultSlip(:,n);
    TotalSlip = SlipModel.TotalSlip(:,n);
    Rake = SlipModel.Rake(:,n);
    SmoothConst = SlipModel.SmoothParam(n);
    RMSE = SlipModel.RMSE(n);
    Solver = SlipModel.Solver;
    LonOrigin = SlipModel.LongitudeOrigin;


    % From fault geometry
    okFault = SlipModel.FaultModel.okFault;
    Strike = SlipModel.FaultModel.Strike;
    Dip = SlipModel.FaultModel.Dip;
    Length = SlipModel.FaultModel.Length;
    Width = SlipModel.FaultModel.Width;
    StrikePatchSize = SlipModel.FaultModel.AlongStrikePatchSize;
    DipPatchSize = SlipModel.FaultModel.AlongDipPatchSize;
    StrikePatchCount = SlipModel.FaultModel.PatchCountStrike;
    DipPatchCount = SlipModel.FaultModel.PatchCountDip;
    Depth = SlipModel.FaultModel.Depth;
    PatchX = SlipModel.FaultModel.PatchX;
    PatchY = SlipModel.FaultModel.PatchY;
    PatchCount = size(okFault,1);
    Xcentroid = okFault(:,1);
    Ycentroid = okFault(:,2);
    Zcentroid = okFault(:,3);
    Slip1 = FaultSlipParse(1:PatchCount);
    Slip2 = FaultSlipParse(PatchCount+1:PatchCount*2);

    % LocalXY back to LL
    LLcentroid = okLocalXYtoLL(LonOrigin,'localx',Xcentroid,'localy',Ycentroid);
    Loncentroid = LLcentroid(:,1);
    Latcentroid = LLcentroid(:,2);
    PatchLL = okLocalXYtoLL(LonOrigin,'localx',PatchX,'localy',PatchY);
    PatchLon = PatchLL(:,1:PatchCount);
    PatchLat = PatchLL(:,PatchCount+1:end);


    % 1. Write the slip inversion report
    fprintf(strcat('***** Make slip solution file',32,SolName,'\n'))
    fid = fopen(SolName,'w');
    fprintf(fid,strcat('*****************************************************','\n'));
    fprintf(fid,strcat('* Standard output of Okada Slip inversion toolbox   *','\n'));
    fprintf(fid,strcat('* Coordinates are the fault patch centroids         *','\n'));
    fprintf(fid,strcat('*                                                   *','\n'));
    fprintf(fid,strcat('* Note that:                                        *','\n'));
    fprintf(fid,strcat('* 1. The Rake angle is very likely to be wrong. It  *','\n'));
    fprintf(fid,strcat('*    is highly recommended to calculate it again    *','\n'));
    fprintf(fid,strcat('*    based on the dip- and strike-slip you get.     *','\n'));
    fprintf(fid,strcat('* 2. Coordinates are the fault patch centroids      *','\n'));
    fprintf(fid,strcat('*                                                   *','\n'));
    fprintf(fid,strcat('*                                                   *','\n'));
    fprintf(fid,strcat('*        University of California, Riverside        *','\n'));
    fprintf(fid,strcat('* Author:           Li-Chieh Lin                    *','\n'));
    fprintf(fid,strcat('* email:           llin148@ucr.edu                  *','\n'));
    fprintf(fid,strcat('*****************************************************','\n'));
    fprintf(fid,strcat('***** Fault geometry *****','\n'));
    fprintf(fid,strcat('# Length:',32,num2str(Length),'\n'));
    fprintf(fid,strcat('# Width (Down-dip):',32,num2str(Width),'\n'));
    fprintf(fid,strcat('# Depth:',32,num2str(Depth),'\n'));
    fprintf(fid,strcat('# Strike:',32,num2str(Strike),'\n'));
    fprintf(fid,strcat('# Dip:',32,num2str(Dip),'\n'));
    fprintf(fid,strcat('# Along strike patch size:',32,num2str(StrikePatchSize),'\n'));
    fprintf(fid,strcat('# Along dip patch size:',32,num2str(DipPatchSize),'\n'));
    fprintf(fid,strcat('# Along strike patch count:',32,num2str(StrikePatchCount),'\n'));
    fprintf(fid,strcat('# Along dip patch size:',32,num2str(DipPatchCount),'\n'));
    fprintf(fid,strcat('# Total patch count:',32,num2str(PatchCount),'\n'));
    fprintf(fid,strcat('***** Slip inversion parameters *****','\n'));
    fprintf(fid,strcat('# Smoothing constant:',32,num2str(SmoothConst),'\n'));
    fprintf(fid,strcat('# Seismic moment:',32,num2str(SeisMoment),'\n'));
    fprintf(fid,strcat('# Moment magnitude (Mw):',32,num2str(SeisMomentMag),'\n'));
    fprintf(fid,strcat('# Residual (RMSE):',32,num2str(RMSE),'\n'));
    fprintf(fid,strcat('# Inversion solver:',32,Solver,'\n'));
    fprintf(fid,strcat('*************************************************************************************','\n'));
    fprintf(fid,strcat('Lon',32,'Lat',32,'Depth',32,'Slip1',32,'Slip2',32,'TotalSlip',32,'Rake','\n'));
    digit = length(num2str(PatchCount));
    progressStr = strcat('*** Writing patch: #',repmat('0',1,digit));
    fprintf(progressStr)
    for i = 1:PatchCount
        Num = strcat(repmat('0',1,digit-length(num2str(i))),num2str(i));
        progressStr = strcat('*** Writing patch: #',Num);
        fprintf([repmat('\b', 1, length(progressStr)), progressStr]);

        fprintf(fid,strcat(num2str(Loncentroid(i)),32,num2str(Latcentroid(i)),32,num2str(Zcentroid(i)),32,num2str(Slip1(i)),32,num2str(Slip2(i)),32,num2str(TotalSlip(i)),32,num2str(Rake(i)),'\n'));
    end
    fclose(fid);

    % 2. Write the GMT plottable format
    % 2.1. Slip 1 GMT file
    disp(' ')
    fprintf(strcat('***** Make GMT plottable file',32,Slip1GMTName,'\n'))
    digit = length(num2str(PatchCount));
    progressStr = strcat('*** Writing patch: #',repmat('0',1,digit));
    fprintf(progressStr)
    fid = fopen(Slip1GMTName,'w');
    for i = 1:PatchCount
        Num = strcat(repmat('0',1,digit-length(num2str(i))),num2str(i));
        progressStr = strcat('*** Writing patch: #',Num);
        fprintf([repmat('\b', 1, length(progressStr)), progressStr]);

        fprintf(fid,strcat('> -Z',num2str(Slip1(i)),'\n'));
        for j = 1:4
            fprintf(fid,strcat(num2str(PatchLon(j,i)),32,num2str(PatchLat(j,i)),'\n'));
        end
    end
    fclose(fid);
    disp(' ')

    % 2.2 Slip 2 GMT file
    fprintf(strcat('***** Make GMT plottable file',32,Slip2GMTName,'\n'))
    digit = length(num2str(PatchCount));
    progressStr = strcat('*** Writing patch: #',repmat('0',1,digit));
    fprintf(progressStr)
    fid = fopen(Slip2GMTName,'w');
    for i = 1:PatchCount
        Num = strcat(repmat('0',1,digit-length(num2str(i))),num2str(i));
        progressStr = strcat('*** Writing patch: #',Num);
        fprintf([repmat('\b', 1, length(progressStr)), progressStr]);

        fprintf(fid,strcat('> -Z',num2str(Slip2(i)),'\n'));
        for j = 1:4
            fprintf(fid,strcat(num2str(PatchLon(j,i)),32,num2str(PatchLat(j,i)),'\n'));
        end
    end
    fclose(fid);
    disp(' ')

    % 2.3 Total Slip GMT file
    fprintf(strcat('***** Make GMT plottable file',32,TotalSlipGMTName,'\n'))
    digit = length(num2str(PatchCount));
    progressStr = strcat('*** Writing patch: #',repmat('0',1,digit));
    fprintf(progressStr)
    fid = fopen(TotalSlipGMTName,'w');
    for i = 1:PatchCount
        Num = strcat(repmat('0',1,digit-length(num2str(i))),num2str(i));
        progressStr = strcat('*** Writing patch: #',Num);
        fprintf([repmat('\b', 1, length(progressStr)), progressStr]);

        fprintf(fid,strcat('> -Z',num2str(TotalSlip(i)),'\n'));
        for j = 1:4
            fprintf(fid,strcat(num2str(PatchLon(j,i)),32,num2str(PatchLat(j,i)),'\n'));
        end
    end
    fclose(fid);
    disp(' ')


elseif strcmp(OutType,'lcurve')
    LcurveName = strcat(OutPrefix,'_Lcurve.txt');
    Lcurve = SlipModel.Lcurve;
    SmoothConst = SlipModel.SmoothParam;

    fprintf(strcat('***** Make L-curve file',32,LcurveName,'\n'))
    fid = fopen(LcurveName,'w');
    fprintf(fid,'%s %s %s\n','Solution Res.','Model Res.','Smoothing Const.');
    fprintf(fid,'%f %f %e\n',transpose([Lcurve,SmoothConst]));
    fclose(fid);

elseif strcmp(OutType,'displ')
    DisplName = strcat(OutPrefix,'_Displ.txt');
    Displ = SlipModel.Displ;
    LocalX = SlipModel.LocalX;
    LocalY = SlipModel.LocalY;
    Pred = SlipModel.ModelPrediction(:,n);
    LonOrigin = SlipModel.LongitudeOrigin;

    LL = okLocalXYtoLL(LonOrigin,'localx',LocalX,'localy',LocalY);
    Lon = LL(:,1);
    Lat = LL(:,2);

    fprintf(strcat('***** Make displacement and model prediction file',32,DisplName,'\n'))
    fid = fopen(DisplName,'w');
    fprintf(fid,'%s %s %s %s %s %s\n','LocalX','LocalY','Lon','Lat','Displ.','Pred.');
    fprintf(fid,'%f %f %f %f %f %f\n',transpose([LocalX,LocalY,Lon,Lat,Displ,Pred]));
    fclose(fid);

elseif strcmp(OutType,'green')
    GreenName = strcat(OutPrefix,'_Green.txt');
    GreenFunc = sum(SlipModel.GreenFunc,2);
    LocalX = SlipModel.LocalX;
    LocalY = SlipModel.LocalY;
    LonOrigin = SlipModel.LongitudeOrigin;

    LL = okLocalXYtoLL(LonOrigin,'localx',LocalX,'localy',LocalY);
    Lon = LL(:,1);
    Lat = LL(:,2);

    fprintf(strcat('***** Make Greens function file',32,GreenName,'\n'))
    fid = fopen(GreenName,'w');
    fprintf(fid,'%s %s %s %s %s\n','LocalX','LocalY','Lon','Lat','Green');
    fprintf(fid,'%f %f %f %f %f\n',transpose([LocalX,LocalY,Lon,Lat,GreenFunc]));
    fclose(fid);
    
end


end





