%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                              Lin,Li-Chieh                               %
%                       Earth and Planetary Sciences                      %
%                   University of California, Riverside                   %
%                               2025.02.26                                %
%                                                                         %
%             ***********************************************             %
%             *** Routine for finite fault slip inversion ***             %
%             ***           okMakeFaultModel.m            ***             %
%             ***********************************************             %
%                                                                         %
% Make discrtized fault patches for running okada slip inversion          %
%                                                                         %
%-------------------------------------------------------------------------%
%                                                                         %
% Input:                                                                  %
% 1. StartXYZ: Numeric. Top-left corner of the fault patch. Z(Depth) <= 0 %
% 2. Length: Numeric. Along-strike length                                 %
% 3. Width: Numeric. Along-dip (down-dip) width                           %
% 4. Strike: Numeric. Strike angle. Follow right-hand rule. Fault dip to  %
%    the right of the fault. North-South stirking fault dips to east will %
%    be 0; East-West stirking fault dips to south will be 90              %
% 5. Dip: Numeric. Dip angle                                              %
% 6. PatchStrike: Numeric. Number of along-strike patches                 %
% 7. PatchDip: Numeric. Numebr of along-dip patches                       %
%                                                                         %
% Example:                                                                %
% An N30E striking fault dips at 45 degrees to the east with 10 patches   %
% along-strike 5 patches along-dip; length=30,000m width=15,000m top-left %
% coordinate is [10 10 0]. All length units should be the same            %
%                                                                         %
% FaultModel = okMakeFaultModel([10 10 0],30000,15000,30,45,10,5)         %
%                                                                         %
%                                                                         %
% Output:                                                                 %
% FaultModel: Structure. Contain the fault geometry                       %
%   FaultModel.okFault: Fault geometry fit for later okada function       %
%      okFault = [X Y Depth Strike Dip PatchLength PatchWidth]            %
%      Coordinates are the patch centroid                                 %
%   FaultModel.Nodes: Coordinates of the top-left corner of the patches   %
%      Nodes = [X Y Depth PatchNumber]                                    %
%   FaultModel.PatchX: X Coordinate of the four corners of the patches for%
%      plotting purposes                                                  %
%   FaultModel.PatchY: Y Coordinate of the four corners of the patches for%
%      plotting purposes                                                  %
%   FaultModel.PatchZ: Z Coordinate of the four corners of the patches for%
%      plotting purposes                                                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FaultModel = okMakeFaultModel(StartXYZ,Length,Width,Strike,Dip,PatchStrike,PatchDip)
% Prepare some lengths and numbers for later calculation
Depth = -Width*sind(Dip) + StartXYZ(3);
DepthInv = (Depth-StartXYZ(3))/PatchDip;
DepthSeg = StartXYZ(3):DepthInv:Depth-DepthInv;
Str = -Strike+90; % To comform with MatLab 
ASPatchSize = Length/PatchStrike;
ADPatchSize = Width/PatchDip;

disp(' ')
disp('******* Constructing fault model okMakeFaultModel.m *******')
disp(strcat('*** Length:',32,num2str(Length)))
disp(strcat('*** Width:',32,num2str(Width)))
disp(strcat('*** Strike angle:',32,num2str(Strike)))
disp(strcat('*** Dip angle:',32,num2str(Dip)))
disp(strcat('*** Along-strike patch count:',32,num2str(PatchStrike)))
disp(strcat('*** Along-dip patch count:',32,num2str(PatchDip)))
disp('*********************************')
disp(strcat('*** Depth (Down-dip):',32,num2str(Depth)))
disp(strcat('*** Along-strike (Length) patch size:',32,num2str(ASPatchSize)))
disp(strcat('*** Along-dip (Width) patch size:',32,num2str(ADPatchSize)))
disp(' ')



% Make Nodes
Nodes = zeros((PatchStrike)*(PatchDip),4);
Xnd = zeros(1,PatchStrike);
Ynd = zeros(1,PatchStrike);
Znd = zeros(1,PatchStrike);
for i = 1:PatchDip
    PatchNum = (PatchStrike)*(i-1)+1:(PatchStrike)*i;
    if i == 1
        % Make the first row at first depth
        Xnd(i) = StartXYZ(1);
        Ynd(i) = StartXYZ(2);
        Znd(i) = StartXYZ(3);
        xinc = Length*cosd(Str)/PatchStrike;
        yinc = Length*sind(Str)/PatchStrike;
        zinc = 0;
        for j = 2:PatchStrike
            Xnd(j) = Xnd(j-1) + xinc;
            Ynd(j) = Ynd(j-1) + yinc;
            Znd(j) = Znd(j-1) + zinc;
        end
    else
        % Subsequent rows propagate from the first row
        % Depth 
        Znd = DepthSeg(i);
    
        % Along-dip coordinate shift
        DepthHori = DepthInv/tand(Dip);
        xdipinc = DepthHori*cosd(Str+90);
        ydipinc = DepthHori*sind(Str+90);
    
        % Calculate the corresponding coordinates
        Xnd = Xnd + xdipinc;
        Ynd = Ynd + ydipinc;
    end
    Nodes((PatchStrike)*(i-1)+1:(PatchStrike)*i,1) = Xnd;
    Nodes((PatchStrike)*(i-1)+1:(PatchStrike)*i,2) = Ynd;
    Nodes((PatchStrike)*(i-1)+1:(PatchStrike)*i,3) = Znd;
    Nodes((PatchStrike)*(i-1)+1:(PatchStrike)*i,4) = PatchNum;
end


% Make Patch and okFault
okFault = zeros((PatchStrike)*(PatchDip),8);
XPatch = zeros(4,(PatchStrike)*(PatchDip));
YPatch = zeros(4,(PatchStrike)*(PatchDip));
ZPatch = zeros(4,(PatchStrike)*(PatchDip));
FaultNum = zeros(1,(PatchStrike)*(PatchDip));
for i = 1:(PatchStrike)*(PatchDip)   
    XincAS = Nodes(2,1) - Nodes(1,1);
    YincAS = Nodes(2,2) - Nodes(1,2);
    XincAD = Nodes(PatchStrike+1,1) - Nodes(1,1);
    YincAD = Nodes(PatchStrike+1,2) - Nodes(1,2);
    ZincAD = Nodes(PatchStrike+1,3) - Nodes(1,3);
    FaultNum(i) = i;
    if (mod(i,PatchStrike) ~= 0) && (i <= PatchStrike*(PatchDip-1))
        % Inner nodes
        XPatch(1,i) = Nodes(i,1);
        XPatch(2,i) = Nodes(i+1,1);
        XPatch(3,i) = Nodes(i+PatchStrike+1,1);
        XPatch(4,i) = Nodes(i+PatchStrike,1);
        YPatch(1,i) = Nodes(i,2);
        YPatch(2,i) = Nodes(i+1,2);
        YPatch(3,i) = Nodes(i+PatchStrike+1,2);
        YPatch(4,i) = Nodes(i+PatchStrike,2);
        ZPatch(1,i) = Nodes(i,3);
        ZPatch(2,i) = Nodes(i+1,3);
        ZPatch(3,i) = Nodes(i+PatchStrike+1,3);
        ZPatch(4,i) = Nodes(i+PatchStrike,3);

    elseif (mod(i,PatchStrike) == 0) && (i <= PatchStrike*(PatchDip-1))
        % Nodes on the far right excluding the bottom
        % Create one node right to it
        XPatch(1,i) = Nodes(i,1);
        XPatch(2,i) = Nodes(i,1) + XincAS;
        XPatch(3,i) = Nodes(i,1) + XincAS + XincAD;
        XPatch(4,i) = Nodes(i+PatchStrike,1);
        YPatch(1,i) = Nodes(i,2);
        YPatch(2,i) = Nodes(i,2) + YincAS;
        YPatch(3,i) = Nodes(i,2) + YincAS + YincAD;
        YPatch(4,i) = Nodes(i+PatchStrike,2);
        ZPatch(1,i) = Nodes(i,3);
        ZPatch(2,i) = Nodes(i,3);
        ZPatch(3,i) = Nodes(i+PatchStrike,3);
        ZPatch(4,i) = Nodes(i+PatchStrike,3);

    elseif (mod(i,PatchStrike) ~= 0) && (i > PatchStrike*(PatchDip-1))
        % Nodes on the bottom excluding the far right ones
        % Create one node underneath it
        XPatch(1,i) = Nodes(i,1);
        XPatch(2,i) = Nodes(i+1,1);
        XPatch(3,i) = Nodes(i+1,1) + XincAD;
        XPatch(4,i) = Nodes(i,1) + XincAD;
        YPatch(1,i) = Nodes(i,2);
        YPatch(2,i) = Nodes(i+1,2);
        YPatch(3,i) = Nodes(i+1,2) + YincAD;
        YPatch(4,i) = Nodes(i,2) + YincAD;
        ZPatch(1,i) = Nodes(i,3);
        ZPatch(2,i) = Nodes(i+1,3);
        ZPatch(3,i) = Nodes(i+1,3) + ZincAD;
        ZPatch(4,i) = Nodes(i,3) + ZincAD;
    elseif (mod(i,PatchStrike) == 0) && (i > PatchStrike*(PatchDip-1))
        % The very bottom-edge node
        % Create one node right and underneath it
        XPatch(1,i) = Nodes(i,1);
        XPatch(2,i) = Nodes(i,1) + XincAS;
        XPatch(3,i) = Nodes(i,1) + XincAS + XincAD;
        XPatch(4,i) = Nodes(i,1) + XincAD;

        YPatch(1,i) = Nodes(i,2);
        YPatch(2,i) = Nodes(i,2) + YincAS;
        YPatch(3,i) = Nodes(i,2) + YincAS + YincAD;
        YPatch(4,i) = Nodes(i,2) + YincAD;

        ZPatch(1,i) = Nodes(i,3);
        ZPatch(2,i) = Nodes(i,3);
        ZPatch(3,i) = Nodes(i,3) + ZincAD;
        ZPatch(4,i) = Nodes(i,3) + ZincAD;
    end

    % okFault
    Xcentroid = Nodes(i,1) + XincAS/2 + XincAD/2;
    Ycentroid = Nodes(i,2) + YincAS/2 + YincAD/2;
    Zcentroid = Nodes(i,3) + ZincAD/2;
    okFault(i,:) = [Xcentroid,Ycentroid,Zcentroid,Strike,Dip,ASPatchSize,ADPatchSize,i];
end

% output to a structure
FaultModel.okFault = okFault;
FaultModel.Nodes = Nodes;
FaultModel.PatchX = XPatch;
FaultModel.PatchY = YPatch;
FaultModel.PatchZ = ZPatch;
FaultModel.StartXYZ = StartXYZ;
FaultModel.Length = Length;
FaultModel.Width = Width;
FaultModel.Strike = Strike;
FaultModel.Dip = Dip;
FaultModel.Depth = Depth;
FaultModel.AlongStrikePatchSize = ASPatchSize;
FaultModel.AlongDipPatchSize = ADPatchSize;
FaultModel.PatchCountStrike = PatchStrike;
FaultModel.PatchCountDip = PatchDip;
FaultModel.TotalPatchCount = PatchStrike*PatchDip;

end





