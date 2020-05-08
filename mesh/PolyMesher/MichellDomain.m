%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = MichellDomain(Demand,Arg)
  BdBox = [0 5 -2 2];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dRectangle(P,BdBox(1),BdBox(2),BdBox(3),BdBox(4));
  d2 = dCircle(P,0,0,BdBox(4)/2);
  Dist = dDiff(d1,d2);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  CircleNodes = find(abs(sqrt(Node(:,1).^2+Node(:,2).^2)-1.0)<eps);
  Supp = ones(size(CircleNodes,1),3);
  Supp(:,1) = CircleNodes;
  MidRightFace = sqrt((Node(:,1)-BdBox(2)).^2+...
                      (Node(:,2)-(BdBox(3)+BdBox(4))/2).^2);
  [foo,MidRightFace] = sort(MidRightFace);
  Load = [MidRightFace(1),0,-1];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [5 0];
%-------------------------------------------------------------------------%