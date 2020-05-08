%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = WrenchDomain(Demand,Arg)
  BdBox = [-0.3 2.5 -0.5 0.5];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1 = dLine(P,0,0.3,0,-0.3);
  d2 = dLine(P,0,-0.3,2,-0.5);
  d3 = dLine(P,2,-0.5,2,0.5);
  d4 = dLine(P,2,0.5,0,0.3);
  d5 = dCircle(P,0,0,0.3);
  d6 = dCircle(P,2,0,0.5);
  douter = dUnion(d6,dUnion(d5,...
           dIntersect(d4,dIntersect(d3,dIntersect(d2,d1)))));
  d7 = dCircle(P,0,0,0.175);
  d8 = dCircle(P,2,0,0.3);
  din = dUnion(d8,d7);
  Dist = dDiff(douter,din);
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  eps = 0.1*sqrt((BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3))/size(Node,1));
  RightCircleNodes = ...
      find(abs(sqrt((Node(:,1)-2).^2+ Node(:,2).^2)-0.3)<eps);   
  Supp = ones(size(RightCircleNodes,1),3);
  Supp(:,1) = RightCircleNodes;
  LeftHalfCircleNodes = ...
      find(abs(max(sqrt(Node(:,1).^2+Node(:,2).^2)-0.175,Node(:,2)))<eps);
  Load = -0.1*ones(size(LeftHalfCircleNodes,1),3);
  Load(:,1) = LeftHalfCircleNodes; Load(:,2) = 0;
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [];
%-------------------------------------------------------------------------%