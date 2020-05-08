%------------------------------ PolyMesher -------------------------------%
% Ref: C Talischi, GH Paulino, A Pereira, IFM Menezes, "PolyMesher: A     %
%      general-purpose mesh generator for polygonal elements written in   %
%      Matlab," Struct Multidisc Optim, DOI 10.1007/s00158-011-0706-z     %
%-------------------------------------------------------------------------%
function [x] = SuspensionDomain(Demand,Arg)
  BdBox = [-2 24 -2 24];
  switch(Demand)
    case('Dist');  x = DistFnc(Arg,BdBox);
    case('BC');    x = BndryCnds(Arg{:},BdBox);
    case('BdBox'); x = BdBox;
    case('PFix');  x = FixedPoints(BdBox);
  end
%----------------------------------------------- COMPUTE DISTANCE FUNCTIONS
function Dist = DistFnc(P,BdBox)
  d1  = dRectangle(P,0,18.885,0,14.56);
  d2  = dLine(P,18.885,1.3030,4,0);
  d3  = dLine(P,3.92,14.56,6.1699,6.88);
  d4  = dLine(P,9.8651,4.0023,18.885,3.70);
  d5  = dLine(P,4,0,0,4);
  d13 = dLine(P,0,14,3.92,14.56);
  d14 = dCircle(P,10,8,4);
  d15 = dLine(P,9.8651,4.0023,6.1699,6.88);
  d   = dDiff(dDiff(dDiff(dDiff(d1,d2),d5),d13),...
        dUnion(dDiff(dIntersect(d3,d4),d15),d14));
  d6  = dCircle(P,2,2,2);
  d7  = dCircle(P,4,2,2);
  d8  = dCircle(P,2,4,2);
  d   = dUnion(d,dUnion(d6,dUnion(d7,d8)));
  d9  = dCircle(P,2,14,2);
  d10 = dCircle(P,2,16,2);
  d11 = dCircle(P,18.885,2.5,1.2);
  d12 = dCircle(P,20,2.5,1.2);
  Dist = dUnion(d,dUnion(d9,dUnion(d10,dUnion(d11,d12))));
%---------------------------------------------- SPECIFY BOUNDARY CONDITIONS
function [x] = BndryCnds(Node,Element,BdBox)
  CornerCircle = sqrt((Node(:,1)-2.0).^2+(Node(:,2)-2.0).^2);
  [foo,CornerCircle] = sort(CornerCircle);
  UpperCircle = sqrt((Node(:,1)- 2.0).^2+(Node(:,2)-16.0).^2);
  [foo,UpperCircle] = sort(UpperCircle);
  Supp = ones(2,3);
  Supp(1,:) = [CornerCircle(1) 1 1];
  Supp(2,:) = [UpperCircle(1)  1 0];
  RightCircle = sqrt((Node(:,1)-20.0).^2+(Node(:,2)-2.5).^2);
  [foo,RightCircle] = sort(RightCircle);
  Load = ones(1,3);
  Load(1,:) = [RightCircle(1) -8 -1];
  x = {Supp,Load};
%----------------------------------------------------- SPECIFY FIXED POINTS
function [PFix] = FixedPoints(BdBox)
  PFix = [ 2   2;
           2  16;
          20 2.5];
%-------------------------------------------------------------------------%