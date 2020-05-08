%------------------ PolyMesher  version: 1.1 (Aug13) ---------------------%
% Ref1: C Talischi, GH Paulino, A Pereira, IFM Menezes,                   %
%      "PolyMesher: A general-purpose mesh generator for polygonal        %
%      elements written in Matlab", Struct Multidisc Optim, 2012,         %
%      DOI 10.1007/s00158-011-0706-z                                      %
%                                                                         %
% Ref2: A Pereira, C Talischi, GH Paulino, IFM Menezes, MS Carvalho,      %
%      "Implementation of fluid flow topology optimization in PolyTop",   %
%      Struct Multidisc Optim, 2013, DOI XX.XXXX/XXXXXX-XXX-XXX-X         %
%-------------------------------------------------------------------------%
function [Node,Element,Supp,Load,P] = PolyMesher(Domain,NElem,MaxIter,P)
if ~exist('P','var'), P=PolyMshr_RndPtSet(NElem,Domain); end
NElem = size(P,1);
Tol=5e-6; It=0; Err=1; c=1.5;
BdBox = Domain('BdBox'); PFix = Domain('PFix');
Area = (BdBox(2)-BdBox(1))*(BdBox(4)-BdBox(3));
Pc = P; figure;
while(It<=MaxIter && Err>Tol)
  Alpha = c*sqrt(Area/NElem);
  P = Pc; %Lloyd's update
  R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha);   %Generate the reflections
  [P,R_P] = PolyMshr_FixedPoints(P,R_P,PFix); % Fixed Points 
  [Node,Element] = voronoin([P;R_P]);           %Construct Voronoi diagram
  [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem);
  Area = sum(abs(A));
  Err = sqrt(sum((A.^2).*sum((Pc-P).*(Pc-P),2)))*NElem/Area^1.5;
  fprintf('It: %3d   Error: %1.3e\n',It,Err); It=It+1;
  if NElem<=2000, PolyMshr_PlotMsh(Node,Element,NElem); end; 
end
[Node,Element] = PolyMshr_ExtrNds(NElem,Node,Element);  %Extract node list
[Node,Element] = PolyMshr_CllpsEdgs(Node,Element,0.1);  %Remove small edges
[Node,Element] = PolyMshr_RsqsNds(Node,Element);        %Reoder Nodes
BC=Domain('BC',{Node,Element}); Supp=BC{1}; Load=BC{2}; %Recover BC arrays
PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load);         %Plot mesh and BCs
%------------------------------------------------- GENERATE RANDOM POINTSET
function P = PolyMshr_RndPtSet(NElem,Domain)
P=zeros(NElem,2); BdBox=Domain('BdBox'); Ctr=0;
while Ctr<NElem  
  Y(:,1) = (BdBox(2)-BdBox(1))*rand(NElem,1)+BdBox(1);
  Y(:,2) = (BdBox(4)-BdBox(3))*rand(NElem,1)+BdBox(3);
  d = Domain('Dist',Y);
  I = find(d(:,end)<0);                 %Index of seeds inside the domain
  NumAdded = min(NElem-Ctr,length(I));  %Number of seeds that can be added
  P(Ctr+1:Ctr+NumAdded,:) = Y(I(1:NumAdded),:);
  Ctr = Ctr+NumAdded;
end
%------------------------------------------------------------- FIXED POINTS
function [P,R_P] = PolyMshr_FixedPoints(P,R_P,PFix)
PP = [P;R_P];
for i = 1:size(PFix,1)
  [B,I] = sort(sqrt((PP(:,1)-PFix(i,1)).^2+(PP(:,2)-PFix(i,2)).^2));
  for j = 2:4
    n = PP(I(j),:) - PFix(i,:); n = n/norm(n);
    PP(I(j),:) = PP(I(j),:)-n*(B(j)-B(1));
  end
end
P = PP(1:size(P,1),:); R_P = PP(1+size(P,1):end,:);
%--------------------------------------------------------- REFLECT POINTSET
function R_P = PolyMshr_Rflct(P,NElem,Domain,Alpha)
eps=1e-8; eta=0.9;
d = Domain('Dist',P);  
NBdrySegs = size(d,2)-1;          %Number of constituent bdry segments
n1 = (Domain('Dist',P+repmat([eps,0],NElem,1))-d)/eps;
n2 = (Domain('Dist',P+repmat([0,eps],NElem,1))-d)/eps;
I = abs(d(:,1:NBdrySegs))<Alpha;  %Logical index of seeds near the bdry
P1 = repmat(P(:,1),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,1)
P2 = repmat(P(:,2),1,NBdrySegs);  %[NElem x NBdrySegs] extension of P(:,2)
R_P(:,1) = P1(I)-2*n1(I).*d(I);  
R_P(:,2) = P2(I)-2*n2(I).*d(I);
d_R_P = Domain('Dist',R_P);
J = abs(d_R_P(:,end))>=eta*abs(d(I)) & d_R_P(:,end)>0;
R_P=R_P(J,:); R_P=unique(R_P,'rows');
%---------------------------------------------- COMPUTE CENTROID OF POLYGON
function [Pc,A] = PolyMshr_CntrdPly(Element,Node,NElem)
Pc=zeros(NElem,2); A=zeros(NElem,1);
for el = 1:NElem
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(Element{el}); 
  vxS=vx([2:nv 1]); vyS=vy([2:nv 1]); %Shifted vertices
  temp = vx.*vyS - vy.*vxS;
  A(el) = 0.5*sum(temp);
  Pc(el,:) = 1/(6*A(el,1))*[sum((vx+vxS).*temp),sum((vy+vyS).*temp)];
end
%------------------------------------------------------- EXTRACT MESH NODES
function [Node,Element] = PolyMshr_ExtrNds(NElem,Node0,Element0)
map = unique([Element0{1:NElem}]);
cNode = 1:size(Node0,1);
cNode(setdiff(cNode,map)) = max(map);
[Node,Element] = PolyMshr_RbldLists(Node0,Element0(1:NElem),cNode);
%----------------------------------------------------- COLLAPSE SMALL EDGES
function [Node0,Element0] = PolyMshr_CllpsEdgs(Node0,Element0,Tol)
while(true)
  cEdge = [];
  for el=1:size(Element0,1)
    if size(Element0{el},2)<4, continue; end;  %Cannot collapse triangles
    vx=Node0(Element0{el},1); vy=Node0(Element0{el},2); nv=length(vx);
    beta = atan2(vy-sum(vy)/nv, vx-sum(vx)/nv);
    beta = mod(beta([2:end 1]) -beta,2*pi);
    betaIdeal = 2*pi/size(Element0{el},2);
    Edge = [Element0{el}',Element0{el}([2:end 1])'];
    cEdge = [cEdge; Edge(beta<Tol*betaIdeal,:)];
  end
  if (size(cEdge,1)==0), break; end
  cEdge = unique(sort(cEdge,2),'rows');
  cNode = 1:size(Node0,1);
  for i=1:size(cEdge,1)
    cNode(cEdge(i,2)) = cNode(cEdge(i,1));
  end
  [Node0,Element0] = PolyMshr_RbldLists(Node0,Element0,cNode);
end
%--------------------------------------------------------- RESEQUENSE NODES
function [Node,Element] = PolyMshr_RsqsNds(Node0,Element0)
NNode0=size(Node0,1); NElem0=size(Element0,1);
ElemLnght=cellfun(@length,Element0); nn=sum(ElemLnght.^2); 
i=zeros(nn,1); j=zeros(nn,1); s=zeros(nn,1); index=0;
for el = 1:NElem0
  eNode=Element0{el}; ElemSet=index+1:index+ElemLnght(el)^2;
  i(ElemSet) = kron(eNode,ones(ElemLnght(el),1))';
  j(ElemSet) = kron(eNode,ones(1,ElemLnght(el)))';
  s(ElemSet) = 1;
  index = index+ElemLnght(el)^2;
end
K = sparse(i,j,s,NNode0, NNode0);
p = symrcm(K);
cNode(p(1:NNode0))=1:NNode0;
[Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode);
%------------------------------------------------------------ REBUILD LISTS
function [Node,Element] = PolyMshr_RbldLists(Node0,Element0,cNode)
Element = cell(size(Element0,1),1);
[foo,ix,jx] = unique(cNode);
if ~isequal(size(jx),size(cNode)), jx=jx'; end % +R2013a compatibility fix
if size(Node0,1)>length(ix), ix(end)=max(cNode); end;
Node = Node0(ix,:); 
for el=1:size(Element0,1)
  Element{el} = unique(jx(Element0{el}));
  vx=Node(Element{el},1); vy=Node(Element{el},2); nv=length(vx);
  [foo,iix] = sort(atan2(vy-sum(vy)/nv,vx-sum(vx)/nv));
  Element{el} = Element{el}(iix);
end
%---------------------------------------------------------------- PLOT MESH
function PolyMshr_PlotMsh(Node,Element,NElem,Supp,Load)
clf; axis equal; axis off; hold on;
Element = Element(1:NElem)';                 %Only plot the first block
MaxNVer = max(cellfun(@numel,Element));      %Max. num. of vertices in mesh
PadWNaN = @(E) [E NaN(1,MaxNVer-numel(E))];  %Pad cells with NaN
ElemMat = cellfun(PadWNaN,Element,'UniformOutput',false);
ElemMat = vertcat(ElemMat{:});               %Create padded element matrix
patch('Faces',ElemMat,'Vertices',Node,'FaceColor','w'); pause(1e-6)
if exist('Supp','var')&&~isempty(Supp) %Plot Supp BC if specified
  plot(Node(Supp(:,1),1),Node(Supp(:,1),2),'b>','MarkerSize',8);
end
if exist('Load','var')&&~isempty(Load) %Plot Load BC if specified
  plot(Node(Load(:,1),1),Node(Load(:,1),2),'m^','MarkerSize',8);
end
%-------------------------------------------------------------------------%
%------------------------ PolyMesher - History ---------------------------%
% version: 1.1 (Aug13)
%
% history: Created:    8-Jan-12   Anderson Pereira & Cameron Talischi
%          Supervised by:         Ivan Menezes & Glaucio Paulino
%
%          Modified:   6-Jun-13   Anderson Pereira
%          Created a new function called "PolyMshr_FixedPoints" that
%          allows to specify the exact location of vertices. For more
%          information see Appendix A of Ref2.
%
%          Modified:  14-Aug-13   Tomas Zegard and Sundararajan Natarajan
%          Fixed the changed behaviour of the unique function in
%          Matlab 2013a
%-------------------------------------------------------------------------%