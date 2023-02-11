%AS = matrix.AS;
AS = MatCurl;
Ndof = size(AS,1);
[bc,mapper] = boundaryEdge3D(node,gdof,[1,1,1,1,1,1]);
isBdEdge = true(NEdof,1);
isBdEdge(mapper) = false;

edge = gdof;
%node = femI.p;
%isBdEdge = matrix.isBdEdge;
NE = NEdof;
edgeVec = node(edge(:,2),:) - node(edge(:,1),:);
    edgeLength = sqrt(sum(edgeVec.^2,2));
i1 = (1:NE)'; j1 = double(edge(:,1)); s1 = ones(size(edge,1),1);
i2 = (1:NE)'; j2 = double(edge(:,2)); s2 = -s1;
N = size(node,1);  
G = sparse([i1(:);i2(:)],...
    [j1(:);j2(:)],...
    [s1(:);s2(:)],NE,N);
E = G'*AS*G;
AP = G'*G;
sAP = sum(AP,2);
sE = sum(E,2);