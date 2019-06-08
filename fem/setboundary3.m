function [bdFlag,face] = setboundary3(node,elem,varargin)
%% SETBOUNDARY3 set type of boundary faces in 3-D.
%
%  bdFlag = setboundary3(node,elem,'Dirichlet') set all boundary edges to
%  Dirichlet type. 
%
%  bdFlag = setboundary3(node,elem,'Neumann') set all boundary edges to
%  Neumann type. 
%
%  bdFlag = setboundary3(node,elem,'Robin') set all boundary edges to
%  Robin type. 
%
%  bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)') set
%  Dirichlet boundary condition on x=1 and x=-1. Other edges are
%  homongenous Neumann boundary condition.
%
%  bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | ...
%  (x==-1)','Neumann','(y==1) | (y==-1)') set
%  Dirichlet boundary condition on x=1 or x=-1 and Neumann boundary
%  condition on y=1 or y=-1.
%
%  bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | ...
%  (x==-1)','Neumann','y==1', 'Robin',' y==-1') set
%  Dirichlet boundary condition on x=1 or x=-1 and Neumann boundary
%  condition on y=1, and Robin boundary condition on y=-1.
%
%  bdFlag = setboundary3(node,elem,'Dirichlet','all','Neumann','y==1') set
%  Neumann boundary condition on y=1 and others are Dirichlet boundary condition.
%
% Example
%
%     node = [-1,-1,-1; 1,-1,-1; 1,1,-1; -1,1,-1; -1,-1,1; 1,-1,1; 1,1,1; -1,1,1]; 
%     elem = [1,2,3,7; 1,6,2,7; 1,5,6,7; 1,8,5,7; 1,4,8,7; 1,3,4,7];
%     bdFlag = setboundary3(node,elem,'Dirichlet','(x==1) | (x==-1)');
%
% See also setboundary
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Find boundary faces
nv = size(elem,2);
if nv == 4 % tetrahedron
    allFace = [elem(:,[2 4 3]);elem(:,[1 3 4]);elem(:,[1 4 2]);elem(:,[1 2 3])];
    nf = 4;
elseif nv == 8 % cube
    allFace = [elem(:,[1 4 3 2]); elem(:,[1 2 6 5]); elem(:,[5 6 7 8]);...
               elem(:,[8 7 3 4]); elem(:,[4 1 5 8]); elem(:,[2 3 7 6])];    
    nf = 6;
end
Nfall = length(allFace);
matlabversion = version;
if str2double(matlabversion(end-5:end-2)) > 2012
    [face, i2, j] = unique(sort(allFace,2),'rows','legacy');
else
    [face, i2, j] = unique(sort(allFace,2),'rows');
end
NT = size(elem,1);
i1(j(Nfall:-1:1)) = Nfall:-1:1; 
i1 = i1';
bdFlag = zeros(Nfall,1,'uint8');
bdFaceidx = i1(i1==i2);

%% Set up boundary faces
nVarargin = size(varargin,2);
if (nVarargin==1)
    bdType = findbdtype(varargin{1});
    bdFlag(bdFaceidx) = bdType;
end
if (nVarargin>=2)
    for i=1:nVarargin/2
        bdType = findbdtype(varargin{2*i-1});
        expr = varargin{2*i};
        if strcmp(expr,'all')
            bdFlag(bdFaceidx) = bdType;
        else
            if nv == 4 % for tetrahedron, a face is triangular
               x = (node(allFace(bdFaceidx,1),1) + ...
                    node(allFace(bdFaceidx,2),1) + ...
                    node(allFace(bdFaceidx,3),1))/3; %#ok<NASGU>
               y = (node(allFace(bdFaceidx,1),2) + ...
                    node(allFace(bdFaceidx,2),2) + ...
                    node(allFace(bdFaceidx,3),2))/3; %#ok<NASGU>
               z = (node(allFace(bdFaceidx,1),3) + ...
                    node(allFace(bdFaceidx,2),3) + ...
                    node(allFace(bdFaceidx,3),3))/3; %#ok<NASGU>
            elseif nv == 8 % for a cube, a face is a quad
               x = (node(allFace(bdFaceidx,1),1) + ...
                    node(allFace(bdFaceidx,2),1) + ...
                    node(allFace(bdFaceidx,3),1) + ...
                    node(allFace(bdFaceidx,4),1))/4; %#ok<NASGU>
               y = (node(allFace(bdFaceidx,1),2) + ...
                    node(allFace(bdFaceidx,2),2) + ...
                    node(allFace(bdFaceidx,3),2) + ...
                    node(allFace(bdFaceidx,4),2))/4; %#ok<NASGU>
               z = (node(allFace(bdFaceidx,1),3) + ...
                    node(allFace(bdFaceidx,2),3) + ...
                    node(allFace(bdFaceidx,3),3) + ...
                    node(allFace(bdFaceidx,4),3))/4; %#ok<NASGU>
            end
           idx = eval(expr);
           bdFlag(bdFaceidx(idx)) = bdType;
        end
    end
end
bdFlag = reshape(bdFlag,NT,nf);
end
%%
function bdType = findbdtype(bdstr)
    switch bdstr
        case 'Dirichlet'
            bdType = 1;
        case 'Neumann'
            bdType = 2;
        case 'Robin'
            bdType = 3;
    end
end