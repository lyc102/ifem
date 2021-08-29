function [bdNode,bdFacet,isBdNode] = findpolyboundary(varargin)
%% FINDPOLYBOUNDARY finds the boundary facet(s) of a polytopal mesh (dim=2 or 3)
%  2D: facet is given as an NEx2 edge array.
%  3D: facet is given as an NFx4 face array.
%  allFacet contains all facets retrieved from elem{:} or elem
%  If the input is elem, elem has be a cell array
%  Current restriction: face has to be triangular or quadrilateral
%  when f-th face is triangular, the face(f,4) = 0.
% 
%  Examples:
%  
%  totalEdge = uint32([elem(:,[1 2]); elem(:,[2 3]); ...
%                    elem(:,[3 4]); elem(:,[4 1])]);
%  [bdNode,bdFacet,isBdNode] = findpolyboundary(allEdge)
%  
%  allFace = [elem(:,[2, 3, 4]);elem(:,[1, 4,3]);elem(:,[1,2,4]);elem(:,[1,3,2])];
%  allFace(:,4) = 0;
%  [bdNode,bdFacet,isBdNode] = findpolyboundary(allFace)
%  
%  the following example only works in 2D elem is a cell array:
%  elem = {[10,9,11,14,15];[15,14,16,17];[3,1,6,5];[5,6,7,12,11];...
%          [4,3,5,11,9];[2,4,9,10];[12,7,8,13];[14,11,12,13,16]};
%  node = [1 0; 0 0; 0.63591 0; 0.36409 0; 0.68734 0.31329;...
%          1  0.36233; 1 0.63767; 1 1; 0.31266 0.31329;...
%          0 0.36233; 0.49815 0.5; 0.68734 0.68671; 0.63591 1;...
%          0.31266 0.68671;  0 0.63767; 0.36409 1; 0 1];
%  [bdNode,bdFacet,isBdNode] = findpolyboundary(elem);
%  showmeshpoly(node,elem);
%  findnode(node,bdNode);
%  
%  To-do: add bdFlag for more streamlined workflow


%%
if iscell(varargin{1})
    elem = varargin{1};
    elemVertexNumber = cellfun('length',elem);
    maxNv = max(elemVertexNumber);
    minNv = min(elemVertexNumber);
    allFacet = cell(maxNv,1);
    for Nv = minNv:maxNv
        idx = (elemVertexNumber == Nv);
        elemNv = cell2mat(elem(idx));
        NT = sum(idx);
        % # of edge = # of vertices in 2D
        locEdge = [1:Nv; circshift(1:Nv,-1)]';
        allFacet{Nv} = zeros(NT*Nv,2);
        for j = 1:Nv
           allFacet{Nv}((j-1)*NT+1:j*NT,:) = elemNv(:,locEdge(j,:));
        end
    end
    if isrow(allFacet); allFacet = allFacet'; end
    allFacet = cell2mat(allFacet);
elseif ismatrix(varargin{1})
    allFacet = varargin{1};
end
        

[NF, dim] = size(allFacet);
allFacet = sort(allFacet, 2);
[uface, i2, j] = myunique(allFacet); %#ok<*ASGLU>
i1(j(NF:-1:1)) = NF:-1:1;
i1 = i1';
bdFacet = allFacet(i1(i1 == i2), :);
switch dim
    case 2
        isBdNode(bdFacet(:)) = true;
    case 4
        isTriFace = bdFacet(:, 4) == 0;
        isBdNode(bdFacet(isTriFace,1:3)) = true;
        isBdNode(bdFacet(~isTriFace, :)) = true;
end

bdNode = find(isBdNode);
