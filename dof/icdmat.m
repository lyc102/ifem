function D = icdmat(elem2sub,elem2subSign)
%% ICDMAT incidences matrix
%
% D = ICDMAT(elem2sub) assmble the incidence matrix. elem2sub is the index
% pointer from element to sub-simplex.
%
% D = ICDMAT(elem2sub,elem2subSign) assmble the signed incidence matrix.
% elem2sub is the index pointer from element to sub-simplex. elem2subSign
% records the inconsistency of the orientation. Note that D represents the
% differential operator for an appropriate elem2sub.
%
% Note: due to a bug in matlab, the input should be in double type.
%
% Example
%
%     [node, elem] = squaremesh([0 1 0 1],0.25);
%     t2p = icdmat(elem);  % incidence matrix between element and nodes
%     [tempvar,edge] = dofedge(elem);
%     G = icdmat(double(edge),[-1 1]);  % gradient matrix
%     [node, elem] = cubemesh([0 1 0 1 0 1],0.25);
%     elem2face = dof3face(elem);
%     D = icdmat(double(elem2face),[1 -1 1 -1]);  % divergence matrix
%
% See also
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

NT = size(elem2sub,1);
N = max(elem2sub(:));
d = size(elem2sub,2);
% transfer the index to double type but it doesn't work
elem2sub = double(elem2sub); 
if exist('elem2subSign','var'), elem2subSign = double(elem2subSign); end
    
if ~exist('elem2subSign','var')
    D = sparse(repmat((1:NT)',d,1), elem2sub(:), 1, NT, N);
elseif size(elem2subSign,1) == 1
    D = sparse(repmat((1:NT)',d,1), elem2sub(:), ones(NT,1)*elem2subSign,NT,N);
elseif size(elem2subSign,1) == NT    
    D = sparse(repmat((1:NT)',d,1), elem2sub(:), elem2subSign(:), NT, N);
end