function node = meshsmoothing(node,elem,step,rho,method,fixedNode)
%% MESHSMOOTHING improves the geometric mesh quality
%
% node = meshsmoothing(node,elem) improves shape regularity of
% triangles while keeping the density of vertices. 
%
% node = meshsmoothing(node,elem,m) performs m steps mesh smoothing.
% The default setting is m=3.
%
% node = meshsmoothing(node,elem,m,rho) accept a non-uniform density
% given by elementwise function rho. The default choice rho = 1/|t|. The
% quasi-uniform grids corresponds to rho=1. The density function rho can be
% given by a user specified function or a posteriori error estimator in the
% setting of adaptive finite element method.
%
% node = meshsmoothing(node,elem,m,rho,method,fixedNode) accept the
% extra input arguments method and fixedNode. The string |method| is to
% specify the smoothing method: either 'CPT' or 'ODT'. The fixedNode array
% will fix certain nodes, e.g., vertices on a interface. The boundary nodes
% of the mesh will be automatically considered as fixed.
%
% The function meshsmoothing will keep the topology of the input mesh,
% i.e., the node index and connectivity of nodes are unchanged. In
% contrast, the edgeswap function will keep the node but change elem.
%
% The algorithm implemented is the simplified version of ODT smoothing; see
% <a href="http://math.uci.edu/~chenlong/CH2008.html">ODTmesh.</a> 
%
% Example:
%     load airfoilperturbmesh
%     meshquality(node,elem);
%     node = meshsmoothing(node,elem);
%     meshquality(node,elem);
%
% See also  optmesh, bdsmoothing, edgeswap, rmisopoint
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('step','var'), step = 3; end
if ~exist('method','var'), method = 'ODT'; end

%% Smooth boundary nodes
% [node,elem,bdNode] = bdsmoothing(node,elem);

%% Find boundary and interior ndoes
[bdNode,bdEdge,isBdNode,isBdElem] = findboundary(elem); %#ok<*ASGLU>
N = size(node,1);  NT = size(elem,1);

%% Compute mesh quality associated to elements and nodes
qt = 1 - simpqual(node,elem);
qtp = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), [qt, qt, qt], NT, N);
qp = max(qtp);
qp = qp' + 0.01*rand(N,1);  % break the equal case
qp(bdNode) = 0;             % do not include bdNode and fixed nodes
if exist('fixedNode','var'), qp(fixedNode) = 0; end

%% Decompose nodes into independent sets
t2p = sparse([1:NT,1:NT,1:NT], elem(1:NT,:), 1, NT, N);
A = t2p'*t2p;  % connectness of nodes 
nodeSet = coloring(A,5,qp);

%% Move nodes one by one
for m = 1:step
    for k = 1:5  
        % find elements containing nodeSet
        movingNode = nodeSet{k};
        [i,j] = find(t2p(:,movingNode)); %#ok<*NASGU>
        idx = false(NT,1);
        idx(i) = true;
        idx = find(idx);
        % compute centers and areas
        switch upper(method)
          case 'CPT'    
               center = (node(elem(idx,1),:)+node(elem(idx,2),:)+node(elem(idx,3),:))/3;
          case 'ODT'
               center = circumcenter(node,elem(idx,:));
               % modification neary boundary elements: using barycenter to replace circumcenter               
               center(isBdElem(idx),:) = (node(elem(idx(isBdElem(idx)),1),:) ...
                                        + node(elem(idx(isBdElem(idx)),2),:) ...
                                        + node(elem(idx(isBdElem(idx)),3),:))/3;
        end
        % compute the weight
        if exist('rho','var') && ~isempty(rho)
            weight = simplexvolume(node,elem(idx,:));
            if (rho~=1)
                weight = rho(idx).*weight;
            end
        end
        % update to new averaging center
        oldpi = node(movingNode,:);     % keep the location before moving
        if exist('rho','var') && ~isempty(rho) % weighted average
            newnode(:,1) = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)],...
                                repmat(weight.*center(:,1),3,1),[N,1]);
            newnode(:,2) = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)],...
                                repmat(weight.*center(:,2),3,1),[N,1]);
            valence = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)], ...
                                 [weight; weight; weight],[N 1]);            
        else % simple average of centers
            newnode(:,1) = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)],...
                                      [center(:,1);center(:,1);center(:,1)],[N,1]);
            newnode(:,2) = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)],...
                                      [center(:,2);center(:,2);center(:,2)],[N,1]);
            valence = accumarray([elem(idx,1);elem(idx,2);elem(idx,3)],...
                                  ones(3*length(idx),1),[N 1]);
        end
        node(movingNode,:) = newnode(movingNode,:)./[valence(movingNode) valence(movingNode)];
        % check if the moving is valid
        ve2 = node(elem(idx,1),:)-node(elem(idx,3),:);
        ve3 = node(elem(idx,2),:)-node(elem(idx,1),:);
        area = 0.5*(-ve3(:,1).*ve2(:,2)+ve3(:,2).*ve2(:,1));
        invalidElem = idx(area<0);
        if any(invalidElem)
            isNotmoving = false(N,1);
            isNotmoving(elem(invalidElem,:)) = true;
            pidx = isNotmoving(movingNode);
            node(movingNode(pidx),:) = oldpi(pidx,:);
        end
    end
end