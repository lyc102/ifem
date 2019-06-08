function [rp,nodenew,elemnew] = recoverP02P1(node,elem,p,recoverMethod)
%% recoverP
% Recover the pressure from barycenter to the node.
% The recover methods
% 1) 'LS': least square
% 2) 'LA': solve the Laplacian problem, only this method
%           the new node and new elem will be constructed.
%
% Created by Lin Zhong April, 2013. 

if ~exist('recoveryMethod','var')
    recoverMethod = 'LA';
end
NT = size(elem,1);
N = size(node,1);
rp = zeros(N,1);
% baycenter
xnode = 1/3*(node(elem(:,1),:) + node(elem(:,2),:) + node(elem(:,3),:));
NT2N = sparse(repmat(1:NT,1,3),[elem(:,1);elem(:,2);elem(:,3)],true,NT,N);
 
if strcmp(recoverMethod,'LS') 
    % least square fit for every node
    for i = 1:N
        tempx = xnode(NT2N(:,i),:);
        tempp = p(NT2N(:,i));
        tempn = size(tempx,1);
        if tempn == 1
            rp(i) = tempp;
        elseif tempn == 2
            rp(i) = mean(tempp);
        else
            X = ones(tempn,3);
            X(:,1:2) = tempx;
            coefficient = (X'*X)\(X'*tempp);
            rp(i) = node(i,:)*coefficient(1:2) +coefficient(3);
        end
    end
end

if strcmp(recoverMethod,'LA')
    T = auxstructure(elem);
    edge = T.edge;
    edge2elem = T.edge2elem;
    %% find the boundary edge   
    isbdEdge = (edge2elem(:,1) == edge2elem(:,2));
    freeEdgeidx = find(~isbdEdge);
    bdEdge = edge(isbdEdge,:);
    freeEdge = edge(freeEdgeidx,:);

    NE = size(edge,1);
    NEbd = size(bdEdge,1);
    NEin = NE - NEbd;
    NTnew = 2*NEin;
    elemnew = zeros(NTnew,3);
    nodenew = [node; xnode];
    %% find the corner point
    acubdNode = accumarray([bdEdge(:,1);bdEdge(:,2)], 1, [N 1]);
    bdNodeIdx = find(acubdNode);
    freeNodeIdx = find(~acubdNode);
    bdNode = node(bdNodeIdx);
    NbdNode = size(bdNode,1);
    
    % new elems corresponding to the boundary edges
%     elemnew(1:NEbd,:) = [bdEdge(:,1) bdEdge(:,2) N+edge2elem(isbdEdge,1)];   

    % new elems corresponding to the interior edges
    elemnew(1:NEin,:) = [freeEdge(:,1) N+edge2elem(freeEdgeidx,1) N+edge2elem(freeEdgeidx,2)];
    elemnew(NEin+1:end,:) = [freeEdge(:,2) N+edge2elem(freeEdgeidx,1) N+edge2elem(freeEdgeidx,2)];
    
%   A = assemblematrix(nodenew,elemnew,1);
    %------- construc stiffness matrix-------------
    A11 = sparse(N,N);
    A12 = sparse(N,NT);
    [Dlambda,area] = gradbasis(nodenew,elemnew);
    % the original node 
    ii = double(elemnew(:,1));
    DiDi =area.*dot(Dlambda(:,:,1),Dlambda(:,:,1),2);
    A11 = A11 +sparse(ii,ii,DiDi,N,N);
    % the original node and the new nodes
    for j = 2:3
        jj = double(elemnew(:,j)-N);
        DiDj =area.*dot(Dlambda(:,:,1),Dlambda(:,:,j),2);
        A12 = A12 +sparse(ii,jj,DiDj,N,NT);
    end
    A = [A11 A12];
    %------------------------------------------
    rp(freeNodeIdx) = (-A(freeNodeIdx,N+1:end)*p)./diag(A(freeNodeIdx,freeNodeIdx));
    
    % use least square to recover the boundary nodes
    for i = 1:NbdNode
        ii = bdNodeIdx(i);
        tempx = xnode(NT2N(:,ii),:);
        tempp = p(NT2N(:,ii));
        tempn = size(tempx,1);
        if tempn == 1
            rp(ii) = tempp;
        elseif tempn == 2
            rp(ii) = mean(tempp);
        else
            X = ones(tempn,3);
            X(:,1:2) = tempx;
            coefficient = (X'*X)\(X'*tempp);
            rp(ii) = node(ii,:)*coefficient(1:2) +coefficient(3);
        end
    end
end  
end


