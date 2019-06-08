function u = SchwarzsmootherDarcy(M,B,f,g,u,elem,edge,elem2edge,freeEdge,itStep)
%% SCHWARZSMOOTHERDARCY Schwarz smoother for Darcy equation
%
% u = SchwarzsmootherDarcy(M,B,f,g,u,elem,edge,elem2edge,freeEdge,itStep)
% implements the overalpping Schwarz smoother for the saddle point system
% arisizing from RT0-P0 discretization of Darcy equation
%
%      |M  B'| |u|  = |f|
%      |B  0 | |p|  = |g|
%
% At each vertex patch, a local problem with prescribed boundary flux is
% solved. Detailed description and convergence analysis can be found in the
% reference.
%
% The negative iteration step, i.e. itStep < 0 means backward Gauss-Seidel
%
% Reference: L. Chen. Multigrid Methods for Constrained Minimization
% Problems and Application to Saddle Point Problems.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
N = max(elem(:));
NT = size(elem,1);
NE = size(edge,1);
% N0 = length(freeNode); % number of interior node
e2v = sparse([1:NE,1:NE], double(edge(1:NE,:)), 1, NE, N);
t2v = sparse([1:NT,1:NT,1:NT], elem(:), 1, NT, N);
oppe2v = sparse(double(elem2edge(:)), elem(:), 1, NE, N);

%% Find free nodes and index map 
isFixedEdge = true(NE,1);
isFixedEdge(freeEdge) = false;
isFreeNode = true(N,1);
isFreeNode(edge(isFixedEdge,:)) = false;
idxMap = zeros(NE,1);
idxMap(freeEdge) = 1:length(freeEdge);

%% Schwarz smoother
for k = 1:abs(itStep)
    if sign(itStep)
        order = 1:N;
    else
        order = N:-1:1;
    end
    for i = order
        if isFreeNode(i) % only solve for freeNode
            starEdgeIdx = idxMap(find(e2v(:,i))); %#ok<*FNDSB>
            starBdEdgeIdx = idxMap(find(oppe2v(:,i)));
            domainBdEdgeIdx = (starBdEdgeIdx == 0);
            starBdEdgeIdx(domainBdEdgeIdx) = [];
            starElemIdx = find(t2v(:,i));
            % solve loc problem
            locM = M(starEdgeIdx,[starEdgeIdx; starBdEdgeIdx]);
            locB = B(starElemIdx,[starEdgeIdx; starBdEdgeIdx]);
            locf = f(starEdgeIdx) - locM*u([starEdgeIdx; starBdEdgeIdx]);
            locg = g(starElemIdx) - locB*u([starEdgeIdx; starBdEdgeIdx]);
            locF = [locf; locg];
            Nf = length(locf);
            Ng = length(locg);
            locM = locM(:,1:Nf);
            locB = locB(:,1:Nf);
            if Nf == Ng  % local problem is unique up to a constant
                activeIdx = 1:(Nf+Ng-1); 
            else
                activeIdx = 1:(Nf+Ng); 
            end
            locL = [locM locB'; locB  zeros(Ng)];
            locep = locL(activeIdx,activeIdx)\locF(activeIdx);
            u(starEdgeIdx) = u(starEdgeIdx) + locep(1:Nf);
            % test inexact solve
%             locLin = [diag(diag(locM)) locB'; locB  zeros(Ng)];
%             locepin = locLin(activeIdx,activeIdx)\locF(activeIdx);
%             locs = locepin(1:Nf);
%             alpha = locf'*locs/(locs'*(locM*locs));
%             u(starEdgeIdx) = u(starEdgeIdx) + alpha*locs;            
        end
    end
end