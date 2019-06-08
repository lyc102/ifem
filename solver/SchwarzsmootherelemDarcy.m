function u = SchwarzsmootherelemDarcy(M,B,f,g,u,elem,edge,elem2edge,freeEdge,itStep)
%% SCHWARZSMOOTHERELEMDARCY Schwarz smoother (elementwise) for Darcy equation
%
% u = SchwarzsmootherelemDarcy(M,B,f,g,u,elem,edge,elem2edge,freeEdge,itStep)
% implements the overalpping Schwarz smoother for the saddle point system
% arisizing from RT0-P0 discretization of Darcy equation
%
%      |M  B'| |u|  = |f|
%      |B  0 | |p|  = |g|
%
% At each coarse triangle (i.e. 4 small triangles of the current mesh), a
% local problem with prescribed boundary flux is solved. Detailed
% description and convergence analysis can be found in the reference.
%
% The negative iteration step, i.e. itStep < 0 means backward Gauss-Seidel
%
% Reference: L. Chen. Multigrid Methods for Constrained Minimization
% Problems and Application to Saddle Point Problems.
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Data structure
% N = max(elem(:));
NT = size(elem,1);
NTc = NT/4;
NE = size(edge,1);
idxMap = zeros(NE,1);
idxMap(freeEdge) = 1:length(freeEdge);

%% Schwarz smoother
% showmesh(node,elem);
for k = 1:abs(itStep)
    if sign(itStep)
        order = 1:NTc;
    else
        order = NTc:-1:1;
    end
    for i = order
        %     3
        %   / 3 \
        %  5  -  4
        % /1\ 4 /2 \
        %1 -  6 -   2
        locElem = [i NTc+i 2*NTc+i 3*NTc+i];        
        % debug
%         findelem(node,elem,locElem); 
        elemEdgeIdx = idxMap(elem2edge(locElem(4),:)); %#ok<*FNDSB>
        elemBdEdgeIdx = idxMap([elem2edge(locElem(1),[2 3]) elem2edge(locElem(2),[1 3]) ...
                                elem2edge(locElem(3),[1 2])]);
        domainBdEdgeIdx = (elemBdEdgeIdx == 0);
        elemBdEdgeIdx(domainBdEdgeIdx) = [];
%         findedge(node,edge,freeEdge(elemBdEdgeIdx));
        % solve loc problem
        locM = M(elemEdgeIdx,[elemEdgeIdx; elemBdEdgeIdx]);
        locB = B(locElem,[elemEdgeIdx; elemBdEdgeIdx]);
        locf = f(elemEdgeIdx) - locM*u([elemEdgeIdx; elemBdEdgeIdx]);
        locg = g(locElem) - locB*u([elemEdgeIdx; elemBdEdgeIdx]);
        locF = [locf; locg];
        Nf = length(locf);
        Ng = length(locg);
        locM = locM(:,1:Nf);
        locB = locB(:,1:Nf);
        locL = [locM locB'; locB  zeros(Ng)];
        % local problem is unique up to a constant
        activeIdx = 1:(Nf+Ng-1);
        locep = locL(activeIdx,activeIdx)\locF(activeIdx);
        u(elemEdgeIdx) = u(elemEdgeIdx) + locep(1:Nf);
    end
end