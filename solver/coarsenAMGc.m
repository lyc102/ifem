function [isC,As] = coarsenAMGc(A,theta)
%% COARSENAMGC coarsen the graph of A.
%
% isC = coarsenAMGc(A) terturns a logical array to make a set of nodes as
% the coarse ndoes based on As, a strong connectness matrix modified from
% A. The strong connectness is slightly different with the standard
% definition.
%
% [isC,As] = coarsenAMGc(A,theta) accepts the parameter theta to define the
% strong connectness. The default setting is theta = 0.025. It also outputs
% the strong connectness matrix As which could be used in the constrction
% of prolongation and restriction.
%
% Example
%   load lakemesh
%   A = assemblematrix(node,elem);
%   [isC,As] = coarsenAMGc(A);
%
% See also: coarsenAMGrs, interpolationAMGs, amg
%
% Reference page in Help browser
%       <a href="matlab:ifem coarseAMGdoc">coarsenAMGdoc</a> 
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details. 

%% Parameters
if ~exist('theta','var'), theta = 0.025; end
N = size(A,1);
isC = false(N,1);       % C: coarse node
N0 = min(sqrt(N),25);   % number of the coarest nodes

%% Generate strong connectness matrix
Dinv = spdiags(1./sqrt(diag(A)),0,N,N);
Am = Dinv*A*Dinv;       % normalize diagonal of A
[im,jm,sm] = find(Am); 
idx = (-sm > theta);    % delete weakly connect off-diagonal and diagonal
As = sparse(im(idx),jm(idx),sm(idx),N,N); % matrix for strong connectness

%%
% The diagonal of Am is 1. The negative off-diagonal measures the
% diffusivity. The positive off-diagonal is filtered out. 

%% Compute degree of vertex
deg = sum(spones(As)); % number of strongly connected neighbors
deg = full(deg');
if sum(deg>0) < 0.25*sqrt(N)   % too few connected nodes e.g. A is mass matrix
    isC(ceil(rand(N0,1)*N)) = true; % randomly chose N0 nodes
    return                    % smoother is a good preconditioner
end           
idx = (deg>0);
deg(idx) = deg(idx) + 0.1*rand(sum(idx),1); % break the equal degree case so
% that every coarsening can chose around half nodes. Otherwise the
% coarsening could be O(N^2). W-cycle is used to make AMG robust.

%% Find an approximate maximal independent set and put to C set
isF = false(N,1);       % F: fine node
isF(deg == 0) = true;   % isolated nodes are added into F set
isU = true(N,1);        % U: undecided set
while sum(isC) < N/2 && sum(isU) > N0
    % Mark all undecided nodes
    isS = false(N,1);   % S: selected set, changing in the coarsening 
    isS(deg>0) = true;  % deg will be set to zero for nodes in C and F
    S = find(isS); 

    % Find marked nodes with local maximum degree
    [i,j] = find(triu(As(S,S),1));  % i,j and i<j: edges of subgraph S
    idx = deg(S(i)) >= deg(S(j));   % compare degree of connected vertices
    isS(S(j(idx))) = false;         % remove vertices with smaller degree 
    isS(S(i(~idx))) = false;        % if degrees are equal, keep the nodes with smaller index
    isC(isS) = true;                % set selected nodes as coarse nodes

    % Remove coarse nodes and neighboring nodes from undecided set
    [i,j] = find(As(:,isC)); %#ok<*NASGU>
    isF(i) = true;        % neighbor of C nodes are F nodes
    isU = ~(isF | isC);   
    deg(~isU) = 0;        % remove current C and F from the graph  
    
    if sum(isU) <= N0     % add small undecided nodes into C nodes
        isC(isU) = true;
        isU = [];         % to exit the while loop;
    end
end
% fprintf('Number of Coarse Nodes: %6.0u\n',sum(isC));