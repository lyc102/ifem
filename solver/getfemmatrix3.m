function [eqn,T,option] = getfemmatrix3(mesh,pde,fem,var)
%% GETFEMMATRIX generate matrices based on finite element method
%
% [eqn,T] = getfemmatrix(mesh,pde,fem,var) returns an eqn structure for
% matrix equations and T structure for the underlying mesh.
%
% Input arguments
%
%
% mesh.shpae: 'cube', 'Lshape', 'bunny', 'oilpump'
% mesh.type:  'uniform', 'adaptive'
% mesh.size:  default is 1e5
%
% ---------- symmetric and positive definite systems -------------------
%   pde         fem                              var
% 'Poisson'   'P1','P2'                        'jump': jump diffusion
%                                              'osc': oscillate diffusion
% 'Maxwell'   'ND0','ND1','ND2'                'planewave': complex number                                      
%
% ---------- saddle point system -------------------
%   pde         fem                              var
% 'HodgeLap'   'ND0','RT0'                               
%
% See also for examples
%
%     comparesolvers3
%   
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Mesh
N0 = 1;
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = mesh.elem;
    mesh.shape = 'no'; % if the mesh is given, no need to generate
    N0 = size(node,1);
end
if ~isfield(mesh,'size')
    mesh.size = 1e5; % default size is 100,000
end
if ~isfield(mesh,'shape')
    mesh.shape = 'square';
end
if ~isfield(mesh,'type')
    mesh.type = 'uniform';
end
% estimate n0 and refinement level
L = floor(log2(mesh.size/N0)/3); 
n = ceil((mesh.size/8^L)^(1/3)); % number of nodes in one direction
% types of meshes
if strcmp(mesh.type,'adaptive')
    L = L-2;
    switch lower(mesh.shape)
        case 'cube'
            load cubeadaptivemesh;
        case 'lshape'
            load Lshapeadaptivemesh3;
        case 'bunny'
            load bunny
    end
else
    switch lower(mesh.shape)
        case 'cube'
            [node,elem] = cubemesh([-1,1,-1,1,-1,1],2/n);
        case 'lshape'
            [node,elem] = cubemesh([-1,1,-1,1,-1,1],1/n);
            % [node,elem] = delmesh(node,elem,'x>0 & y<0');
            [node,elem] = delmesh(node,elem,'x<0 & y<0 & z>0');
        case 'bunny'
            load bunny;
        case 'oilpump'
            load oilpump;
    end    
end
bdFlag = setboundary3(node,elem,'Dirichlet');
% showboundary3(node,elem);
option.N0 = size(node,1);

%% FEM
option.solver = 'none';
switch upper(pde)
%% Poisson equation
    case 'POISSON'
        if ~exist('var','var')
            Poissonpde = sincosdata;
        else
            switch lower(var)
                case 'jump'
                    Poissonpde = jumpmgdata3;
                case 'osc'
                    Poissonpde = oscdiffdata3;
            end
        end
        switch upper(fem)
            case 'P1'
                for k = 1:L
                       [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
                end
                [soln,eqn,info] = Poisson3(node,elem,bdFlag,Poissonpde,option);
            case 'P2'
                for k = 1:L-1
                       [node,elem,bdFlag] = uniformrefine3(node,elem,bdFlag);
                end                
                [soln,eqn,info] = Poisson3P2(node,elem,bdFlag,Poissonpde,option);                
        end
%% Darcy equation: saddle point system
%     case 'DARCY'
%         if ~exist('var','var')
%             Darcypde = Darcydata1;  % 'isotropic'
%         else
%             switch lower(var)
%                 case 'jump'
%                     Darcypde = Darcydata3;
%                 case 'anisotropic'
%                     Darcypde = Darcydata2;
%             end
%         end
%         switch upper(fem)
%             case 'RT0'
%                 [p,u,eqn,info] = DarcyRT0(node,elem,bdFlag,Darcypde,option);
% %             case 'BDM1'
% %                 [p,u,eqn,info] = DarcyBDM1(node,elem,bdFlag,pde,option);
%         end
%% Maxwell equation        
    case 'MAXWELL'
        if ~exist('var','var')
            Maxwelldata = Maxwelldata2;
        else
            switch  lower(var)
                case 'planewave'
                    Maxwelldata = planewavedata; % plane wave with real coefficients and complex solution
            end
        end
        N0 = size(node,1);
        HB = zeros(N0,4);
        HB(1:N0,1:3) = repmat((1:N0)',1,3); 
        for k = 1:L-1
            [node,elem,bdFlag,HB] = uniformrefine3(node,elem,bdFlag,HB);
        end                                    
        switch upper(fem)
            case 'ND0'
                [u,edge,eqn] = Maxwell(node,elem,HB,Maxwelldata,bdFlag,option); 
            case 'ND1'
                [u,edge,eqn] = Maxwell1(node,elem,HB,Maxwelldata,bdFlag,option);                 
            case 'ND2'
                [u,T,eqn] = Maxwell2(node,elem,HB,Maxwelldata,bdFlag,option); 
                edge = T.edge;
        end
%% Hodge Laplace equation      
    case 'HODGELAP'
        HodgeLapdata = HodgeLaplacian3Edata1;
        switch upper(fem)
            case 'ND0'
                [sigma,u,eqn,info] = HodgeLaplacian3E(node,elem,bdFlag,HodgeLapdata,option);                
            case 'RT0'
                [sigma,u,eqn,info] = HodgeLaplacian3F(node,elem,bdFlag,HodgeLapdata,option);                
        end
%     case 'HELMHOLTZ'
end
%% Output the triangulation
T = struct('node',node,'elem',elem,'bdFlag',bdFlag);
if exist('edge','var'),T.edge = edge; end
if exist('HB','var'),T.HB = HB; end