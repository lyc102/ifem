function [eqn,T] = getfemmatrix(mesh,pde,fem,var)
%% GETFEMMATRIX generate matrices based on finite element method
%
% [eqn,T] = getfemmatrix(mesh,pde,fem,pdedata) returns an eqn structure for
% matrix equations and T structure for the underlying mesh.
%
% Input arguments
%
%
% mesh.shpae: 'square', 'Lshape', 'circle', 'lake'
% mesh.type:  'uniform', 'adaptive'
% mesh.size:  default is 1e5
%
% ---------- symmetric and positive definite systems -------------------
%   pde         fem                              var
% 'Poisson'   'P1','P2','P3'                   'jump': jump diffusion
%                                              'osc': oscillate diffusion
% 'Maxwell'   'ND0'                                      
%
% ---------- saddle point system -------------------
%   pde         fem                              var
% 'Stokes'    'P2P0','P2P1','CRP0','P1BP1' 
%             'isoP2P0','isoP2P1'
%             'RTP0','BDMP0'                               
% 'Darcy'     'RT0','BDM1'                     'jump': random jump diffusion
%                                              'anisotropic':  tensor
% 'Biharmonic'  'P1','P2','P3'                   
%
%
% See also for examples
%
%     comparesolvers
%   
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

%% Mesh
if isfield(mesh,'node') && isfield(mesh,'elem')
    node = mesh.node;
    elem = mesh.elem;
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
L = floor(log2(mesh.size)/2); 
n = ceil(sqrt(mesh.size/4^L)); % number of nodes in one direction
% types of meshes
if strcmp(mesh.type,'adaptive')
    switch lower(mesh.shape)
        case 'square'
            load('squareadaptivemesh','node','elem');
        case 'lshape'
            load('Lshapeadaptivemesh','node','elem');
        case 'circle'
            load('circleadaptivemesh','node','elem');
        case 'lake'
            load('lakemesh','node','elem');
    end
else
    switch lower(mesh.shape)
        case 'square'
            [node,elem] = squaremesh([0 1 0 1], 1/n);
        case 'lshape'
            [node,elem] = squaremesh([-1,1,-1,1],2/n);
            [node,elem] = delmesh(node,elem,'x>0 & y<0');
        case 'circle'
            [node,elem] = circlemesh(0,0,1,3/n);
        case 'lake'
            load('lakemesh','node','elem');
    end    
end
bdFlag = setboundary(node,elem,'Dirichlet');
showmesh(node,elem);

%% FEM
option.solver = 'none';
switch upper(pde)
%% Poisson equation
    case 'POISSON'
        if ~exist('pdedata','var')
            Poissonpde = sincosdata;
        else
            switch lower(var)
                case 'jump'
                    Poissonpde = Kelloggdata;
                case 'osc'
                    Poissonpde = oscdiffdata;
            end
        end
        switch upper(fem)
            case 'P1'
                for k = 1:L
                       [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end
                [soln,eqn,info] = Poisson(node,elem,bdFlag,Poissonpde,option);
            case 'P2'
                for k = 1:L-1
                       [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end                
                [soln,eqn,info] = PoissonP2(node,elem,bdFlag,Poissonpde,option);                
            case 'P3'
                for k = 1:L-2
                       [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end
                [soln,eqn,info] = PoissonP3(node,elem,bdFlag,Poissonpde,option);                                
        end
%% Stokes equation: saddle point system        
    case 'STOKES'
        if ~exist('pdedata','var')
            Stokespde = Stokesdata0;
        else
            switch lower(var)
                case 'jump'
                    Stokespde = Stokesdata0; % to add an interface one
            end
        end
        for k = 1:L-1
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        end                            
        switch upper(fem)
            % all elements lead to system with similar sizes
            case 'P2P1'
                [soln,eqn,info] = StokesP2P1(node,elem,bdFlag,Stokespde,option);
            case 'P2P0'
                [soln,eqn,info] = StokesP2P0(node,elem,bdFlag,Stokespde,option);
            case 'ISOP2P1'
                [soln,eqn,info] = StokesisoP2P1(node,elem,bdFlag,Stokespde,option);
            case 'ISOP2P0'
                [soln,eqn,info] = StokesisoP2P0(node,elem,bdFlag,Stokespde,option);
            case 'CRP0'
                [soln,eqn,info] = StokesCRP0(node,elem,bdFlag,Stokespde,option);
            case 'CRP1'
                [soln,eqn,info] = StokesCRP1(node,elem,bdFlag,Stokespde,option);
            case 'P1BP1'
                [soln,eqn,info] = StokesP1bP1(node,elem,bdFlag,Stokespde,option);
            case 'RTP0'
                [soln,eqn,info] = StokesRT0(node,elem,bdFlag,Stokespde,option);        
            case 'BDMP0'
                [soln,eqn,info] = StokesBDM1B(node,elem,bdFlag,Stokespde,option);                    
        end
%% Darcy equation: saddle point system
    case 'DARCY'
        if ~exist('pdedata','var')
            Darcypde = Darcydata1;  % 'isotropic'
        else
            switch lower(var)
                case 'jump'
                    Darcypde = Darcydata3;
                case 'anisotropic'
                    Darcypde = Darcydata2;
            end
        end
        switch upper(fem)
            case 'RT0'
                [p,u,eqn,info] = DarcyRT0(node,elem,bdFlag,Darcypde,option);
%             case 'BDM1'
%                 [p,u,eqn,info] = DarcyBDM1(node,elem,bdFlag,pde,option);
        end
%% Maxwell equation        
    case 'MAXWELL'
        if ~exist('pdedata','var')
            Maxwelldata = eddycurrentdata1;
        else
            Maxwelldata = eddycurrentdata1; % add more cases
        end
        for k = 1:L-1
            [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
        end                                    
        switch upper(fem)
            case 'ND0'
                [u,edge,eqn] = eddycurrent(node,elem,bdFlag,Maxwelldata,option);
%             case 'BDM1'
%                 [p,u,eqn,info] = DarcyBDM1(node,elem,bdFlag,pde,option);
        end
      
%% Biharmonic equation: saddle point system
    case 'BIHARMONIC'
        if ~exist('pdedata','var')
            Biharmonicpde = biharmonicdata;
        else
            Biharmonicpde = biharmonicdata; % more cases            
        end
        switch upper(fem)
            case 'P1'
                for k = 1:L
                    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end
                [soln,eqn,info] = biharmonicP1(node,elem,bdFlag,Biharmonicpde,option);
            case 'P2'
                for k = 1:L-1
                    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end                
                [soln,eqn,info] = biharmonicP2(node,elem,bdFlag,Biharmonicpde,option);                
            case 'P3'
                for k = 1:L-2
                    [node,elem,bdFlag] = uniformrefine(node,elem,bdFlag);
                end
                [soln,eqn,info] = biharmonicP3(node,elem,bdFlag,Biharmonicpde,option);                                
        end        
%     case 'HODGELAP'
%     case 'HELMHOLTZ'
end
%% Output the triangulation
T = struct('node',node,'elem',elem,'bdFlag',bdFlag);