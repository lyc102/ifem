function [p,t,q,u] = odtmesh2d(fd,fh,h0,hbox,pfix,dispOption,varargin)
%% ODTMESH2D mesh generator based on Optimal Delaunay Triangulation
%
%  odtmesh2d is build upon the distmesh package. Please download a copy
%  from http://www-math.mit.edu/~persson/mesh/ first and add distmesh into
%  the path.
%
%  odtmesh2d shares the same interface with distmesh2d. The following is
%  modified from the M-Lint of distmesh2d
%
%  [p,t,q,u] = odtmesh2d(fd,fh,h0,hbox,pfix,dispOption,varargin)
%      P:         Node positions (Nx2)
%      T:         Triangle indices (NTx3)
%      FD:        Distance function d(x,y)
%      FH:        Scaled edge length function h(x,y)
%      H0:        Initial edge length
%      BBOX:      Bounding box [xmin,ymin; xmax,ymax]
%      PFIX:      Fixed node positions (NFIXx2)
%      OPTION:    plot dispOption
%      FPARAMS:   Additional parameters passed to FD and FH
%
% Example: (Flow past a cyliner. Refined at circle boundary)
%
% % Domain = (0,2.2)*(0,0.41) - {(x,y) | (x-0.2)^2+(y-0.2)^2<=0.05^2}
%
% fd=inline('ddiff(drectangle(p,0,2.2,0,0.41),dcircle(p,0.2,0.2,0.05))','p');
% fh=inline('min(sqrt(sum((p-0.2).^2,2))-(0.05-0.0125),0.346)','p');
% [p,t]=odtmesh2d(fd,fh,0.0035,[0,0;2.2,0.41],[0,0;2.2,0;2.2,0.41;0,0.41]);
%
% More examples can be found in odtmeshdemo2d.
%
%   See also: odtmeshdemo2d, odtmeshopt
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

if ~exist('dispOption','var'), dispOption = 0; end

%% Initial mesh
p = initdistmesh(fd,fh,h0,hbox,pfix,varargin{:});
N = size(p,1);
level = max(min(round(log2(N/300)/2),8),0);
h0 = h0*2^level;
[p,t] = initdistmesh(fd,fh,h0,hbox,pfix,varargin{:});

if dispOption == 1
    showmesh(p,t); pause(0.1)
    fprintf('1. smoothing ')
end
[p,t] = distmeshsmoothing(p,t,fd,fh,pfix,100,varargin{:});
[p,t] = cleanup(p,t,fd,pfix,varargin{:});

if dispOption == 1
    clf; showmesh(p,t);
    meshquality(p,t,fh);
    pause(0.1) 
    fprintf('2. refinement')
end
for k = 1:level
    [p,t] = uniformrefine(p,t);
end
if dispOption == 1
    clf; showmesh(p,t);
    meshquality(p,t,fh);
    pause(0.1)
    fprintf('3. smoothing ')
end
[p,t] = distmeshsmoothing(p,t,fd,fh,pfix,20,varargin{:});
[p,t] = cleanup(p,t,fd,pfix,varargin{:});
if dispOption == 1
    clf; showmesh(p,t);
    meshquality(p,t,fh);
    pause(0.1)
    fprintf('4. odtmeshopt')
end
[p,t] = odtmeshopt(p,t,fd,fh,pfix,3,'odt',varargin{:});
[p,t] = cleanup(p,t,fd,pfix,varargin{:});
if dispOption == 1
    clf; showmesh(p,t); 
end
if dispOption == 1
    [q,u] = meshquality(p,t,fh);
end