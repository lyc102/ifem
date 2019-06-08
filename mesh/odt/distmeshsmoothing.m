function [p,t] = distmeshsmoothing(p,t,fd,fh,pfix,step,varargin)
%% DISTMESHSMOOTHING mesh smoothing using distmesh
%
% Please refer to distmesh

% parameters
N = size(p,1);  h0 = 1.25/sqrt(N);
dptol=.001; deltat = 0.25; ttol = 0.1; geps = 0.0001*h0; 
Fscale = 1.2; deps=sqrt(eps)*h0;

pold = inf;                                            % For first iteration
for k = 1:step
    % Retriangulation by the Delaunay algorithm
    if (max(sqrt(sum((p-pold).^2,2))/h0)>ttol)            % Any large movement?
        pold = p;                                          % Save current positions
        t = delaunayn(p);                                  % List of triangles
        pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
        t = t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
        t = fixorder(p,t);
        bars = [t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
        matlabversion = version;
        if str2double(matlabversion(end-5:end-2)) > 2012
            bars = unique(sort(bars,2),'rows','legacy');                % Bars as node pairs
        else        
            bars = unique(sort(bars,2),'rows');                % Bars as node pairs
        end
%         clf; showmesh(p,t); pause(0.1)
    end

    % Move mesh points based on bar lengths L and forces F
    barvec = p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
    L = sqrt(sum(barvec.^2,2));                          % L = Bar lengths
    hbars = feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
    L0 = hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
    F = max(L0-L,0);                                     % Bar forces (scalars)
    Fvec = F./L*[1,1].*barvec;                           % Bar forces (x,y components)
    Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:size(pfix,1),:) = 0;                          % Force = 0 at fixed points
    p = p + deltat*Ftot;                                   % Update node positions

    % Bring outside points back to the boundary
    area = abs(simplexvolume(p,t));
    patchArea = sparse(t,ones(size(t,1),3),area*[1,1,1],N,1);
    ideaHeight = sqrt(sqrt(3)*patchArea/6);
    d = feval(fd,p,varargin{:});
    ix = ((d>0) | (d+0.25*ideaHeight)>0);      % nodes out side of domain or close to boundary
%     d = feval(fd,p,varargin{:}); 
%     ix = d>0;               % Find points outside (d>0)
    dgradx = (feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
    dgrady = (feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; %    gradient
    p(ix,:) = p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];     % Project back to boundary

    if max(sqrt(sum(deltat*Ftot(d<-geps,:).^2,2))/h0)<dptol, break; end

end
t = delaunayn(p);                                  % List of triangles
pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
t = t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
t = fixorder(p,t);