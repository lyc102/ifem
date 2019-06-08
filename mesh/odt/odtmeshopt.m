function [p,t] = odtmeshopt(p,t,fd,fh,pfix,step,method,varargin)
%% ODTMESHOPT mesh optimization based on ODT 
%
% Please refer to distmesh for the usage
%
% Copyright (C) Long Chen. See COPYRIGHT.txt for details.

N = size(p,1);  h0 = 1/sqrt(N);  
dptol = 0.001; ttol = 0.1; geps = 0.001*h0; 
deps = sqrt(eps)*h0; alpha = 0.85; 

[bdNode,tempvar,isBdNode] = findboundary(t);
% N = size(node,1); 
% NT = size(t,1);
isInNode = ~isBdNode;
pold = p;

for k = 1:step
    area = abs(simplexvolume(p,t));
    center = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:))/3;
    switch upper(method)
        case 'CPT'    
%            center = (p(t(:,1),:) + p(t(:,2),:) + p(t(:,3),:))/3;
           density = 1./fh(center).^2;
        case 'ODT'
           density = 1./fh(center).^3;
           center = circumcenter(p,t);
           isBdElem = isBdNode(t(:,1)) | isBdNode(t(:,2)) | isBdNode(t(:,3));
           center(isBdElem,:) = (p(t(isBdElem,1),:) + p(t(isBdElem,2),:) ...
                               + p(t(isBdElem,3),:))/3;
    end
    density = density/sum(density);
    area = area.*density;
    force1 = (p(t(:,1),:) - center).*[area area];
    force2 = (p(t(:,2),:) - center).*[area area];
    force3 = (p(t(:,3),:) - center).*[area area];
    dE = zeros(N,2);
    dE(:,1) = accumarray(t(:),[force1(:,1);force2(:,1);force3(:,1)],[N 1]);
    dE(:,2) = accumarray(t(:),[force1(:,2);force2(:,2);force3(:,2)],[N 1]);
    dE(bdNode,:) = 0;
    dE(1:size(pfix,1),:) = 0;                     

    A = sparse(N,N);
    for i = 1:3
        for j= i:3
            if (j==i)
                A = A + sparse(t(:,i),t(:,j),area,N,N);
            else
                A = A + sparse([t(:,i);t(:,j)],[t(:,j);t(:,i)],...
                               -1/2*[area; area],N,N);        
            end        
        end
    end
    dx = zeros(N,2);
    dx(isInNode,:) = -A(isInNode,isInNode)\dE(isInNode,:);
    % 
    % omega = accumarray(t(:),[area;area;area],[N 1]);
    % dx = - dE./[omega omega];   % scale of 2/(d+1)
    p = p + alpha*dx;
%     showmesh(p,t);
%     findnode(p,find(isnan(p)));

    % Bring outside points back to the boundary
    d = feval(fd,p,varargin{:});
    ix = (d>0);               % Find points outside (d>0)
    if find(ix)
    dgradx = (feval(fd,[p(ix,1)+deps,p(ix,2)],varargin{:})-d(ix))/deps; % Numerical
    dgrady = (feval(fd,[p(ix,1),p(ix,2)+deps],varargin{:})-d(ix))/deps; %    gradient
    p(ix,:) = p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];     % Project back to boundary
    end

    % Retriangulation by the Delaunay algorithm
    if (max(sqrt(sum((p-pold).^2,2))/h0)>ttol)            % Any large movement?
        pold = p;
        t = delaunayn(p);                                  % List of triangles
        pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
        t = t(feval(fd,pmid,varargin{:})<-geps,:);         % Keep interior triangles
        t = fixorder(p,t);
%         clf; showmesh(p,t); pause(0.1)
    end

    if max(sqrt(sum(alpha*dx(d<-geps,:).^2,2))/h0)<dptol, break; end

end