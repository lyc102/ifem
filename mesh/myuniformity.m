function u = myuniformity(node,elem,fh,varargin)

N = size(node,1);
NT = size(elem,1);
if size(node,2) == 3
    trep = triangulation(elem,node(:,1),node(:,2),node(:,3));
elseif size(node,2) == 2
    trep = triangulation(elem,node(:,1),node(:,2));
end
[pc,r] = circumcenters(trep);

if ~isnumeric(fh) 
    hc = feval(fh,pc,varargin{:});
else
    switch length(fh)
        case NT
            hc = fh;
        case N
            hc = (fh(elem(:,1))+fh(elem(:,2))+fh(elem(:,3)))/3;
        case 1
            hc = fh*ones(NT,1);
    end
end
sz = r./hc;
u = std(sz)/mean(sz);