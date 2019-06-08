function u = PoissonVEM(node,elem,f)

N = size(node,1);
% NT = size(elem,1);
elemVertexNumber = cellfun('length',elem);
nnz = sum(elemVertexNumber.^2);
ii = zeros(nnz,1);
jj = zeros(nnz,1);
ss = zeros(nnz,1);
b = zeros(N,1);
for Nv = min(elemVertexNumber):max(elemVertexNumber)
    % find polygons with Nv vertices
    idx = find(elemVertexNumber == Nv); % index of elements with k vertices
    NT = length(idx);
    % vertex index and coordinates
    vertex = cell2mat(elem(idx));
    x1 = node(vertex,1);
    y1 = node(vertex,2);
    x2 = circshift(x1,[0 -1]);
    y2 = circshift(y1,[0 -1]);
    % Compute geometry quantity: edge, normal, area, center
    bdIntegral = y1.*x2 - x1.*y2;
    area = sum(bdIntegral,[],2)/2;
    h = sqrt(abs(area));
    cx = sum(x1)/Nv;
    cy = sum(y1)/Nv;
    normVecx = y2 - y1; % normal vector is a rotation of edge vector
    normVecy = x1 - x2; 
    % matrix B, D, I - P
    Bx = (normVecx + circshift(normVecx,1))./(2*h);
    By = (normVecy + circshift(normVecy,1))./(2*h);
    Dx = (x - repmat(cx,1,Nv))./h;
    Dy = (y - repmat(cy,1,Nv))./h;    
    c1 = (1 - (sum(Dx,[],2).*Bx + sum(Dy,[],2).*By))/Nv;
    IminusP = zeros(NT,Nv,Nv);
    for i = 1:Nv
        for j = 1:Nv
            IminusP(:,i,j) = c1(:,j) + Dx(:,i).*Bx(:,j) + Dy(:,i).*By(:,j);
            IminusP(:,i,j) = ones(NT,1) - IminusP(:,i,j);
        end
    end
    % assemble the matrix
    for i = 1:Nv
        for j = 1:Nv
            ii(index+1:index+NT) = vertex(:,i);
            jj(index+1:index+NT) = vertex(:,j);
            ss(index+1:index+NT) = Bx(:,i).*Bx(:,j) + By(:,i).*By(:,j) ...
                                 + IminusP(:,i,:).*IminusP(:,:,j);
            index = index + NT;
        end
    end
    % compute the right hand side
    ft = f(cx,cy)/Nv;
    b = b + accumarry(vertex(:),repmat(ft,Nv,1),[N 1]);
end
A = sparse(ii,jj,ss,N,N);
u = zeros(N,1);
u(freeNode) = A(freeNode,freeNode)\b(freeNode);
