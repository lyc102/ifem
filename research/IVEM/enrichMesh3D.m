function mesh = enrichMesh3D(mesh,detail)

%% Usage: mesh structure of a uniform cubic partition
%
% INPUTS:
% mesh --- basic mesh info, contains p and t
% option --- mesh detail level, possible value: 1,2,3.
%   If detail = 1 (sufficient for standard IFE w/o PP), it will contain
%            mesh.e : ne-by-2 matrix stores edges.
%            mesh.f : nt-by-4: 4 faces surrounding each elem
%            mesh.t_e : nt-by-6 matrix stores 6 edges surrounding each cell
%            mesh.t_f : nt-by-4 matrix stores 4 faces surrounding each cell
%   If detail = 2, (sufficient for PPIFE), include all from detail 1, plus
%            mesh.f_t : nf-by-2: two elements sharing each face
%            mesh.f_norm : nf-by-3: unit normal vector of each face
%            mesh.f_e : nf-by-3: edges surrounding each face
% OUTPUTS:
% mesh --- a struct data contains mesh information.
%
% Last Modified: 08/07/2020 by Xu Zhang

%% 0. Check Inputs
if nargin == 1
    detail = 3; % highest level detail
end
p = mesh.p; t = mesh.t; nt = size(t,1);

%% 1. Generate e and t_e (Detail Level = 1)
eAll = [t(:,[1,2]);t(:,[1,3]);t(:,[1,4]);t(:,[2,3]);t(:,[2,4]);t(:,[3,4])];
% the local edge order is (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
[e,~,ic] = unique(sort(eAll,2),'rows','stable');
t_e = reshape(ic,nt,6);
% Generate f and t_f
fAll = [t(:,[1,2,3]);t(:,[1,2,4]);t(:,[1,3,4]);t(:,[2,3,4])];
[f,~,ic] = unique(sort(fAll,2),'row','stable');
t_f = reshape(ic,nt,4);

mesh.e = e; mesh.f = f; mesh.t_e = t_e; mesh.t_f = t_f;

if detail >= 2
    %% 2. Generate f_t (Detail Level = 2)
    nf = size(f,1); one = (1:nt)';
    f_t1 = zeros(nf,1);
    f_t1(t_f(:,1)) = f_t1(t_f(:,1)) + one;
    f_t1(t_f(:,2)) = f_t1(t_f(:,2)) + one;
    f_t1(t_f(:,3)) = f_t1(t_f(:,3)) + one;
    f_t1(t_f(:,4)) = f_t1(t_f(:,4)) + one;
    
    f_t2 = zeros(nf,1);
    f_t2(t_f(:,1)) = f_t2(t_f(:,1)) + one;
    id = f_t2(t_f(:,2)) == 0;
    f_t2(t_f(id,2)) = one(id);
    id = f_t2(t_f(:,3)) == 0;
    f_t2(t_f(id,3)) = one(id);
    id = f_t2(t_f(:,4)) == 0;
    f_t2(t_f(id,4)) = one(id);
    f_t = [f_t2,f_t1-f_t2];
    
    %% 3. Generate f_norm
    n = cross(p(f(:,1),:) - p(f(:,3),:),p(f(:,2),:) - p(f(:,3),:));
    f_norm = n./(sqrt(n(:,1).^2+n(:,2).^2+n(:,3).^2));
    % Adjust normal direction so it is outward to the 1st element in f_t
    ctF = 1/3*(p(f(:,1),:) + p(f(:,2),:) + p(f(:,3),:));
    t1 = f_t(:,1);
    ctT = 1/4*(p(t(t1,1),:) + p(t(t1,2),:) + p(t(t1,3),:) + p(t(t1,4),:));
    flag = sum((ctF-ctT).*f_norm,2);
    f_norm((flag < 0),:) = -1*f_norm((flag < 0),:);
    
    %% 4. Generate f_e
    feAll = [t_e(:,[1,2,4]);t_e(:,[1,3,5]);t_e(:,[2,3,6]);t_e(:,[4,5,6])];
    [f_e,~,~] = unique(sort(feAll,2),'row','stable');
    
    mesh.f_t = f_t; mesh.f_e = f_e; mesh.f_norm = f_norm;
end

if detail >= 3
    %% 5. may add more mesh detail later.
end
