%  QUIVER3d *True* 3D quiver plot.
%
%     This function provides an improved technique for visualizing 3D
%     vector fields.  It relies on the efficient use of a single patch
%     call and lighting effects to provide depth cueing.
%
%     QUIVER3D(X,Y,Z,U,V,W) plots velocity vectors as arrows with components
%       (u,v,w) at the points (x,y,z).  The matrices X,Y,Z,U,V,W must all be
%       the same size and contain the corresponding position and velocity
%       components.  QUIVER3D automatically scales the arrows to fit.
%  
%     QUIVER3D(X,Y,Z,U,V,W,COLOR) provides an input argument for coloring
%       the vecotors.  The COLOR argument can be a single color ([1x3] or
%       [3x1]), and indexed color ([1xN] or [Nx1]), or a true color ([3xN or
%       Nx3).  Where N is equal to the number of elements in X or Y or Z or...
%
%     QUIVER3D(X,Y,Z,U,V,W,COLOR,S) automatically
%       scales the arrows to fit and then stretches them by S.
%       Use S=0 to plot the arrows without the automatic scaling.  S=1
%       provides simple autoscaling.
%  
%     QUIVER3D(X,Y,Z,U,V,W,COLOR,S,N) provides an input that defines the
%       density of the tessellation used in generating the arrows.  If N is
%       high (>15) then the rendering time can be long and interactivity of
%       the rendering can be slow.  If N<10 the rendering is generally quick
%       for ~100s of vectors.  The default is N=10.
%
%     H = QUIVER3D(...) returns a patch object of the arrow cluster.
%  
%     Example:
%       quiver3d;  % Generates a sample output of the function.
%
% DBE 10/17/04
% DBE 10/18/04 Fixed bug in arrow diameter scaling.  Fixed bug in arrow length.
function p=quiver3d(X,Y,Z,U,V,W,arrow_color,scaling,N)
% Check the input data size
if nargin==0   % This just generates an example 
  [X,Y] = meshgrid(-2:.2:2,-1:.15:1);
  Z=X.*exp(-X.^2-Y.^2);
  [U,V,W]=surfnorm(X,Y,Z);
  arrow_color=[1 0 0]';  %Z(:)';
  scaling=1;
  N=8;
  figure; hold on; 
    cameratoolbar('Show');
    cameratoolbar('SetMode','orbit','SetCoordSys','none');
    axis equal; axis vis3d;
    axis off; 
    set(gcf,'Color',[0 0 0]);
    l=light; 
    surf(X,Y,Z);
elseif nargin<6
  error('Minimum number of input arguments is 6. MUST include X,Y,Z,U,V, and W data.');
elseif nargin==6
  if ~isequal(size(X),size(Y),size(Z),size(U),size(V),size(W)), error('Dimensions of X,Y,Z,U,V, and W must be the same.'); end
  arrow_color=[1 0 0];
  scaling=1;
  N=10;
elseif nargin==7
  scaling=1;
  N=10;
elseif nargin==8
  N=10;
end
% Linearize the input data...
  X=X(:); Y=Y(:); Z=Z(:);
  U=U(:); V=V(:); W=W(:);
% Input color can be single color vector (1x3 or 3x1) *OR* indexed color
% vector (1xN or Nx1, where N=length(X)) *OR* true color matrix (3xN or Nx3)
% Regardless, the input color has to be repeated for each arrow
single_color =isequal([3 1],size(arrow_color(:)));                                                                                         % Single  color therefore [3x1] or [1x3]...
indexed_color=(size(arrow_color,1)==1 & size(arrow_color,2)==length(X(:))) | (size(arrow_color,1)==length(X(:)) & size(arrow_color,2)==1); % Indexed color therefore [1xN] or [Nx1]...
true_color   =(size(arrow_color,1)==3 & size(arrow_color,2)==length(X(:))) | (size(arrow_color,1)==length(X(:)) & size(arrow_color,2)==3); % True    color therefore [3xN] or [Nx3]...
if single_color
  arrow_color=repmat(arrow_color(:)',[size(X,1) 1]);
elseif indexed_color
  arrow_color=arrow_color(:);                                 
elseif true_color 
  if     isequal([3 3],size(arrow_color)), warning('Ambiguous TRUE COLOR definition intended results *may* require transpose of color matrix.');
  elseif size(arrow_color,1)==3, arrow_color=arrow_color';
  end
else 
  error('Color argument did not appear to be a SINGLE color, nor INDEXED color, nor TRUE color.');
end
[xat,yat,zat,faces0,vertices0]=gen_template_arrow(X,Y,Z,N);
D=0.5*mean((U.^2+V.^2+W.^2).^0.5);   % Calculate an arrow body diameter base on the mean vector magnitude
vertices=[]; faces=[]; fc=[];
for k=1:size(X,1)
    % Scale the normalized arrow data by the vector's length...then autoscale the data
    xa=xat*sqrt(U(k,1).^2+V(k,1).^2+W(k,1).^2);
    ya=yat*D;
    za=zat*D;
    if scaling
      A=get_autoscale_value(X,Y,Z,U,V,W,scaling);
      xa=A*xa;
      ya=A*ya;
      za=A*za;
    end      
  % Generate orthogonal basis for the rotation
    Evct(:,1)=[U(k,1); V(k,1); W(k,1)]/norm([U(k,1); V(k,1); W(k,1)]);  % First unit vector along vector direction
    Evct(:,2)=cross(Evct(:,1),[1 0 0])/norm(cross(Evct(:,1),[1 0 0]));
    Evct(:,3)=cross(Evct(:,1),Evct(:,2));
  % Rotate the template arrow
    XYZ=Evct*[xa(:)'; ya(:)'; za(:)'];
    [xa,ya,za]=deal(reshape(XYZ(1,:),size(xa)),reshape(XYZ(2,:),size(ya)),reshape(XYZ(3,:),size(za)));
  
  % Translate the template arrow
    xa=xa+X(k);  ya=ya+Y(k);  za=za+Z(k);  % x,y,z are the tesselated surface points...
  % Triangulate the surface points
  vertices0=[xa(:),ya(:),za(:)];   
  fc0=repmat(arrow_color(k,:),[size(faces0,1) 1]);  
  
  % Concatenate the patch surfaces for each glyph
  vertices=[vertices; vertices0                     ];
  faces   =[faces;    faces0+(k-1)*size(vertices0,1)];
  fc      =[fc;       fc0                           ];
  
end
% Draw all the glyphs
  p=patch('Vertices',vertices,'Faces',faces,'FaceVertexCData',fc,'FaceColor','Flat');
    set(p,'EdgeColor','None');
    set(p,'BackFaceLighting','lit','FaceLighting','Phong');
    set(p,'SpecularColorReflectance',0.1,'SpecularExponent',1);
    set(p,'DiffuseStrength',0.75);
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SUBROUTINE FUNCTIONS........ %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xa,ya,za,faces,vertices]=gen_template_arrow(X,Y,Z,N);
% Generate a normalized arrow geometry
% L=1;             % Arrow length...
% R=1/10;          % Arrow radius
% S=3/10;          % Arrow slope
L=3/4;             % Arrow length...
R=3/40;            % Arrow radius...
S=1/4;             % Arrow slope...
% Generate the FAUX surface data for tesselation.
% ...can't easily tesselate the true surface data because of the
% ...sudden change in radius at the Tip-Body interface
[xt,yt,zt]=cylinder(linspace(R,0,2),max(N));  % Arrow Tip
[xb,yb,zb]=cylinder([R 2*R],max(N));          % Arrow Body
zb=zb*L;                                      % Scale the body
zt=S*zt+2*L;                                  % Scale and displace the arrow Tip
xx=[zt zb];                                    % Combine the body and top coordinates
yy=[yt yb];                                    % Combine the body and top coordinates
zz=[xt xb];                                    % Combine the body and top coordinates
vertices=[xx(:),yy(:),zz(:)];
faces=convhulln(vertices);
% Generate the REAL surface data for rendering
[xt,yt,zt]=cylinder(linspace(0.7*S,0,2),max(N));  % Arrow tip
[xb,yb,zb]=cylinder(R,max(N));                    % Arrow body
% Scale the data and shift the tip to the end of the arrow body...
zb=zb*L;
zt=S*zt+L;
xa=[zt zb];  % Final arrow coordinates (tip+body)
ya=[yt yb];
za=[xt xb];
return
% The get_auto_scale function uses code that was borrowed from The
% Mathwork's QUIVER3 function.
function autoscale=get_autoscale_value(x,y,z,u,v,w,autoscale)
  % Base autoscale value on average spacing in the x and y
  % directions.  Estimate number of points in each direction as
  % either the size of the input arrays or the effective square
  % spacing if x and y are vectors.
  if min(size(x))==1, n=sqrt(prod(size(x))); m=n; else [m,n]=size(x); end
  delx = diff([min(x(:)) max(x(:))])/n; 
  dely = diff([min(y(:)) max(y(:))])/m;
  delz = diff([min(z(:)) max(y(:))])/max(m,n);
  del = sqrt(delx.^2 + dely.^2 + delz.^2);
  if del>0
    len = sqrt((u/del).^2 + (v/del).^2 + (w/del).^2);
    maxlen = max(len(:));
  else
    maxlen = 0;
  end
  
  if maxlen>0
    autoscale = autoscale*0.9 / maxlen;
  else
    autoscale = autoscale*0.9;
  end
return