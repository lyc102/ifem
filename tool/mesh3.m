function mesh3(x,y,z,alpha)
%% MESH3 plot 3-dimensional tensor product grids
%
%
% Example
%  [z,x,y] = ndgrid(0:.2:1, 0:.25:2, 0:0.5:3);
%  mesh3(x,y,z);

if ~exist('alpha','var'), alpha = 1e-3; end

col = [1-alpha 1-alpha 1-alpha];

startx = x(:,:,1);
endx = x(:,:,end);
starty = y(:,:,1);
endy = y(:,:,end);
startz = z(:,:,1);
endz = z(:,:,end);
plot3([startx(:) endx(:)]',[starty(:) endy(:)]',[startz(:) endz(:)]','color',col);
hold on
plot3([startx(end,:); endx(end,:)],[starty(end,:); endy(end,:)],[startz(end,:); endz(end,:)],'k','linewidth',1);
plot3([startx(:,1) endx(:,1)]',[starty(:,1) endy(:,1)]',[startz(:,1) endz(:,1)]','k','linewidth',1);

startx = squeeze(x(1,:,:));
endx = squeeze(x(end,:,:));
starty = squeeze(y(1,:,:));
endy = squeeze(y(end,:,:));
startz = squeeze(z(1,:,:));
endz = squeeze(z(end,:,:));
plot3([startx(:) endx(:)]',[starty(:) endy(:)]',[startz(:) endz(:)]','color',col);
plot3([startx(1,:); endx(1,:)],[starty(1,:); endy(1,:)],[startz(1,:); endz(1,:)],'k','linewidth',1);
plot3([startx(:,1) endx(:,1)]',[starty(:,1) endy(:,1)]',[startz(:,1) endz(:,1)]','k','linewidth',1);

startx = squeeze(x(:,1,:));
endx = squeeze(x(:,end,:));
starty = squeeze(y(:,1,:));
endy = squeeze(y(:,end,:));
startz = squeeze(z(:,1,:));
endz = squeeze(z(:,end,:));
plot3([startx(:) endx(:)]',[starty(:) endy(:)]',[startz(:) endz(:)]','color',col);
plot3([startx(end,:); endx(end,:)],[starty(end,:); endy(end,:)],[startz(end,:); endz(end,:)],'k','linewidth',1);
plot3([startx(:,1) endx(:,1)]',[starty(:,1) endy(:,1)]',[startz(:,1) endz(:,1)]','k','linewidth',1);
hold off
axis equal