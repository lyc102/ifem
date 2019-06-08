function Ty = gradmap(y0,y1,gamma,My,xistar)
%% GRADMAP generates a graded mesh 


if ~exist('xistar','var'), xistar = 3/4; end
xi = (0:My)'/My;

%% Mapping
Ty = zeros(My+1,1);
Fbound = (y1-y0)/(1-xistar);
if gamma <= Fbound
    Ty = xi.^gamma*(y1-y0) + y0;
else
    ystar = 1/(1+(1-xistar)/xistar*gamma);
    idx = (xi < xistar);
    Ty(idx) = ystar*xi(idx).^gamma/xistar^gamma*(y1-y0) + y0;
    Ty(~idx) = (ystar + (xi(~idx)-xistar)*(1-ystar)/(1-xistar))*(y1-y0) + y0;
end