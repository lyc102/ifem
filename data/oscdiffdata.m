function pde = oscdiffdata

pde = struct('f',1,'d',@d);

    function K =  d(p)  % Diffusion constant
    x = p(:,1); y = p(:,2);        
    K = 2*(2+sin(10*pi*x).*cos(10*pi*y));
    end
end