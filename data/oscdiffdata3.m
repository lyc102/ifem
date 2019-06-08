function pde = oscdiffdata3

pde = struct('f',1,'d',@d);

    function K =  d(p)  % Diffusion constant
    x = p(:,1); y = p(:,2); z = p(:,3);      
    K = 2*(2+sin(10*pi*x).*cos(10*pi*y)*sin(20*pi*z));
    end
end