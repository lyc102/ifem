function pde = simpledata

pde = struct('f',@f,'g_D',@g_D);

    function z = f(p)
        z = ones(size(p,1),1);
    end

    function z = g_D(p)
        z = zeros(size(p,1),1);
    end
end