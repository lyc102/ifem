function surface = twospherephi

c = [0.51 0.1 0.0;-0.49 -0.2 0.0];
r = [0.6; 0.6];
number = 2;
surface = struct('phi',@phi,'phi_i',@phi_i);
surface.number = number;

    function rhs = phi(p)
        N = size(p, 1);
        vals = zeros(N, number);
        for i = 1:number
            vals(:, i) = phi_i(p, i);
        end
        rhs= min(vals, [], 2);
    end

    function rhs = phi_i(p, i)
        pc = bsxfun(@minus, p, c(i, :));
        rhs = sqrt(sum(pc.^2, 2)) - r(i);
    end
end