classdef MolecularSurface
    properties
        center
        radius
        number
    end
    methods
        function obj = MolecularSurface(filename)
            data = importdata(filename, ' ');
            obj.center = data(:, 1:3);
            obj.radius = data(:, 4);
            obj.number = size(data, 1);
        end
        function val = phi(obj, p)
            N = size(p, 1);
            vals = zeros(N, obj.number);
            for i = 1:obj.number
               vals(:, i) = obj.phi_i(p, i);
            end
            val= min(vals, [], 2);
        end
        function val = phi_i(obj, p, i)
            pc = bsxfun(@minus, p, obj.center(i, :));
            val = sqrt(sum(pc.^2, 2)) - obj.radius(i);
        end
    end
end
