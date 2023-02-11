function is = IsSumOnes(v)

c = ones(size(v));

for i = 1:length(v)
    A = nchoosek(1:length(v),i);
    for j = 1:size(A,1)
        Aj = A(j,:);
        cj = c;
        cj(Aj) = -1;
        vj = cj.*v;
        if abs(sum(vj))<10^(-10)
            stp = 1;
        end
    end
end
