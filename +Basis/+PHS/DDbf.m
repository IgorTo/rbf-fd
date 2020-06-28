function res = DDbf(r,p,x)
    res = p.*r.^(p-2) + p*(p-2).*r.^(p-4).*x.^2;
end
