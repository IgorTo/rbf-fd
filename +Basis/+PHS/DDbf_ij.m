function res = DDbf_ij(r,p,x,y)
    res = p*(p-2)*r.^(p-4).*x.*y;
end
