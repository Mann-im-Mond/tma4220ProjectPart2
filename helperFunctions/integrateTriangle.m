function I = integrateTriangle(p1,p2,p3,f)
    M=transpose([p2-p1 ; p3-p1]);

    phi = @(x,y,i) M(i,1).*x + M(i,2).*y + p1(i);

    fbasis = @(x, y) abs(det(M)) * f(phi(x,y,1),phi(x,y,2));
    frenormed = @(x, y) (1-x).*fbasis(x,y.*(1-x));
    I = integral2(frenormed,0,1,0,1);
end