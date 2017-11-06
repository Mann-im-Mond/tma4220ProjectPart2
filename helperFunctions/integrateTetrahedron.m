function I = integrateTetrahedron(p0,p1,p2,p3,f)
M=transpose([p1-p0; p2-p0 ; p3-p0]);

phi = @(x,y,z,i) M(i,1).*x + M(i,2).*y +M(i,3).*z + p1(i);

fbasis = @(x, y, z) abs(det(M)) * f(phi(x,y,z,1),phi(x,y,z,2),phi(x,y,z,3));
frenormed = @(x, y, z) (1-x).*(1-y).*(1-x).*fbasis(x,y.*(1-x),z.*(1-x).*(1-y));
I = integral3(frenormed,0,1,0,1,0,1);
end