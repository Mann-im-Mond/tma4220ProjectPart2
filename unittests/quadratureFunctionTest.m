path(pathdef)
addpath('../utils/')
addpath('../helperFunctions')
quadGen = QuadratFunctionGenerator();
quadr1 = quadGen.getQuadratureFunction(1);
quadr2 = quadGen.getQuadratureFunction(2);
quadr3 = quadGen.getQuadratureFunction(3);

g = @(x) exp(x);
val1 = quadr1([1;2],g);
exactVal1 = integral(g,1,2);

f = @(x) log (x(1)+x(2));
f_ex = @(x,y) log(x + y);
p1=[1,0];
p2=[3,1];
p3=[3,2];
val2 = quadr2([p1;p2;p3], f);
exactVal2 = integrateTriangle(p1,p2,p3,f_ex);

h = @(x) exp(x(1));
h_ex = @(x,y,z) exp(x);
P0=[0,0,0];
P1=[0,2,0];
P2=[0,0,2];
P3=[2,0,0];
val3 = quadr3([P0;P1;P2;P3], h);
exactVal3 = integrateTetrahedron(P0,P1,P2,P3,h_ex);

success = true;

if not(equalUpTo(val1,exactVal1,10e-3))
    disp('The quadradture function for 1D failed!')
 	success = false;
end
if not(equalUpTo(val2,exactVal2,10e-3))
    disp('The quadradture function for 2D failed!')
 	success = false;
end
if not(equalUpTo(val3,exactVal3,10e-3))
    disp('The quadradture function for 3D failed!')
 	success = false;
end

if success
    disp([pad('[unittest/quadratureFunctionTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/quadratureFunctionTest]',40), 'failed!'])
end