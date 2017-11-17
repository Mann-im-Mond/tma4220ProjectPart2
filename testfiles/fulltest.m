addpath('../utils/')
addpath('../helperFunctions/')
%{
Work with the following easy mesh.
 (0,2)  ___________ (2,2) 
        |\     - /|
        | \  -  / |
 (1,1)  |  X---X  | (1.5,1)
        | /  -  \ |
 (0,0)  |/_____-_\| (2,0)
%}
points = [0,0;0,2;1,1;1.5,1;2,0;2,2];
triangles = [1,2,3;1,3,5;2,3,6;3,4,5;3,4,6;4,5,6];
success = true;
neumannIdentifier = @(x) (x(1)>4);
mesh = FullMesh(triangles, points, neumannIdentifier);
alpha = @(x) x(1)*0 + 1;
timeInterval = TimeInterval(0,10,.1,10);
u0 = @(x) x(1)*0;
gD = @(x,t) norm(x(:)'*[0;1],1)*0 +1 +0*t^2;
gN = @(x,t) norm(x(:)'*[1;1],1)*0 +0 +0*t;

%to display the mesh uncommand the following
%triplot(triangles,points(:,1),points(:,2));

solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN);

solver.initializeSystem();
A = solver.stiffnessMatrix
M = solver.odeMatrix
V = solver.getSummedVectors
u_0 = solver.initialValues
tI = solver.timeInterval;
solver.neumannVector;
solver.dirichletVector;
tspan = [0 10];
%[t,y] = ode45(@(t,y) -A*y+V, tspan, u_0)
%u = solver.solve();
%u = euler(tI,u_0,A,M,V)
%{
for i =1:length(u(1,:))
    u(:,i)
end
%}

function u = euler(tI,u_0,A,M,V)
    t = tI.t_0;
    h = tI.h;
    u=u_0(:);
    i=1;
    while t<tI.t_max
        u = [u,u(:,i)+h*inv(M)*(-A*u(:,i)+V)];
        i=i+1;
        t = t+h;
    end    
end


