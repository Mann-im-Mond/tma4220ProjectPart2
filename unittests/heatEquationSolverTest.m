path(pathdef)
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
neumannIdentifier = @(x) (x(1)>.5);
mesh = FullMesh(triangles, points, neumannIdentifier);
alpha = @(x) x(1)*0 +1;
timeInterval = 0;
u0 = @(x) x(1)*0 +0;
gD = @(x,t) norm(x(:)'*[0;1],1) +1 +1*t^2;
gN = @(x,t) norm(x(:)'*[1;1],1) +0 +0*t;

%to display the mesh uncommand the following
%triplot(triangles,points(:,1),points(:,2));

solver = HeatEquationSolver(mesh,timeInterval,alpha,u0,gD,gN);

% --- Check Stiffness Matrix ---
stiffnessMatrix = solver.getStiffnessMatrix();
% I can't give a 100% guarantee for the correctness of the following
% stiffness matrix:
realStiffnessMatrix = [     1.0000,         0,   -1.0000,         0,         0,         0;...
                                 0,    1.0000,   -1.0000,         0,         0,         0;...
                           -1.0000,   -1.0000,    5.5000,   -3.0000,   -0.2500,   -0.2500;...
                                 0,         0,   -3.0000,    6.0000,   -1.5000,   -1.5000;...
                                 0,         0,   -0.2500,   -1.5000,    1.3750,    0.3750;...
                                 0,         0,   -0.2500,   -1.5000,    0.3750,    1.3750...
                                 ];
if(not(isequal(stiffnessMatrix,realStiffnessMatrix)))
    disp('The stiffness matrices are not calculated correctly!')
    success=false; 
end

% --- Check ODE Matrix ---
quadFunc = QuadratFunctionGenerator().getQuadratureFunction(2);
odeMatrix = solver.getODEMatrix();
realODEMatrix = zeros(length(points));
basicTriangle = BasicTriangle(2);
for cornerIndices = triangles'
    corners = points(cornerIndices,:);
    triangle = Triangle(corners);
    invTrafo = triangle.getTrafoToBasis();
    for i=1:3
        for j=1:3
            phi_i = basicTriangle.getBasisFunction(i);
            phi_j = basicTriangle.getBasisFunction(j);
            eval = @(x) phi_i(invTrafo(x)).*phi_j(invTrafo(x));
            c_i = cornerIndices(i);
            c_j = cornerIndices(j);
            realODEMatrix(c_i,c_j)=realODEMatrix(c_i,c_j) + ...
                quadFunc(corners,eval);
        end
    end
end
if(not(equalUpTo(odeMatrix,realODEMatrix,1e-6)))
    disp('The ODE matrices are not calculated correctly!')
    success=false; 
end

% --- Check Neumann Vector ---
neumannVector = solver.getNeumannVector(@(x) solver.neumannFunction(x,0));
realNeumannVector = [0;0;4/3;5/3];
if(not(equalUpTo(neumannVector,realNeumannVector,1e-6)))
    disp('The Neumann Vector is not calculated correctly!')
    success=false; 
end

% --- Check Dirichlet Vector ---
dirichletVector = solver.getDirichletVector(@(x) solver.dirichletFunction(x,0), ...
@(x) solver.derivedDirichletFunction(x,0));
realDirichletVector = [4;0;0;0];
if(not(equalUpTo(dirichletVector,realDirichletVector,1e-6)))
    disp('The Dirichlet Vector is not calculated correctly!')
    success=false; 
end

% --- Check Remove Dirichlet Boundary ---
% I decided not to check this, as the chek would look the same, as the
% function.
solver.initializeSystem();
solver.stiffnessMatrix;
solver.odeMatrix;
solver.neumannVector;
solver.dirichletVector;

% --- Check Reinsert Dirichlet Boundary ---
uWithBoundary = solver.reinsertDirichletBoundary([2;3;5;7]);
realuWithBoundary = [1;3;2;3;5;7];
if(not(equalUpTo(uWithBoundary,realuWithBoundary,1e-6)))
    disp('The Boundary is not inserted correctly afterward!')
    success=false; 
end

if success
    disp([pad('[unittest/heatEuationSolverTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/heatEuationSolverTest]',40), 'failed!'])
end