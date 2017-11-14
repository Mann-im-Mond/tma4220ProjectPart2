addpath('../utils/')
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
gD = @(x) x(1)*0 +0;
gN = @(x) x(1)*0 +0;

%to display the mesh uncommand the following
triplot(triangles,points(:,1),points(:,2));

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



if success
    disp([pad('[unittest/heatEuationSolverTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/heatEuationSolverTest]',40), 'failed!'])
end