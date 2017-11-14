addpath('../utils/')
%{
Work with the following easy mesh.
 (0,2)  ___________ (2,2) 
        |\       /|
        |  \   /  |
        |    X    | (1,1)
        |  /   \  |
 (0,0)  |/_______\| (2,0)
%}
points = [0,0;2,0;1,1;0,2;2,2];
triangles = [1,2,3;1,3,4;2,3,5;3,4,5];
success = true;
neumannIdentifier = @(x) (x(1)>.5);
mesh = FullMesh(triangles, points, neumannIdentifier);
boundaryEdges = mesh.boundaryEdges();
realBoundaryEdges = [1,2;1,4;2,5;4,5];
realDirichletNodes = [1;4];
dirichletNodes = mesh.dirichletBoundaryNodes;
realNeumannFaces = [2,5];
neumannFaces = mesh.neumannBoundaryFaces;

% TODO we can summarize the following in a more general test function.

%Test whether we found the correct boundary edges.
for edge=boundaryEdges'
    C=intersect(perms(edge),realBoundaryEdges,'rows');
    if isempty(C)
        disp('We found wrong boundary edges!')
        success = false;
    end
end
for edge=realBoundaryEdges'
    C=intersect(perms(edge),boundaryEdges,'rows');
    if isempty(C)
        disp('We did not find all boundary edges!')
        success = false;
    end
end

%Test whether we found the correct Neumann faces.
for edge=neumannFaces'
    C=intersect(perms(edge),realNeumannFaces,'rows');
    if isempty(C)
        disp('We found wrong Neumann faces!')
        success = false;
    end
end
for edge=realNeumannFaces'
    C=intersect(perms(edge),neumannFaces,'rows');
    if isempty(C)
        disp('We did not find all Neumann faces!')
        success = false;
    end
end

%Test whether we found the correct Dirichlet nodes.
if(not(and(length(dirichletNodes) == length(realDirichletNodes), ...
    length(intersect(dirichletNodes,realDirichletNodes) ...
    ==length(realDirichletNodes)))))
    disp('The dirichlet nodes are not found correctly.')
    success = false;
end


if success
    disp([pad('[unittest/fullMeshTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/fullMeshTest]',40), 'failed!'])
end