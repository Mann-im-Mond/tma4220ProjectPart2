classdef HeatEquationSolver
   
    properties
        mesh
        timeInterval
        alpha
        initialValues
        dirichletFunction
        neumannFunction
        
        stiffnessMatrix
        odeMatrix
        neumannVector
        dirichletVector
    end
    
    methods
        function obj = HeatEquationSolver(mesh,timeInterval,alpha,initialValues,dirichletFunction,neumannFunction)
            obj.mesh = mesh;
            obj.timeInterval = timeInterval;
            obj.alpha = alpha;
            obj.initialValues = initialValues;
            obj.dirichletFunction = dirichletFunction;
            obj.neumannFunction = neumannFunction;
        end
        
        function u = solve(obj)
            % TODO this has still to be implemented!
        end 
    end
    
    methods
        function  [] = initializeSystem(obj)
            obj.stiffnessMatrix = obj.getStiffnessMatrix();
            obj.odeMatrix = obj.getODEMatrix();
            obj.neumannVector = obj.getNeumannVector();
            obj.dirichletVector = obj.getDirichletVector();
            obj.removeDirichletBoundary();
        end
        
        function A = getStiffnessMatrix(obj)
            dim = obj.mesh.dim();
            triangles = obj.mesh.triangulation;
            points = obj.mesh.points;
            A = zeros(length(points));
            basicTriangle = BasicTriangle(dim);
            B = basicTriangle.getPDEOfbasisFunction();                      % matrix of partial differential equations

            for cornerIndices = triangles'
                corners = obj.mesh.corners(cornerIndices);
                triangle = Triangle(corners);
                J = triangle.getJacobiFromBasis();
                G = inv(transpose(J))*B;                                    % only a helpingt step towards our submatrix
                subMatrix = transpose(G) * G * abs(det(J))/factorial(dim);  % this is our submatrix
                subMatrix = subMatrix * obj.alpha(triangle.centerOfMass()); % here we scale the submatrix by the termal factor alpha, which belongs to this triangle
                A(cornerIndices,cornerIndices) = ...
                    A(cornerIndices,cornerIndices) + subMatrix;
            end
        end
        
        function M = getODEMatrix(obj)
        
        end
        
        function N = getNeumannVector(obj)
            
        end
        
        function D = getDirichletVector(onj)
            
        end
        
        function [] = removeDirichletBoundary(obj)
            
        end
        
        function u = reinsertDirichletBoundary(obj, uNoBoundary)
            
        end
    end
end
