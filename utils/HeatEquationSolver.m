classdef HeatEquationSolver
   
    properties
        mesh
        timeInterval
        alpha
        initialValues
        dirichletFunction
        derivedDirichletFunction
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
            % TODO we have to calculate this derivative function.
            obj.derivedDirichletFunction = @(x) x(1)*0;
            obj.neumannFunction = neumannFunction;
        end
        
        function u = solve(obj)
            % TODO this has still to be implemented!
        end 
    end
    
    methods
        function  obj = initializeSystem(obj)
            obj.stiffnessMatrix = obj.getStiffnessMatrix();
            obj.odeMatrix = obj.getODEMatrix();
            obj.neumannVector = obj.getNeumannVector();
            obj.dirichletVector = obj.getDirichletVector();
            obj = obj.removeDirichletBoundary();
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
                subMatrix = transpose(G) * G * abs(det(J))/factorial(dim)...  % this is our submatrix
                    * obj.alpha(triangle.centerOfMass()); % here we scale the submatrix by the termal factor alpha, which belongs to this triangle
                A(cornerIndices,cornerIndices) = ...
                    A(cornerIndices,cornerIndices) + subMatrix;
            end
        end
        
        function M = getODEMatrix(obj)
            dim = obj.mesh.dim();
            triangles = obj.mesh.triangulation;
            points = obj.mesh.points;
            M = zeros(length(points));
            basicTriangle = BasicTriangle(dim);
            
            productIntegrals = basicTriangle.productBasisFunctionIntegralMatrix();
            
            for cornerIndices = triangles'
                corners = obj.mesh.corners(cornerIndices);
                triangle = Triangle(corners);
                J = triangle.getJacobiFromBasis();
                alphafactor = obj.alpha(triangle.centerOfMass());
                M(cornerIndices,cornerIndices) = ...
                    M(cornerIndices,cornerIndices) + ...
                    abs(det(J))*productIntegrals*alphafactor;
            end
        end
        
        function N = getNeumannVector(obj)
            points = obj.mesh.points;
            N = zeros(length(points),1);
            neumannFaces = obj.mesh.neumannBoundaryFaces;
            if isempty(neumannFaces)
                return
            end
            
            subdim = length(neumannFaces(1,:))-1;
            basicTriangle = BasicTriangle(subdim);
            quad = QuadratFunctionGenerator().getQuadratureFunction(subdim);
            basicCorners = basicTriangle.getCornerPoints();

            for face = neumannFaces'
                corners = points(face',:);
                triangle = Triangle(corners);
                trafo = triangle.getTrafoFromBasis();
                transformedNeumann = @(x) obj.neumannFunction(trafo(x));
                volumeFace = triangle.getVolume();
                volumeBasic = 1/(factorial(subdim));
                quot = abs(volumeFace/volumeBasic);
                alphafactor = obj.alpha(triangle.centerOfMass());
                
                % add the amound of the integral over this basis element

                for i=1:(subdim+1)
                    phi = basicTriangle.getBasisFunction(i);
                    idx = face(i);
                    N(idx) = N(idx) + ...
                        quad(basicCorners,@(x) transformedNeumann(x).*...
                        phi(x).*quot.*alphafactor);
                end
            end
        end
        
        function D = getDirichletVector(obj)
            A = obj.getStiffnessMatrix();
            M = obj.getODEMatrix();
            points = obj.mesh.points;
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            dirichletValues = zeros(length(points),1);
            dirichletDerivatives = zeros(length(points),1);
            for idx = dirichletNodes
                coord = points(idx,:);
                dirichletValues(idx) = obj.dirichletFunction(coord);
                dirichletDerivatives(idx) = obj.derivedDirichletFunction(coord);
            end
            D = -A*dirichletValues + M*dirichletDerivatives;
        end
        
        function obj = removeDirichletBoundary(obj)
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            obj.stiffnessMatrix(dirichletNodes,:)   = [];
            obj.stiffnessMatrix(:,dirichletNodes)   = [];
            obj.odeMatrix(dirichletNodes,:)         = [];
            obj.odeMatrix(:,dirichletNodes)         = [];
            obj.neumannVector(dirichletNodes)       = [];
            obj.dirichletVector(dirichletNodes)     = [];
        end
        
        function u = reinsertDirichletBoundary(obj, uNoBoundary)
            points = obj.mesh.points;
            u = zeros(length(points),1);
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            nonDirichletNodes = 1:length(points);
            nonDirichletNodes(dirichletNodes) = [];
            u(nonDirichletNodes)=uNoBoundary;
            for node = dirichletNodes
                u(node) = obj.dirichletFunction(points(node,:));
            end
        end
    end
end

