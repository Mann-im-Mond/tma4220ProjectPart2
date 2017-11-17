classdef HeatEquationSolver < handle
   
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
            syms x t;
            h = matlabFunction(diff(dirichletFunction, t),'vars', [x t]);
            reset(symengine);
            obj.derivedDirichletFunction = h;
            obj.neumannFunction = neumannFunction;
        end
        
        function u = solve(obj)
            % TODO this has still to be implemented!
        end 
    end
    
    methods
        function  initializeSystem(obj)
            obj.stiffnessMatrix = obj.getStiffnessMatrix();
            obj.odeMatrix = obj.getODEMatrix();
            syms x t;
            if(diff(obj.neumannFunction, t) == 0)
                obj.neumannVector = obj.getNeumannVector(@(x) obj.neumannFunction(x,0));
            else
                obj.neumannVector = @(t) obj.getNeumannVector(@(x) obj.neumannFunction(x,t));
            end
            if(diff(obj.dirichletFunction, t) == 0)
                obj.dirichletVector = obj.dirichletVector(@(x) obj.dirichletFunction(x,0), ...
                    @(x) obj.derivedDirichletFunction(x,0));
            else
                obj.dirichletVector = @(t) obj.dirichletVector(@(x) obj.dirichletFunction(x,t), ...
                    @(x) obj.derivedDirichletFunction(x,t));
            end
            reset(symengine);
            obj.removeDirichletBoundaryFromMatrices();
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
        
        function N = getNeumannVector(obj, timeConstantNeumann)
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
                transformedNeumann = @(x) timeConstantNeumann(trafo(x));
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
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            N(dirichletNodes)       = [];
        end
        
        function D = getDirichletVector(obj, timeConstantDirichlet, timeConstantDerivative)
            A = obj.getStiffnessMatrix();
            M = obj.getODEMatrix();
            points = obj.mesh.points;
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            dirichletValues = zeros(length(points),1);
            dirichletDerivatives = zeros(length(points),1);
            for idx = dirichletNodes
                coord = points(idx,:);
                dirichletValues(idx) = timeConstantDirichlet(coord');
                dirichletDerivatives(idx) = timeConstantDerivative(coord');
            end
            D = -A*dirichletValues + M*dirichletDerivatives;
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            D(dirichletNodes)     = [];
        end
        
        function removeDirichletBoundaryFromMatrices(obj)
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            obj.stiffnessMatrix(dirichletNodes,:)   = [];
            obj.stiffnessMatrix(:,dirichletNodes)   = [];
            obj.odeMatrix(dirichletNodes,:)         = [];
            obj.odeMatrix(:,dirichletNodes)         = [];
        end
        
        function u = reinsertDirichletBoundary(obj, uNoBoundary)
            points = obj.mesh.points;
            u = zeros(length(points),1);
            dirichletNodes = obj.mesh.dirichletBoundaryNodes;
            nonDirichletNodes = 1:length(points);
            nonDirichletNodes(dirichletNodes) = [];
            u(nonDirichletNodes)=uNoBoundary;
            for node = dirichletNodes
                u(node) = obj.dirichletFunction(points(node,:),0);
            end
        end
    end
end

