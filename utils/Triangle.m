classdef Triangle

    properties
        cornerPoints
        dimTriangle
        dimSpace
    end
    
    methods
        function obj = Triangle(corners)
            obj.cornerPoints = corners;
            obj.dimTriangle = length(corners(:,1))-1;
            obj.dimSpace = length(corners(1,:));
        end
        
        %returns the jacobi matrix belonging to the transformation from the
        %basis triangle to this one
        function J = getJacobiFromBasis(obj)
        	J=zeros(obj.dimSpace,obj.dimTriangle);
        	for i=1:obj.dimTriangle
                J(:,i) = (obj.cornerPoints(i,:) ...
                    - obj.cornerPoints(obj.dimTriangle+1,:) )';
        	end
        end
        
        %returns the transformation function from the basis Triangle to
        %this one
        function trafo = getTrafoFromBasis(obj)
            J = obj.getJacobiFromBasis();
            trafo = @(x) (J*(x)' ...
                + obj.cornerPoints(obj.dimTriangle+1,:)')';
        end
        
        %returns the transformation function from this one to basis Triangle
        % TODO implement this for general triangles - by now the matrix is
        % not necessarily quradratic.
        function trafo = getTrafoToBasis(obj)
            J = obj.getJacobiFromBasis();
            trafo = @(x) (inv(J)*((x)' ...
                - obj.cornerPoints(obj.dimTriangle+1,:)'))';
        end    
        
        %returns the volume of this triangle
        function vol = getVolume(obj)
          if obj.dimTriangle==1
            vol = 1;
            return
          end
          J = obj.getJacobiFromBasis();
          O = orth(J);
          edges_in_new_basis = O\J;
          vol = abs(det(edges_in_new_basis))/factorial(obj.dimTriangle);
        end
        
        function com = centerOfMass(obj)
            com = zeros(obj.dimSpace,1);
            for corner = obj.cornerPoints(:,:)'
                com = com + corner;
            end
            com = com / (obj.dimTriangle +1);
        end
    end
    
end

