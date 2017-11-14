classdef BasicTriangle
    
  properties (Access = private)
    dim
  end
  
  methods
    function obj = BasicTriangle(dim_in)
     obj.dim = dim_in;
    end

    %returns the points marking the corners of the basic triangle.
    % TODO this might be optimized by saving them for each object.
    function tri = getCornerPoints(obj)
      tri = zeros(obj.dim+1,obj.dim);
      tri(1:obj.dim,:) = eye(obj.dim);
    end

    %basis function belonging to the ith corner point.
    % TODO this might be optimized by saving them for each object.
    function phi = getBasisFunction(obj,i)
      if (i<obj.dim +1)
        phi = @(x) x(i);
      elseif (obj.dim+1 == i)
        phi = @(x) 1 - norm(x,1);
      else
        error('');
      end
    end

    %returns the Matrix containing all the partial derivatives of the basis
    %functions.
    % TODO this might be optimized by saving them for each object.
    function B = getPDEOfbasisFunction(obj)
      B = (obj.getCornerPoints())';
      B(:,obj.dim+1) = -ones(obj.dim,1);
    end
    
    %returns a matrix, where the entry (i,j) contains:
    %integral (phi_i * phi_j)
    function I = productBasisFunctionIntegralMatrix(obj)
       I = zeros(obj.dim +1);
       for i = 1:(obj.dim+1)
           for j = 1:(obj.dim+1)
               I(i,j) = obj.productBasisFunctionIntegral(i,j);
           end
       end
    end
    
    
    % TODO this can propablz be simplified, as the case d+1 seems to be not
    % necessary, but I lack the theory for that yet.
    
    % returns integral (phi_i * phi_j)
    function I = productBasisFunctionIntegral(obj,i,j)
        d = obj.dim;
        if d<=0
            error('invalid dimension');
        end
        if(or(i==d+1,j==d+1))
            if(i~=j)
                I = (1/d-1/(d+1))*1/(factorial(d-1)) ...
                    - obj.productBasisFunctionIntegral(1,1);
                if(d>1)
                    I = I - (d-1)*obj.productBasisFunctionIntegral(1,2);
                end
                return
            else
                I = 1/(factorial(d)) - 2*d*(1/d-1/(d+1))*1/(factorial(d-1)) ...
                    + d*obj.productBasisFunctionIntegral(1,1);
                if(d>1)
                   I = I + d*(d-1)*obj.productBasisFunctionIntegral(1,2); 
                end
                return
            end
        else
            if(i~=j)
                I = (1/(d+1)-1/(d+2))*(1/(d-1)-1/d)*1/(factorial(d-2));
                return
            else
                I = (1/(d+2)+1/d-2/(d+1))*1/(factorial(d-1));
                return
            end
        end
    end
  end
end