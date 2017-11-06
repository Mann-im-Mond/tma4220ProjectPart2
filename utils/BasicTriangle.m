classdef BasicTriangle
  properties (Access = private)
    dim
  end
  methods
    function obj = BasicTriangle(dim_in)
     obj.dim = dim_in;
    end

    %returns the points marking the corners of the basic triangle.
    % TODO this might be optimized by saving them in the class.
    function tri = getCornerPoints(obj)
      tri = zeros(obj.dim+1,obj.dim);
      tri(1:obj.dim,:) = eye(obj.dim);
    end

    %basis function belonging to the ith corner point.
    % TODO this might be optimized by saving them in the class.
    function phi = getBasisFunction(obj,i)
      if (i<obj.dim +1)
        phi = @(x) x(i);
      end
      if (obj.dim+1 == i)
        phi = @(x) 1 - norm(x,1);
      end
    end

    %returns the Matrix containing the PDE's of the basis functions.
    % TODO this might be optimized by saving them in the class.
    function B = getPDEOfbasisFunction(obj)
      B = (obj.getCornerPoints())';
      B(:,obj.dim+1) = -ones(obj.dim,1);
    end
  end
end