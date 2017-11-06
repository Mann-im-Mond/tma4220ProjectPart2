classdef QuadratFunctionGenerator
    
  properties
  end
  
  methods
      
    function obj = QuadratFunctionGenerator()
    end
    
    function quad = getQuadratureFunction(obj,dim)
      switch dim
        case 1
          quad = @(x,g) obj.quadrature1D(x(1,:),x(2,:),4,g);
        case 2
          quad = @(x,g) obj.quadrature2D(x(1,:),x(2,:),x(3,:),4,g);
        case 3
          quad = @(x,g) obj.quadrature3D(x(1,:),x(2,:),x(3,:),x(4,:),5,g);
      end
    end
  end
  methods(Static)
    function I = quadrature1D(a, b, Nq, g)
      length = norm((b - a));
      midpoint = (b + a) / 2;
      halfinterval = (b - a) / 2;
      transformed_g = @(x) g(x .* halfinterval + midpoint) * length/2;
      switch Nq
        case 1
          points = 0;
          density = 2;
        case 2
          points = [-sqrt(1 / 3), sqrt(1 / 3)];
          density = [1, 1];
        case 3
          points = [-sqrt(3/5), 0, sqrt(3/5)];
          density = [5/9, 8/9, 5/9];
        case 4
          points = [-sqrt((3 + 2*sqrt(6/5))/(7)), -sqrt((3 - 2*sqrt(6/5))/(7)), ...
                    sqrt((3 - 2*sqrt(6/5))/(7)), sqrt((3 + 2*sqrt(6/5))/(7))];
          density = [(18 - sqrt(30))/(36), (18 + sqrt(30))/(36), ...
                    (18 + sqrt(30))/(36), (18 - sqrt(30))/(36)];
      end
      evaluated = zeros(Nq,1);
      for i = 1:Nq
        evaluated(i) = transformed_g(points(i));
      end
        I = dot(evaluated, density);
    end
    
    function I = quadrature2D(p1, p2, p3, Nq, g)
        switch Nq
            case 1
              barycentic = [[1/3, 1/3, 1/3]];
              density = [1];
            case 3
              barycentic = [[1/2, 1/2, 0]; ...
                        [1/2, 0, 1/2]; ...
                        [0, 1/2, 1/2]];
              density = [1/3, 1/3, 1/3];
            case 4
              barycentic = [[1/3, 1/3, 1/3]; ...
                        [3/5, 1/5, 1/5]; ...
                        [1/5, 3/5, 1/5]; ...
                        [1/5, 1/5, 3/5]];
              density = [-9/16, 25/48, 25/48, 25/48];
          end
          evaluationPoints = barycentic*[p1;p2;p3];
          evaluated = zeros(1,Nq);
          for i = 1:Nq
            evaluated(i) = g(evaluationPoints(i,:));
          end
          V=abs(det([p2-p1;p3-p1]))/2;
          I = dot(evaluated, density) * V;
    end
    
    function I = quadrature3D(p0, p1, p2, p3, Nq, g)
      switch Nq
        case 1
          barycentic = [[1/4, 1/4, 1/4, 1/4]];
          density = [1];
        case 4
          lambda1 = 1/(3*sqrt(5)-5);
          lambda2 = 1/(5+sqrt(5));
          barycentic = [[lambda1, lambda2, lambda2, lambda2]; ...
                        [lambda2, lambda1, lambda2, lambda2]; ...
                        [lambda2, lambda2, lambda1, lambda2]; ...
                        [lambda2, lambda2, lambda2, lambda1]];
          density = [1/4, 1/4, 1/4, 1/4];
        case 5
          barycentic = [[1/4, 1/4, 1/4, 1/4]; ...
                        [1/2, 1/6, 1/6, 1/6]; ...
                        [1/6, 1/2, 1/6, 1/6]; ...
                        [1/6, 1/6, 1/2, 1/6]; ...
                        [1/6, 1/6, 1/6, 1/2]];
          density = [-4/5, 9/20, 9/20, 9/20, 9/20];
      end
      evaluationPoints = barycentic*[p0;p1;p2;p3];
      evaluated = zeros(1,Nq);
      for i = 1:Nq
        evaluated(i) = g(evaluationPoints(i,:));
      end
      V=abs(det([p1-p0;p2-p0;p3-p0]))/6;
      I = dot(evaluated, density) * V;
    end
  end
end