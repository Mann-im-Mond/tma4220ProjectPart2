classdef Mesh
    
  properties
    neumannBoundaryFaces
    dirichletBoundaryNodes
    triangulation
    points
  end
  
  methods
    %tri_in is a list of triangles, where the points are given by the index
    % of the point in the list of points p_in
    function obj = Mesh(tri_in,p_in,neumannIdentifier)
      obj.triangulation = tri_in;
      obj.points = p_in;
      obj.dirichletBoundaryNodes = obj.dirichletNodes(neumannIdentifier);
      obj.neumannBoundaryFaces = obj. neumannFaces(neumannIdentifier);
    end
    
    %returns the cornerpoints of a triangle (or higher dim equivalent)
    % given by the indices
    function cor = corners(obj,triangle)
      cor = obj.points(triangle,:);
    end

    %returns all edges, which lay on the boundary of the mesh
    function be = boundaryEdges(obj)
      TR = triangulation(obj.triangulation, obj.points);
      be = freeBoundary(TR);
    end

    %returns all neumann faces (they are of dimension dim-1)
    % the neumann identifier is a function, which returns true, iff
    % the input point is in the neumann Boundary and false if not
    % WARNING: neumannIdentifier can do anything with points, that are not
    % on the boundary at all!
    function ne = neumannFaces(obj, neumannIdentifier)
      be = obj.boundaryEdges();
      ne = [];
      for edge=be'
        edge = edge';
        tester = true;
        for v = edge
          p = obj.points(v,:);
          tester = and(tester,neumannIdentifier(p));
        end
        if(tester)
          ne = [ne;edge];
        end
      end
    end
    
    %returns all nodes on the dirichlet boundary
    % the neumann identifier is a function, which returns true, iff
    % the input point is in the neumann Boundary and false if not
    % WARNING: neumannIdentifier can do anything with points, that are not
    % on the boundary at all!
    function dn = dirichletNodes(obj, neumannIdentifier)
      dn = [];
      points = obj.points;
      be = obj.boundaryEdges();
      for edge=be'
        for v=edge'
          if(not(neumannIdentifier(points(v,:))))
            if not(ismember(v,dn))
              dn = [dn,v];
            end
          end
        end
      end
    end
    function d = dim(obj)
      d = length(obj.points(1,:));
    end
  end
end