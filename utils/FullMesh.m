classdef FullMesh
    
  properties
    neumannBoundaryFaces
    dirichletBoundaryNodes
    allBoundaryNodes
    triangulation
    points
    rod
    neumannFaceRod
  end
  
  methods
    %tri_in is a list of triangles, where the points are given by the index
    % of the point in the list of points p_in
    function obj = FullMesh(tri_in,p_in,rod,neumannIdentifier)
      obj.triangulation = tri_in;
      obj.points = p_in;
      obj.allBoundaryNodes = obj.getBoundaryNodes();
      obj.dirichletBoundaryNodes = obj.dirichletNodes(neumannIdentifier);
      obj.neumannBoundaryFaces = obj. neumannFaces(neumannIdentifier);
      obj.rod = rod;
      obj.neumannFaceRod = obj.getNeumannFaceRod();
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
    
    function nfr = getNeumannFaceRod(obj)
       nf = obj.neumannBoundaryFaces;
       if isempty(nf)
           nfr = [];
           return
       end
       nfr = zeros(length(nf),1);
       for i = 1:length(nf(:,1))
           for j = 1:length(obj.triangulation)
                if(all(ismember(nf(i,:),obj.triangulation(j,:))))
                    nfr(i) = obj.rod(j);
                    break
                end
           end
       end
    end
    
    function bn = getBoundaryNodes(obj)
      bn = [];
      be = obj.boundaryEdges();
      for edge=be'
        for v=edge'
          if not(ismember(v,bn))
            bn = [bn,v];
          end
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