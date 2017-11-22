classdef Plain
%PLAIN(R,v_0) creates the plain, that is rotated by R and
%shifted by v_0 from the plain, where x=0.

% R:            Rotation Matrix, that rotates the plain
% R_inverse:    Inverse Matrix of R, so the matrix that rotates back to the
%               standart plain where x=0
% v_0:          Vector, that moves the plain, so the plain does not need to
%               include the origin [0,0,0].
% normalVector: The normal vector to the plain.
% plain:        Function handle, that applies R and v_0 to a vector in the
%               standart plain x=0
% plain_inverse:Function handle, that applies -v_0 and R_inverse on a
%               vector in R^3
    
    properties
        R
        R_inverse
        v_0
        normalVector
        plain
        plain_inverse
    end
    
    methods
        function obj=Plain(R,v_0)
            %Standard constructor
            obj.R=R;
            obj.R_inverse=R^(-1);
            obj.v_0=v_0;
            obj.normalVector=R*[1;0;0];
            obj.plain=@(y,z) R*[0;y;z]+v_0;
            obj.plain_inverse=@(v) obj.R_inverse*(v-v_0);
        end
        
        function bool=closeTo(obj,v,eps)
        %closeto(v,eps) returns if true if v is at most eps away from obj.
            v_transposed=obj.plain_inverse(v);
            bool=(abs(v_transposed(3))<=eps);
        end
        
        function projection=project(obj,v)
        %project(v) returns the rectengular projection of v into the plain.
            tmp=obj.plain_inverse(v);
            projection=obj.plain(tmp(:,2),tmp(:,3));
        end
    end
    
end

