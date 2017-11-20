classdef Plain
    %PLAIN(R,v_0) creates the plain, that is rotated by R and
    %shifted by v_0 from the plain, where x=0.
    
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
            obj.R=R;
            obj.R_inverse=R^(-1);
            obj.v_0=v_0;
            obj.normalVector=R*[1;0;0];
            obj.plain=@(y,z) R*[0;y;z]+v_0;
            obj.plain_inverse=@(v) obj.R_inverse*(v-v_0);
        end
        
        function bool=closeTo(obj,v,eps)
        %closeto(v) returns if true if v is at most eps away from obj.
            v_transposed=obj.plain_inverse(v);
            bool=(abs(v_transposed(3))<=eps);
        end
    end
    
end

