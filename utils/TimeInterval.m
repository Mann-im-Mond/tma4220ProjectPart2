classdef TimeInterval
% TimeSolver(t_0,t_max,h) generates the interval [t_0,t_max] with stepwidth
% h.
%
% See also 
    properties (Constant)
        eps=10^(-9);
    end
    properties
        t_0
        t_max
        h
    end
    
    methods
        function N = getNumberOfSteps(obj)
            N=ceil((obj.t_0-obj.t_max)/obj.h-TimeInterval.eps);
        end
        function interval = getDescreteInterval(obj)
            interval=obj.t_0:obj.h:obj.t_max;
            if(not(equalUpTo(interval(end),obj.t_max,TimeInterval.eps)))
                interval=[interval,obj.t_max];
            end
        end
        function obj=setStepwidth(obj,N)
            obj.h=(obj.t_0-obj.t_max)/N;
        end
        function obj=extendIntervalToEqualStepwidth(obj)
            obj.t_max=obj.t_0+obj.getNumberOfSteps()*obj.h;
        end
    end
end