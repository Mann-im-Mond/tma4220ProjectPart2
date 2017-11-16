classdef TimeInterval < handle
% TimeSolver(t_0,t_max,h) generates the interval [t_0,t_max] with stepwidth
% h.
%
% See also 
    properties (Access=private,Constant)
        eps=10^(-9);
    end
    properties
        t_0
        t_max
        h
        descreteInterval
    end
    
    methods
        function obj=TimeInterval(t_0,t_max,h)
            obj.t_0=t_0;
            obj.t_max=t_max;
            obj.h=h;
            obj.setDescreteInterval();
        end
        function N = getNumberOfSteps(obj)
            N=ceil(abs(obj.t_max-obj.t_0)/obj.h-TimeInterval.eps);
        end
        function setDescreteInterval(obj)
            obj.descreteInterval=obj.t_0:obj.h:obj.t_max;
            if(not(equalUpTo(obj.descreteInterval(end),obj.t_max,TimeInterval.eps)))
                obj.descreteInterval=[obj.descreteInterval,obj.t_max];
            end
        end
        function setStepWidthFromNumberOfSteps(obj,N)
            obj.h=abs(obj.t_max-obj.t_0)/N;
            obj.setDescreteInterval();
        end
        function extendIntervalToEqualStepWidth(obj)
            obj.t_max=obj.t_0+obj.getNumberOfSteps()*obj.h;
            obj.setDescreteInterval();
        end
        function equalizeStepWidth(obj)
            obj.h=(obj.t_max-obj.t_0)/obj.getNumberOfSteps();
            obj.setDescreteInterval();
        end
    end
end