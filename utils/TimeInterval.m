classdef TimeInterval < handle
% TimeSolver(t_0,t_max,h) generates the interval [t_0,t_max] with stepwidth
% h.
%
% t_0:      Left boundary of interval
% t_max:    Right boundary of interval
% h:        Descrete stepwidth
% n_to_plot:
%           Helper variable that safes after how many timeSteps a solution
%           should be saved (actually does not belong here).
% descreteInterval:
%           A vector, that includes each timestep, including t_0 and t_1
%
% See also TimeSolver
    properties (Access=private,Constant)
        eps=10^(-9);
    end
    properties
        t_0
        t_max
        h
        n_to_plot
        descreteInterval
    end
    
    methods
        function obj=TimeInterval(t_0,t_max,h,n_to_plot)
        %Strandart constructor 
            obj.t_0=t_0;
            obj.t_max=t_max;
            obj.h=h;
            obj.n_to_plot = n_to_plot;
            obj.setDescreteInterval();
        end
        
        function N = getNumberOfSteps(obj)
        %Get the number of time steps or with other words, the number of
        %subintervals.
            N=ceil(abs(obj.t_max-obj.t_0)/obj.h-TimeInterval.eps);
        end
        
        function setDescreteInterval(obj)
        %Set the descreteInterval, should be called after changing h, t_min
        %or t_max
            obj.descreteInterval=obj.t_0:obj.h:obj.t_max;
            if(not(equalUpTo(obj.descreteInterval(end),obj.t_max,TimeInterval.eps)))
                obj.descreteInterval=[obj.descreteInterval,obj.t_max];
            end
        end
        
        function setT_0(obj, t_0)
        %Sets t_0 and updates descreteInterval
            obj.t_0=t_0;
            obj.setDescreteInterval();
        end
        
        function setT_max(obj, t_max)
        %Sets t_max and updates descreteInterval
            obj.t_max=t_max;
            obj.setDescreteInterval();
        end
        
        function setH(obj, h)
        %Sets h and updates descreteInterval
            obj.h=h;
            obj.setDescreteInterval();
        end
        
        function setStepWidthFromNumberOfSteps(obj,N)
        %Sets h so that the number of subintervals is N
            obj.setH(abs(obj.t_max-obj.t_0)/N);
        end
        function extendIntervalToEqualStepWidth(obj)
        %In case the length of the interval is not devisable by h, the last
        %subinterval is shorter then the others. This function sets the
        %length of the last subinterval to h by eventually extending the
        %length of the interval.
            obj.setT_max(obj.t_0+obj.getNumberOfSteps()*obj.h);
        end
        
        function equalizeStepWidth(obj)
        %In case the length of the interval is not devisable by h, the last
        %subinterval is shorter then the others. This function changes h,
        %so that each interval has same length.
            obj.setH((obj.t_max-obj.t_0)/obj.getNumberOfSteps());
        end
    end
end