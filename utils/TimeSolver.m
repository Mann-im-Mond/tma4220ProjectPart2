classdef TimeSolver < handle
% TimeSolver(u_0,interval,M,A,V) generates the system of ODEs Mu'=Au+V with
% initial values u_0 on the interval. The system can then be solved by the
% method 'solve'.
%
% u_0:      N-dimensional Vektor that is the initial function (u(t_0)=u_0)
% interval: TimeInterval Object, that includes the descretization.
% M:        Invertible NxN Matrix
% A:        NxN Matrix
% V:        N-dimensional Vektor or N-dimensional Functionhandle depending
%           on t
% See also TimeInterval

    properties
        u_0
        interval
        M
        A
        V
    end
    
    properties (Access=private)
        rightSideVector
        iteratorLeftHandSide
        iteratorRightHandSide
    end
    
    methods
        function obj=TimeSolver(u_0,interval,M,A,V)
            obj.u_0=u_0;
            obj.interval=interval;
            obj.M=M;
            obj.A=A;
            obj.V=V;
        end
            
        function u = solve(obj,varargin)
            %Extend the time interval, so it consits only of time steps ...
            %with length h
            obj.interval.extendIntervalToEqualStepWidth();
            %Check whether V is a vector or a function handle
            if(isa(obj.V,'function_handle'))
                VisFunctionHandle=true;
            elseif(isa(obj.V,'double'))
                VisFunctionHandle=false;
            else
                error('V has to be either a function handle or a vector of doubles');
            end
            
            method='';
            stoppingCondition=@(u) false;
            stoppingConditionIsSet= false;
            saveSolutionEveryTimeSteps=obj.interval.getNumberOfSteps();
            saveSolutionEveryTimeStepsIsSet= false;
            skipNextLoop=false;
            for j=1:length(varargin)
                if(skipNextLoop)
                    skipNextLoop=false;
                    continue
                end
                switch varargin{j}
                    case 'method'
                        if j<length(varargin)
                            method=varargin{j+1};
                        end
                        skipNextLoop=true;
                    case 'stoppingCondition'
                        if j<length(varargin)
                            if(isa(varargin{j+1},'function_handle'))
                                stoppingCondition=varargin{j+1};
                                stoppingConditionIsSet=true;
                            end
                        end
                        skipNextLoop=true;
                    case 'safeSolutionEvery'
                        if j<length(varargin)
                            if(isa(varargin{j+1},'double'))
                                saveSolutionEveryTimeSteps=varargin{j+1};
                                saveSolutionEveryTimeStepsIsSet=true;
                            end
                        end
                        skipNextLoop=true;
                end
            end
            switch method
                case {'euler','forwardEuler','Euler','forward euler','forward Euler','Forward Euler'}
                    if(VisFunctionHandle)
                        obj.rightSideVector=@(t) obj.interval.h.*obj.V(t);
                        obj.iteratorRightHandSide=@(v,t) v+obj.interval.h.*obj.A*v+obj.rightSideVector(t);
                    else
                        obj.rightSideVector=obj.interval.h.*obj.V;
                        obj.iteratorRightHandSide=@(v,t) v+obj.interval.h.*obj.A*v+obj.rightSideVector;
                    end                        
                    u=obj.explicitMethod(stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet);
                case {'improvedEuler','improved euler','improved Euler','Improved Euler'}%Not yet finnished
                    if(VisFunctionHandle)
                        obj.rightSideVector=@(t) (obj.interval.h/2).*(obj.V(t)+obj.V(t+obj.interval.h));
                        obj.iteratorRightHandSide=@(v,t) v+obj.interval.h.*obj.A*v+((obj.interval.h)^2/2).*TimeSolver.linearSolver(obj.M,(obj.A*v+obj.V(t)))+obj.rightSideVector(t);
                    else
                        obj.rightSideVector=obj.interval.h.*obj.V;
                        obj.iteratorRightHandSide=@(v,t) v+obj.interval.h.*obj.A*v+obj.rightSideVector;
                    end                        
                    u=obj.explicitMethod(stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet);
                case {'backwardsEuler','backwards euler','backwards Euler','Backwards Euler'}
                    obj.iteratorLeftHandSide=obj.M-obj.interval.h.*obj.A;
                    if(VisFunctionHandle)
                        obj.rightSideVector=@(t) obj.interval.h.*obj.V(t);
                        obj.iteratorRightHandSide=@(v,t) v+obj.rightSideVector(t);
                    else
                        obj.rightSideVector=obj.interval.h.*obj.V;
                        obj.iteratorRightHandSide=@(v,t) v+obj.rightSideVector;
                    end                        
                    u=obj.implicitOneStepMethod(stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet);
                otherwise %Crank-Nicolson Method
                    if(not(ismember(method,{'','crankNicolson','Crank-Nicolson','crank-cicolson','Crank Nicolson','crank nicolson'})))
                        warning(['Method ',method,' not recognized, use Crank-Nicolson as fallback']);
                    end
                    obj.iteratorLeftHandSide=obj.M-(obj.interval.h/2).*obj.A;
                    if(VisFunctionHandle)
                        obj.rightSideVector=@(t) obj.interval.h.*obj.V(t);
                        obj.iteratorRightHandSide=@(v,t) v+(obj.interval.h/2).*obj.A*v+obj.rightSideVector(t);
                    else
                        obj.rightSideVector=obj.interval.h.*obj.V;
                        obj.iteratorRightHandSide=@(v,t) v+(obj.interval.h/2).*(obj.A*v)+obj.rightSideVector;
                    end                        
                    u=obj.implicitOneStepMethod(stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet);
            end
        end
    end
    methods (Access=private,Static)
        function u=linearSolver(A,b)
            u=A\b;
        end
    end
        
    methods (Access=private)
        function [U,t]=explicitMethod(obj,stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet)
            %initiate the interval
            K=obj.interval.getNumberOfSteps();
            %initiate u
            u=zeros(length(obj.u_0),2);
            u(:,1)=obj.u_0;
            swap=1;
            
            if(saveSolutionEveryTimeStepsIsSet)
                %initialize return matrix
                j=0;
                U=zeros(length(u),floor(K/saveSolutionEveryTimeSteps));
                if(stoppingConditionIsSet)
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=obj.iteratorRightHandSide(u(:,3-swap),t);
                        if((mod(k-1,saveSolutionEveryTimeSteps))==0)
                            j=j+1;
                            U(:,j)=obj.M\u(:,swap);
                        end
                        if(stoppingCondition(u(:,swap)))
                            if(j==0)
                                U=u(:,swap);
                            else
                                U=U(:,j);
                            end
                            return;
                        end
                    end
                else
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=obj.iteratorRightHandSide(u(:,3-swap),t);
                        if((mod(k-1,saveSolutionEveryTimeSteps))==0)
                            j=j+1;
                            U(:,j)=obj.M\u(:,swap);
                        end
                    end
                end
            else
                if(stoppingConditionIsSet)
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=obj.iteratorRightHandSide(u(:,3-swap),t);
                        if(stoppingCondition(u(:,swap)))
                            U=obj.M\u(:,swap);
                            return;
                        end
                    end
                else
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=obj.iteratorRightHandSide(u(:,3-swap),t);
                    end
                end
                U=obj.M\u(:,swap);
            end
        end
        
        function [U,t]=implicitOneStepMethod(obj,stoppingConditionIsSet,stoppingCondition,saveSolutionEveryTimeSteps,saveSolutionEveryTimeStepsIsSet)
            %initiate the interval
            K=obj.interval.getNumberOfSteps();
            %initiate u
            u=zeros(length(obj.u_0),2);
            u(:,1)=obj.u_0;
            swap=1;
            
            if(saveSolutionEveryTimeStepsIsSet)
                %initialize return matrix
                j=0;
                U=zeros(length(u),floor(K/saveSolutionEveryTimeSteps));
                if(stoppingConditionIsSet)
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=TimeSolver.linearSolver(obj.iteratorLeftHandSide,obj.iteratorRightHandSide(u(:,3-swap),t));
                        if((mod(k-1,saveSolutionEveryTimeSteps))==0)
                            j=j+1;
                            U(:,j)=u(:,swap);
                        end
                        if(stoppingCondition(u))
                            U=U(:,j);
                            return;
                        end
                    end
                else
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=TimeSolver.linearSolver(obj.iteratorLeftHandSide,obj.iteratorRightHandSide(u(:,3-swap),t));
                        if((mod(k-1,saveSolutionEveryTimeSteps))==0)
                            j=j+1;
                            U(:,j)=u(:,swap);
                        end
                    end
                end
            else
                if(stoppingConditionIsSet)
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=TimeSolver.linearSolver(obj.iteratorLeftHandSide,obj.iteratorRightHandSide(u(:,3-swap),t));
                        if(stoppingCondition(u))
                            U=u(:,swap);
                            return;
                        end
                    end
                else
                    for k = 2:K
                        swap=3-swap;
                        t=obj.interval.descreteInterval(k);
                        u(:,swap)=TimeSolver.linearSolver(obj.iteratorLeftHandSide,obj.iteratorRightHandSide(u(:,3-swap),t));
                    end
                end
                U=u(:,swap);
            end
        end
    end
end