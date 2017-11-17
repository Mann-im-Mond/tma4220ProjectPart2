path(pathdef)
addpath('../utils/')
addpath('../helperFunctions')

dim=3;
outputEverySteps=100;
numberOfOutputs=10;
eps=10^(-2);

success = true;

t_0=0;
t_max=1;
h=(t_max-t_0)/(outputEverySteps*numberOfOutputs);
timeInterval=TimeInterval(t_0,t_max,h,outputEverySteps);
u_0=zeros(dim,1);
M=sparse(eye(dim));
A=M;
V=@(t) exp(t.*ones(dim,1));
timeSolver=TimeSolver(u_0,timeInterval,M,A,V);
u_exact=@(t) t.*exp(t);

for method={'Forward Euler','Improved Euler','Backwards Euler','Crank-Nicolson'}
    u=timeSolver.solve('method',method{1},'saveSolutionEvery',outputEverySteps,'stoppingCondition');
    [~,counter]=size(u);
    for i=2:counter
        if(not(equalUpTo((u_exact(((i-1)*outputEverySteps).*timeSolver.interval.h)-u(1,i))/u_exact(((i-1)*outputEverySteps).*timeSolver.interval.h),0,eps)))
          disp([method{1},' failed at t=',num2str((i-1)*outputEverySteps)]);
          success = false;
        end
    end
end
for method={'Forward Euler','Improved Euler','Backwards Euler','Crank-Nicolson'}
    u=timeSolver.solve('method',method{1},'saveSolutionEvery',outputEverySteps,'stoppingCondition',@(u) (norm(u,1)/dim>=exp(1)/2));
    [~,counter]=size(u);
    if(not(counter==7))
        disp('Stopping Condition does not work');
    end
    for i=2:counter
        if(not(equalUpTo((u_exact(((i-1)*outputEverySteps).*timeSolver.interval.h)-u(1,i))/u_exact(((i-1)*outputEverySteps).*timeSolver.interval.h),0,eps)))
          disp([method{1},' failed at t=',num2str((i-1)*outputEverySteps),', with stopping coindition set']);
          success = false;
        end
    end
end
for method={'Forward Euler','Improved Euler','Backwards Euler','Crank-Nicolson'}
    u=timeSolver.solve('method',method{1},'stoppingCondition',@(u) (norm(u,1)/dim>=exp(1)));
    if(not(equalUpTo((u_exact((t_max))-u(1))/u_exact(t_max),0,eps)))
        disp([method{1},' failed with stopping condition set']);
        success = false;
    end
end
for method={'Forward Euler','Improved Euler','Backwards Euler','Crank-Nicolson'}
    u=timeSolver.solve('method',method{1});
    if(not(equalUpTo((u_exact(t_max)-u(1))/u_exact(t_max),0,eps)))
        disp([method{1},' failed']);
        success = false;
    end
end

if success
    disp([pad('[unittest/timeSolverTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/timtSolverTest]',40), ' failed!'])
end
