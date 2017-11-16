path(pathdef)
addpath('../utils/')
addpath('../helperFunctions')
%h should not be a multiple of t_max-t_0, this case will be tested below.
t_0=0;
t_max=20;
h=6;
timeInterval=TimeInterval(t_0,t_max,h);

success = true;

plot(timeInterval.descreteInterval,'o');
answer=input(['Do you see ',num2str(timeInterval.getNumberOfSteps()+1)',' circles in the plot? (yes/no)'],'s');
while(not(ismember(answer,{'yes','no'})))
    answer=input(['Do you see ',num2str(timeInterval.getNumberOfSteps()+1)',' circles in the plot? (yes/no)'],'s');
end
if(isequal(answer,'no'))
    disp('Number of steps not calculated correctly.');
    success=false;
end
if not(isequal(timeInterval.descreteInterval(1), t_0))
  disp('The left boundary of the interval does not fit.')
  success = false;
end
if not(isequal(timeInterval.descreteInterval(end), t_max))
  disp('The right boundary of the interval does not fit.')
  success = false;
end
timeInterval.extendIntervalToEqualStepWidth();
if not(isequal(timeInterval.descreteInterval(end)-timeInterval.descreteInterval(end-1),h))
  disp('The extension to gain equal step width did not work.')
  success = false;
end
timeInterval.t_max=t_max;
timeInterval.equalizeStepWidth();
h_tmp_1=timeInterval.descreteInterval(2)-timeInterval.descreteInterval(1);
for i=2:length(timeInterval.descreteInterval)-1
    h_tmp=timeInterval.descreteInterval(i+1)-timeInterval.descreteInterval(i);
    if not(equalUpTo(h_tmp_1,h_tmp,10^(-9)))
      disp('The equalization of step width did not work.')
      success = false;
    end
end

timeInterval.t_max=t_max;
numberOfSteps=10;
timeInterval.setStepWidthFromNumberOfSteps(numberOfSteps);
plot(timeInterval.descreteInterval(),'o');
answer=input(['Do you see ',num2str(timeInterval.getNumberOfSteps()+1),' circles in the plot? (yes/no)'],'s');
while(not(ismember(answer,{'yes','no'})))
    answer=input(['Do you see ',num2str(numberOfSteps+1),' circles in the plot? (yes/no)'],'s');
end
if(isequal(answer,'no'))
    disp('Width not correctly set from given number of steps');
    success=false;
end
if not(isequal(timeInterval.descreteInterval(1), t_0))
  disp('The left boundary of the interval does not fit.')
  success = false;
end
if not(isequal(timeInterval.descreteInterval(end), t_max))
  disp('The right boundary of the interval does not fit.')
  success = false;
end
close;

if success
    disp([pad('[unittest/timeIntervalTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/timtIntervalTest]',40), ' failed!'])
end
