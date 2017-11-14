function b = equalUpTo(x,y, eps)
    b = all((abs(x-y) < eps));
end