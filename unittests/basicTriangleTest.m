addpath('../utils/')
tri = BasicTriangle(3);
M = [1,0,0;0,1,0;0,0,1;0,0,0];
D = [1,0,0;0,1,0;0,0,1;-1,-1,-1]';
cornerPoints = tri.getCornerPoints();
partialDerivative = tri.getPDEOfbasisFunction();
success = true;
if not(isequal(cornerPoints, M))
  disp('The corner poits are not read correctly!')
  success = false;
end
if not(isequal(partialDerivative, D))
  disp('The partial derivatives are not builded correctly!')
  success = false;
end
for i=1:4
    phi = tri.getBasisFunction(i);
    for j=1:4
        cp = cornerPoints(j,:);
        if not(isequal(phi(cp),kronecker(i,j)))
            disp('The Basis Functions are not working correcty!')
            success = false;
        end
    end
end

if success
    disp([pad('[unittest/basicTriangleTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/basicTriangleTest]',40), ' failed!'])
end


function d=kronecker(i,j)
    if i==j
        d = 1;
        return
    end
    d = 0;
    return
end
