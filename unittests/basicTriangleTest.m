addpath('../utils/')
addpath('../helperFunctions')
dim = 3;
tri = BasicTriangle(dim);
M = [1,0,0;0,1,0;0,0,1;0,0,0];
D = [1,0,0;0,1,0;0,0,1;-1,-1,-1]';
cornerPoints = tri.getCornerPoints();
partialDerivative = tri.getPDEOfbasisFunction();

quadFunc = QuadratFunctionGenerator().getQuadratureFunction(dim);
integrals = tri.productBasisFunctionIntegralMatrix();
realIntegrals = zeros(dim+1);

for i=1:dim+1
    for j=1:dim+1
        phi_i = tri.getBasisFunction(i);
        phi_j = tri.getBasisFunction(j);
        eval = @(x) phi_i(x).*phi_j(x);
        realIntegrals(i,j)=quadFunc(cornerPoints,eval);
    end
end

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

if not(equalUpTo(integrals, realIntegrals, 1e-10))
  disp('The integrals of the basis functions are not computed correctly!')
  success = false;
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
