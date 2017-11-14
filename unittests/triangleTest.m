path(pathdef)
addpath('../utils/')
addpath('../helperFunctions')
tri1 = Triangle([0,0,0;1,1,0;1,1,1]);
tri2 = Triangle([0,0,0;1,1,0;1,1,1;2,0,0]);
J1 = tri1.getJacobiFromBasis();
J2 = tri2.getJacobiFromBasis();
realJ1 = [-1,0;-1,0;-1,-1];
realJ2 = [-2,-1,-1;0,1,1;0,0,1];
trafo1 = tri1.getTrafoFromBasis();
trafo2 = tri2.getTrafoFromBasis();
invTrafo1 = tri1.getTrafoToBasis();
invTrafo2 = tri2.getTrafoToBasis();
btri1 = BasicTriangle(2);
btri2 = BasicTriangle(3);
corners1 = tri1.cornerPoints;
corners2 = tri2.cornerPoints;
bcorners1 = btri1.getCornerPoints();
bcorners2 = btri2.getCornerPoints();
vol1 = tri1.getVolume();
vol2 = tri2.getVolume();
realVol1 = sqrt(2)/2;
realVol2 = 1/3;
success = true;
com1 = tri1.centerOfMass();
com2 = tri2.centerOfMass();
realCom1 = [2/3;2/3;1/3];
realCom2 = [1;1/2;1/4];

%check the Jacobi matrices
if not(isequal(J1,realJ1))
    disp('The Jaco matrices are not correct!')
 	success = false;
end
if not(isequal(J2,realJ2))
    disp('The Jaco matrices are not correct!')
 	success = false;
end

%check the transformations
for i=1:3
    if not(isequal(trafo1(bcorners1(i,:)),corners1(i,:)))
        disp('The trafo function is not working correct!')
        success = false;
    end
end
for i=1:4
    if not(isequal(trafo2(bcorners2(i,:)),corners2(i,:)))
        disp('The trafo function is not working correct!')
        success = false;
    end
end
%{
for i=1:3
    if or(not(isequal(trafo1(invTrafo1(corners1(i,:))),corners1(i,:))), ...
            not(isequal(invTrafo1(trafo1(bcorners1(i,:))),bcorners1(i,:))))
        disp('The invers trafo function is not working correct!')
        success = false;
    end
end
%}
for i=1:4
    if or(not(isequal(trafo2(invTrafo2(corners2(i,:))),corners2(i,:))), ...
            not(isequal(invTrafo2(trafo2(bcorners2(i,:))),bcorners2(i,:))))
        disp('The invers trafo function is not working correct!')
        success = false;
    end
end

%check the volumes
if not(equalUpTo(vol1,realVol1,10e-6))
    disp('The volumes are not correct!')
 	success = false;
end
if not(equalUpTo(vol2,realVol2,10e-6))
    disp('The volumes are not correct!')
 	success = false;
end

%check the center of mass
if not(equalUpTo(com1,realCom1,1e-6))
    disp('The centers of mass are not correct!')
 	success = false;
end
if not(equalUpTo(com2,realCom2,1e-6))
    disp('The centers of mass are not correct!')
 	success = false;
end

if success
    disp([pad('[unittest/triangleTest]',40), 'succeeded!'])
else
    disp([pad('[unittest/triangleTest]',40), 'failed!'])
end
