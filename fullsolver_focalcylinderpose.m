function [R,f] = fullsolver_focalcylinderpose(l,vpd)
% function [R,f] = fullsolver_focalcylinderpose(l,vpd)
%
% input: four image lines l:3x4; four 3D-line directions vpd: 3x4
% output: Rotations 9xnsol; focal lengths f: 1xnsol

nl = 4;
lille = 1e-9;
boffa = zeros(nl,9);
for iii = 1:nl
    di = l(:,iii)*vpd(:,iii)';
    boffa(iii,:)=di(:)';
end
dd = null(boffa);

sols = solver_focalcylinderpose(dd(:));
sols(:,sum(abs(imag(sols)))>lille)=[];

nsol = size(sols,2);
R = zeros(9,nsol);
f = zeros(1,nsol);

badids = [];
for iii = 1:nsol
    A = reshape(dd*[1;sols(:,iii)],3,3);
    A = A/norm(A(3,:));
    n1 = norm(A(1,:));
    n2 = norm(A(2,:));
    if ((n1-n2)/(n1+n2))<.1
        f(iii) = 0.5*(norm(A(1,:))+norm(A(2,:)));
        Ri = A;
        Ri(1:2,:)=Ri(1:2,:)/f(iii);
        if det(Ri)<0
            Ri = -Ri;
        end
        R(:,iii) = Ri(:);
    else
        badids = [badids iii];
    end
end
R(:,badids)=[];
f(badids)=[];

