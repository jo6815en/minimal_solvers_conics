function [ll1,ll2] = linesfromquadric(CC,VP,P)

if numel(CC)==16
    CC = CC([1 2 3 4 6 7 8  11 12 16]).';
end

nc = size(VP,2);
if size(VP,1)==3
    VP = [VP;zeros(1,nc)];
end

cl = class(CC);

vp = P*VP;

C1 = CC(1,:);
C2 = CC(2,:);
C3 = CC(3,:);
C4 = CC(4,:);
C5 = CC(5,:);
C6 = CC(6,:);
C7 = CC(7,:);
C8 = CC(8,:);
C9 = CC(9,:);
C10 = CC(10,:);

p1 = P(1);
p2 = P(2);
p3 = P(3);
p4 = P(4);
p5 = P(5);
p6 = P(6);
p7 = P(7);
p8 = P(8);
p9 = P(9);
p10 = P(10);
p11 = P(11);
p12 = P(12);

cc = zeros(6,nc,cl);
cc(1,:) = p1*(C1*p1 + C2*p4 + C3*p7 + C4*p10) + p4*(C2*p1 + C5*p4 + C6*p7 + C7*p10) + p7*(C3*p1 + C6*p4 + C8*p7 + C9*p10) + p10*(C4*p1 + C7*p4 + C9*p7 + C10*p10);
cc(2,:) = p1*(C1*p2 + C2*p5 + C3*p8 + C4*p11) + p4*(C2*p2 + C5*p5 + C6*p8 + C7*p11) + p7*(C3*p2 + C6*p5 + C8*p8 + C9*p11) + p10*(C4*p2 + C7*p5 + C9*p8 + C10*p11);
cc(3,:) = p1*(C1*p3 + C2*p6 + C3*p9 + C4*p12) + p4*(C2*p3 + C5*p6 + C6*p9 + C7*p12) + p7*(C3*p3 + C6*p6 + C8*p9 + C9*p12) + p10*(C4*p3 + C7*p6 + C9*p9 + C10*p12);
cc(4,:) = p2*(C1*p2 + C2*p5 + C3*p8 + C4*p11) + p5*(C2*p2 + C5*p5 + C6*p8 + C7*p11) + p8*(C3*p2 + C6*p5 + C8*p8 + C9*p11) + p11*(C4*p2 + C7*p5 + C9*p8 + C10*p11);
cc(5,:) = p2*(C1*p3 + C2*p6 + C3*p9 + C4*p12) + p5*(C2*p3 + C5*p6 + C6*p9 + C7*p12) + p8*(C3*p3 + C6*p6 + C8*p9 + C9*p12) + p11*(C4*p3 + C7*p6 + C9*p9 + C10*p12);
cc(6,:) = p3*(C1*p3 + C2*p6 + C3*p9 + C4*p12) + p6*(C2*p3 + C5*p6 + C6*p9 + C7*p12) + p9*(C3*p3 + C6*p6 + C8*p9 + C9*p12) + p12*(C4*p3 + C7*p6 + C9*p9 + C10*p12);

[ll1,ll2] = linesfromconic(cc,vp);


