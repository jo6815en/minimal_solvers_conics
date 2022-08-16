function [ll1,ll2] = linesfromconic(cc,vp)

if numel(cc) == 9
   cc = cc([1 2 3 5 6 9]).';
end

nc = size(vp,2);
if size(vp,1) == 2
    vp = [vp;ones(1,nc)];
end

c1 = cc(1,:);
c2 = cc(2,:);
c3 = cc(3,:);
c4 = cc(4,:);
c5 = cc(5,:);
c6 = cc(6,:);

vpx = vp(1,:);
vpy = vp(2,:);
vpz = vp(3,:);

C2 = c4.*vpx.^2 - 2*c2.*vpx.*vpy + c1.*vpy.^2;
C1 = 2*c3.*vpy.^2 - 2*c5.*vpx.*vpy + 2*c4.*vpx.*vpz - 2*c2.*vpy.*vpz;
C0 = c6.*vpy.^2 - 2*c5.*vpy.*vpz + c4.*vpz.^2;

lx1 = -(C1 - sqrt(C1.^2 - 4*C0.*C2))./(2*C2);
lx2 = -(C1 + sqrt(C1.^2 - 4*C0.*C2))./(2*C2);

ly1 = -(vpz + lx1.*vpx)./vpy;
ly2 = -(vpz + lx2.*vpx)./vpy;

ll1 = [lx1;ly1;ones(1,nc)];
ll2 = [lx2;ly2;ones(1,nc)];

