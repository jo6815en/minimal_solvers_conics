function slask=rita3(bild,st)
if nargin == 1,
 st='-';
end;
plot3(bild(1,:),bild(2,:),bild(3,:),st);
slask=[];
