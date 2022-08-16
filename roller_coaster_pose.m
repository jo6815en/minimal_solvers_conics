
im = imread('images/roller_coaster.jpg');
load('images/roller_coaster_data.mat');
nrcyl = length(CCs);
Pgtn = K\Pgt;


%% reprojections using gt cameras
figure(1)
clf
imshow(im)
hold on

for jjj = 1:nrcyl
    [l1,l2] = linesfromquadric(CCs(:,:,jjj),VPs(:,jjj),Pgt);
    rital(ll(:,2*jjj-1),'b-');
    rital(l1,'r--');
    rital(ll(:,2*jjj),'b-');
    rital(l2,'r--');
end
legend({'Detected lines','Reprojected lines'})
title('Reprojections using ground truth cameras')



%% test pose with focal

fe = 1000; % scale problem for numerical stability
Ke = K;
Ke(1) = fe;
Ke(5) = fe;

bnd = 0.01;
bnd2 = 0.005;
lille = 1e-9;


lln = Ke'*ll;
lln = lln./sqrt(sum(lln.^2));


vpds = [VPs(1:3,:);VPs(1:3,:)];
vpds = reshape(vpds,3,[]);

CDs = cell(1,2*nrcyl);
for iii = 1:nrcyl
    CDs{2*iii-1} = CCs(:,:,iii)/norm(CCs(:,:,iii));
    CDs{2*iii} = CCs(:,:,iii)/norm(CCs(:,:,iii));
end

warning('off');


ins = 0;
fbest = nan;

% test all quads of lines
cids = nchoosek(1:2*nrcyl,4);
for ci = 1:size(cids,1)
    if ~any(diff(cids(ci,:)')==1)
        
        
        l = lln(:,cids(ci,:));
        vpd = vpds(:,cids(ci,:));
        [R,f] = fullsolver_focalcylinderpose(l,vpd);
        for idi = 1:length(f)
            ins = [ins;sum(abs(diag(lln'*diag([f(idi) f(idi) 1])*reshape(R(:,idi),3,3)*vpds))<bnd)];
            if ins(end)>max(ins(1:end-1))
                Rbest = reshape(R(:,idi),3,3);
                fbest = f(idi);
                idsbest = cids(ci,:);
            end
            
        end
    end
end

ins2 = 0;
% check translation
Pe = Rbest;
K2 = diag([fbest fbest 1]);
lln2 =   K2'*lln;

alldata = zeros(10,2*nrcyl);
for jjj = 1:2*nrcyl
    alldata(:,jjj) = getdata_translation_cylinderpose(CDs{jjj},lln2(:,jjj),Pe)';
end

qbest = rot2quat(Rbest);
data = zeros(40,1);
data(10:13) = qbest;
% test all triplets of lines
cids = nchoosek(1:2*nrcyl,3);
for ci = 1:size(cids,1)
    
    data(1:9) = lln2(:,cids(ci,:));
    data(14:22) = -CDs{cids(ci,1)}([1 2 3 4 6 7 8 11 12])/CDs{cids(ci,1)}(16);
    data(23:31) = -CDs{cids(ci,2)}([1 2 3 4 6 7 8 11 12])/CDs{cids(ci,2)}(16);
    data(32:40) = -CDs{cids(ci,3)}([1 2 3 4 6 7 8 11 12])/CDs{cids(ci,3)}(16);
    
    %data = alldata(:,cids(ci,:));
    
    
    %sols = solver_translation_cylinderpose(data(:));
    sols = solver_conics_translation(data);
    sols = sols(:,sum(abs(imag(sols)))<lille);
    nsols = size(sols,2);
    for jjj = 1:nsols
        t1 = sols(1,jjj);
        t2 = sols(2,jjj);
        t3 = sols(3,jjj);
        vv = [t1^2, t1*t2, t1*t3, t1, t2^2, t2*t3, t2, t3^2, t3, 1].';
        res = alldata'*vv;
        insi = sum(abs(res)<bnd2);
        ins2 = [ins2;insi];
        if ins2(end)>max(ins2(1:end-1))
            tbest = [t1;t2;t3];
        end
    end
end
bestP = [Rbest tbest];
bestf = fbest*fe;

disp(Pgtn)
disp(bestP)

disp(K(1))
disp(bestf)


%% reprojections using estimated cameras
figure(2)
clf
Ke = K;
Ke(1) = bestf;
Ke(5) = bestf;

Pe = Ke*bestP;
imshow(im)
hold on


for jjj = 1:nrcyl
    [l1,l2] = linesfromquadric(CCs(:,:,jjj),VPs(:,jjj),Pe);
    rital(ll(:,2*jjj-1),'b-');
    rital(l1,'r--');
    rital(ll(:,2*jjj),'b-');
    rital(l2,'r--');
end
legend({'Detected lines','Reprojected lines'})
title('Reprojections using estimated cameras')





