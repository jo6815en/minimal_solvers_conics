function [Roterr, transerr] = conics_synt_parallell(noise, plots, seed, dbug)
%Parallell trees, synthetic data, possible noise.

% noise - standard deviation. Level of added noise < 0.025
% plots - plot level, 0 - no plots, 1 - resulting plots
% seed - control random generator
% dbug - the dbug mode 0 - off, 1 - on, all prints and plots are shown.
%
% Roterr - rotation error in radians between estimated rotation and ground
% truth
% transerr - translation error (pixels?) between estimated translation and
% gound truth

if nargin < 4;  dbug = 0; end
if nargin > 2;  rng(seed); end
if nargin < 2;  plots = 1; end
if nargin < 1;  noise = 0; end


clearvars -except noise plots seed dbug; close all;
parallell = 1;
bnd = 0.001;

get_paths

% Load synthetic data.
trees = load_cylinders(parallell, plots, dbug);
n_conics = trees.n_conics;
n_lines = 2*n_conics;
iC1 = trees.iC1;
iC2 = trees.iC2;
pts = trees.pts;

% Simulate camera.
t1 = [-2.5;0.5;0.5] + 0.1*randn(3,1);
R1 = randrot([2.5;0;0]+ 0.1*randn(3,1));
P1 = [R1 (-R1*t1)];

if plots || dbug
    figure(30);
    hold on
    pose = rigid3d(R1,t1');
    plotCamera('AbsolutePose',pose,'Opacity',0);
    hold off;
end

% Get lines in image.
[l1, l2] = lines_from_conics(iC1, iC2, P1, n_conics);

% Plot conics and lines.
if dbug
    impts = cell(1, n_conics);
    figure(20); clf;
    for i = 1:n_conics
        impts{i} = pflat(P1*pts{i}); % Project cylindrs into camera
        rita(impts{i},'.'); % Draw cylinders
        hold on;
    end
    axis manual;
    rital([l1 l2]); % Tangent planes (lines in camera)
end
if noise
    l1 = addnoise(l1, noise, dbug);
    l2 = addnoise(l2, noise, dbug);
    if dbug
        rital([l1 l2]);
    end
end

% Find infinity point and rectify.
ll = [l1 l2];
[~,~,V] = svd(ll');
vp = V(:,end); % vanishing point

r2p = vp/norm(vp);
rr = null(r2p');

% Make sure the system is right orthonormal.
Rvp1 = [rr(:,2)';r2p';rr(:,1)'];
Rvp2 = -[rr(:,1)';r2p';rr(:,2)'];
if det(Rvp1)<0
    Rvp1 = [rr(:,1)';r2p';rr(:,2)'];
    Rvp2 = -[rr(:,2)';r2p';rr(:,1)'];
end

llvp1 = Rvp1*ll;
llvp2 = Rvp2*ll;

if dbug
    figure(21)
    title('First rectified sys')
    hold on
    rital(llvp1);
    
    figure(22)
    title('Second rectified sys')
    hold on
    rital(llvp2);
    hold off
end

linenr = randperm(n_lines,3);

% Find xsols, two solvers due to rectification.
[res1, x1] = solution_gen_parallell(llvp1, iC1, linenr);
[res2, x2] = solution_gen_parallell(llvp2, iC1, linenr);

% Residuals for all conics and lines.
res = [res1; res2];

x = [x1 x2];
in_logic = res<bnd;
inliersTemp = sum(in_logic,2);
maxinl = max(inliersTemp);
maxidx = find(inliersTemp==maxinl);

res_inl = res.*in_logic;
[~, min_idx] = min(sum(res_inl(maxidx,:),2));
best_xsol = x(:,maxidx(min_idx));
sys_best = maxidx(min_idx)>8;

a = best_xsol(1); b = best_xsol(2); tz = best_xsol(3);
scale = sqrt(a^2+b^2);
Pestpre1 = [a 0 b 1; 0 scale 0 0; -b 0 a tz]/scale;
Pestpre2 = [a 0 b 1; 0 -scale 0 0; -b 0 a tz]/-scale;

if sys_best == 0
    Rvp = Rvp1;
else
    Rvp = Rvp2;
end

% Four solutions to evaluate.
Porigsys11 = (Rvp1'*Pestpre1);
Porigsys12 = (Rvp1'*Pestpre2);
Porigsys21 = (Rvp2'*Pestpre1);
Porigsys22 = (Rvp2'*Pestpre2);

votum11 = 0;
votum12 = 0;
votum21 = 0;
votum22 = 0;

for idx = 1:n_conics
    c = trees.Tsave{idx}\[trees.c(:,idx); 1];
    proj11 = Porigsys11*c;
    proj12 = Porigsys12*c;
    proj21 = Porigsys21*c;
    proj22 = Porigsys22*c;
    votum11 = votum11 + sign(proj11(end));
    votum12 = votum12 + sign(proj12(end));
    votum21 = votum21 + sign(proj21(end));
    votum22 = votum22 + sign(proj22(end));
end

if sys_best == 0
    [~, bestcamidx] = max([votum11,votum12]);
    if bestcamidx == 1
        bestPest = Porigsys11;
    else % bestcamidx == 2
        bestPest = Porigsys12;
    end
else % sys_best == 1
    [~, bestcamidx] = max([votum21,votum22]);
    if bestcamidx == 1
        bestPest = Porigsys21;
    else % bestcamidx == 2
        bestPest = Porigsys22;
    end
end

% Evaluate solution on non-rectified system.
P_origsys = bestPest;

[l1new, l2new] = lines_from_conics(iC1, iC2, P_origsys, n_conics);
llnew = [l1new l2new];

if dbug
    figure(23);
    clf;
    title('estimated projection')
    imptsnew = cell(1, n_conics);
    for i = 1:n_conics
        imptsnew{i} = pflat(P_origsys*pts{i}); % Project cylindrs into camera
        rita(imptsnew{i},'.'); % Draw cylinders
        hold on;
    end
    axis manual;
    rital(llnew);
end

%Plot new lines into old picture.
if dbug
    figure(20);
    hold on
    rital(llnew);
end


Roterr = acosd((trace(P_origsys(:,1:3)*P1(:,1:3)')-1)/2);
a1 = null(P_origsys);
a2 = null(P1);
a1 = pflat(a1);
a2 = pflat(a2);
transerr = norm(a1([1 3])-a2([1 3]));

% Plot camera centers
if plots || dbug
    cest = pflat(null(P_origsys));
    ctrue = pflat(null(P1));
    figure(30)
    hold on
    plot3(cest(1),cest(2),cest(3),'*k')
    plot3(ctrue(1),ctrue(2),ctrue(3),'*k')
end

% Plot camera
if plots || dbug
    figure(30);
    hold on
    Ros = P_origsys(1:3,1:3);
    tos = P_origsys(1:3,4);
    pose = rigid3d(Ros,(-Ros'*tos)');
    plotCamera('AbsolutePose',pose,'Opacity',0, 'Color', [0, 0, 1]);
    hold off;
    axis equal
end
end
