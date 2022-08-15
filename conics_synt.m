function [Roterr, transerr] = conics_synt(noise, plots, seed, dbug)

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

clearvars -except noise plots seed dbug; 
close all
get_paths

parallell = 0;

%Load synthetic data.
trees = load_cylinders(parallell, plots, dbug); % Struct with everyting needed
[lines, trueP] = load_camera(trees, plots, noise, dbug);

disp('Solving rotation')
linenr = randperm(trees.n_conics*2,trees.n_conics);
[qsols] = solution_gen_rotation([lines.l1 lines.l2], trees.NN, linenr);

% Translation part.
disp('Solving translation')
[ts] = solution_gen_translation([lines.l1 lines.l2], trees.iC1, qsols, linenr);

if isempty(qsols) || isempty(ts)
    disp('no real solutions')
    Roterr = NaN;
    transerr = NaN;
    return;
end

% Picking best solution.
nqs = size(qsols,2);
bestPest = [];
besterr = inf;
for i = 1:nqs
   for j = 1:size(ts{i},2)
       Pest = [quat2rot(qsols(:,i)) ts{i}(:,j)];
       err = norm(Pest-trueP,'fro');
       if err < besterr
           bestPest = Pest;
           besterr = err;
       end
   end
end

disp('Estimated P:')
disp(bestPest)
disp('True P:')
disp(trueP)

% Evaluate.
Roterr = acosd((trace(bestPest(:,1:3)*trueP(:,1:3)')-1)/2);
a1 = null(bestPest); 
a2 = null(trueP);
a1 = pflat(a1); 
a2 = pflat(a2);
transerr = norm(a1-a2);

if plots || dbug
    figure(30);
    hold on
    Ros = bestPest(1:3,1:3);
    tos = bestPest(1:3,4);
    pose = rigid3d(Ros,(-Ros'*tos)');
    plotCamera('AbsolutePose',pose,'Opacity',0, 'Color', [0, 0, 1]);
    hold off;
end

end

function [lines, trueP] = load_camera(trees, plots, noise, dbug)
% Simulate camera
t1 = [-7;0;0] + 0.1*randn(3,1);
R1 = randrot([2.5;0;0]+ 0.1*randn(3,1));
P1 = [R1 (-R1*t1)]; % Camera matrix

if plots || dbug
    figure(30);
    hold on
    pose = rigid3d(R1,t1');
    plotCamera('AbsolutePose',pose,'Opacity',0);
    hold off;
end

[l1, l2] = lines_from_conics(trees.iC1, trees.iC2, P1, trees.n_conics); %(iC1, iC2, P, n_trees)
trueP = P1;
impts = cell(1, trees.n_conics);
if dbug
    figure(20);
    for i = 1:trees.n_conics
        impts{i} = pflat(P1*trees.pts{i}); % Project cylindrs into camera.
        rita(impts{i},'.'); % Draw cylinders.
        hold on;
    end
    axis manual;
    rital([l1 l2]); % Tangent planes (lines in camera)
end

lines.l1 = l1;
lines.l2 = l2;

if noise
    disp('adding noise')
    lines.l1 = addnoise(l1,noise, dbug);
    lines.l2 = addnoise(l2,noise, dbug);
    
    if dbug
        rital([lines.l1 lines.l2],'--');
    end
end
end
