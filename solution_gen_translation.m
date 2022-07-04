function x = solution_gen_translation(ll, iC1, qsols, linenr)

n_conics = length(iC1);
modlinenr = mod(linenr-1,n_conics) + 1;

% Take 3 lines.
l1var = ll(:,linenr(1));
l2var = ll(:,linenr(2));
l3var = ll(:,linenr(3));

% Get corresponing representations.
iC1_1  = iC1{modlinenr(1)};
iC1_2  = iC1{modlinenr(2)};
iC1_3  = iC1{modlinenr(3)};

%Restructure the c matrices.
c1 = [iC1_1(1,:) iC1_1(2,2:4)  iC1_1(3,3:4)]';
c2 = [iC1_2(1,:) iC1_2(2,2:4)  iC1_2(3,3:4)]';
c3 = [iC1_3(1,:) iC1_3(2,2:4)  iC1_3(3,3:4)]';

% Using the solutions for the rotation, generate solutions for the
% translation.
nqs = size(qsols,2);
x = cell(1,nqs);
for i = 1:nqs
    data1 = [l1var; l2var; l3var; qsols(:,i); c1; c2; c3];
    tsols = solver_conics_translation(data1);
    tsols = tsols(:,sum(abs(imag(qsols)))<1e-8);
    x{i} = real(tsols);
end
