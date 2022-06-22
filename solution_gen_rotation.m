function x = solution_gen_rotation(ll, NN, linenr)

n_conics = size(NN,2);
modlinenr = mod(linenr-1,n_conics) + 1;

l1var = ll(:,linenr(1));
l2var = ll(:,linenr(2));
l3var = ll(:,linenr(3));
r1var = NN(1:3,modlinenr(1));
r2var = NN(1:3,modlinenr(2));
r3var = NN(1:3,modlinenr(3));

data1 = [l1var l2var l3var r1var r2var r3var];

qsols = solver_conics_rotation(data1(:));

qsols = qsols(:,sum(abs(imag(qsols)))<1e-8);
qsols = real(qsols);
% Normalize quarternions
qsols = [qsols; ones(1,size(qsols,2))]; % Get full quart
% length_quart = vecnorm(qsols); % Calculating the norm of "pre-quarterninon"
qsols = qsols./sqrt(sum(qsols.^2)); %./repmat(length_quart,[4 1]);  % q = q/norm(q);  norm(q)=1

x = qsols;