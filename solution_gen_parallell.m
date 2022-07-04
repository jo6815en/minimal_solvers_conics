function [err, x] = solution_gen_parallell(llvp, iC1, linenr)

n_conics = length(iC1);
modlinenr = mod(linenr-1,n_conics) + 1;

% Normalize
llvp = llvp./llvp(1,:);

% Take 3 lines.
l1 = llvp(3,linenr(1));
l2 = llvp(3,linenr(2));
l3 = llvp(3,linenr(3));

% Get corresponing representations.
iC1_1  = iC1{modlinenr(1)};
iC1_2  = iC1{modlinenr(2)};
iC1_3  = iC1{modlinenr(3)};

%Restructure the c matrices.
c1 = [iC1_1(1,[1 3 4]) iC1_1(3,[3 4])]';
c2 = [iC1_2(1,[1 3 4]) iC1_2(3,[3 4])]';
c3 = [iC1_3(1,[1 3 4]) iC1_3(3,[3 4])]';


%Sort data correctly
data = [l1; l2; l3; c1; c2; c3];

% Generate solutions.
x = solver_parallell_conics(data);

% Calculate the error given the 8 solutions.
err = zeros(8,size(llvp,2));
tmax = length(llvp);
for s = 1:8 
    xx = x(:,s);
    a = xx(1); b = xx(2); tz = xx(3);
    scale = sqrt(a^2+b^2);
    P = [a b 1; -b a tz]/scale;
    
    for t = 1:tmax
        cinv = iC1{mod(t-1,tmax/2)+1}([1 3 4], [1 3 4]);
        l = llvp([1 3],t);
        l = l/norm(l); 
        err(s,t) = abs(l'*P*cinv/norm(cinv)*P'*l); % Does pi meet the dual function.
        
    end
end
end
