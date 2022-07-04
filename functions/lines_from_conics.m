function [l1, l2] = lines_from_conics(iC1, iC2, P, n_trees)

% Fix lines
l1 = zeros(3,n_trees);
l2 = zeros(3,n_trees);
%Project 3D conics to 2D conics
for i = 1:n_trees 
    ic1 = P*iC1{i}*P';
    ic2 = P*iC2{i}*P';
    [~,~,v]=svd(ic2); % s(2:3) = 0
    b1 = v(:,2);% b1 och b2 är en bas för nollrummet till oändlighetspunkten = alla linjer som går igenom oändlighetspunkten
    b2 = v(:,3);
    cc = [b1'*ic1*b1 2*b1'*ic1*b2 b2'*ic1*b2]; % (b1'x + b2)C(b1'x + b2)'=0 för linjer som tangerar cylindern
    ll = roots(cc);
    l1(:,i) = b1*ll(1)+b2;
    l2(:,i) = b1*ll(2)+b2;
end


end