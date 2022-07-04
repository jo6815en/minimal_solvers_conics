function trees = load_cylinders(parallell, plots, dbug)

% Simulate 3D cylinders
n_conics = 4;
c = rand(3,n_conics); % x,y,z coordinates for n "center"-points
r = (0.1*rand(2,n_conics)+0.2);

for i = 1:n_conics
    if parallell
        R = [0 0 1; 1 0 0 ; 0 1 0];
    else
        [R,~,~]=svd(randn(3,3)); % Getting a random rotation matrix
    end
    T = [R zeros(3,1);zeros(1,3) 1]; % Avbildningsmatris m rot
    r1 = r(1,i); % Radius
    r2 = r(2,i);
    ir1 = 1/(r1^2);
    ir2 = 1/(r2^2);
    c1 = c(1,i);
    c2 = c(2,i);
    C1{i}=[diag([ir1 ir2 0]) (-[c1*ir1;c2*ir2;0]); ...
        (-[c1*ir1;c2*ir2;0]') (c1^2*ir1 + c2^2*ir2 -1)];
    
    C2{i}=diag([0 0 1 0]); % The plane z=0
    iC1tmp = zeros(4,4);
    iC1tmp([1 2 4],[1 2 4]) = inv(C1{i}([1 2 4],[1 2 4])); % Samma sak som att titta på problemet utan z-axel. Där är problemet ej degenererat och dualen är C^-1.
    iC1{i} = iC1tmp;
    iC2{i} = C2{i};
    C1{i} = T'*C1{i}*T; % The C-matrix containing the restrictions of tree-points.
    %OBS! inv(T)=T' jfr T=inv(A) C = inv(A')*T*inv(A) <=> C1 '=' T'*C1*T
    C2{i} = T'*C2{i}*T;
    iC1{i} = inv(T)*iC1{i}*inv(T'); % Dual of C, containing the tangent planes
    iC2{i} = inv(T)*iC2{i}*inv(T');
    
    [th,z]=meshgrid(0:0.1:(2*pi),-1:0.1:3);
    X = r1*cos(th) + c1; % Rotation och translation
    Y = r2*sin(th) + c2;
    Z = z;
    A = cos(th)/r1;
    B = sin(th)/r2;
    C = zeros(size(z));
    D = - ( A.*X + B.*Y + C.*Z );
    trees.pts{i} = T'*[X(:) Y(:) Z(:) ones(size(Z(:)))]'; % Cylinder skrew. jfr x = A*x_t
    trees.planes{i} = inv(T)*([A(:) B(:) C(:) D(:)]'); % Planes for every point jfr pi = C*x (?)
    trees.Tsave{i} = T;
    
    [~,~,V]=svd(iC2{i});
    NN(:,i)=V(:,1);
    
    if plots || dbug
        figure(30); title('Cylinders in 3D')
        hold on;
        rita3(trees.pts{i},'*'); hold on;
        
        if false
            for k = 1:size(trees.pts{i},2)
                U = trees.pts{i}(:,k);
                Pi = trees.planes{i}(:,k);
                blubb1(k)=U'*C1{i}*U;     % Ekvationen för punkterna
                blubb2(k)=Pi'*iC1{i}*Pi;  % Ekvationen för tangentplanen
                blubb3(k)=Pi'*iC2{i}*Pi;  % Ekvationen för tangentplanen
                blubb4(k)=Pi'*U;          % Punkten ligger i planet
            end
            norm(blubb1);
            norm(blubb2);
            norm(blubb3);
            norm(blubb4);
        end
    end
end

trees.n_conics = n_conics;
trees.c = c;
trees.NN = NN;
trees.iC1 = iC1;
trees.iC2 = iC2;

end
