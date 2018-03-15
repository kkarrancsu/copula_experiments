%% Weak Bivariate, Strong Tri-Variate RV Dependency Generation
% Experiment #1
%   Dep(X1,X2) = weak
%   Dep(X1,X3) = weak
%   Dep(X2,X3) = weak
%   Dep(X1,X2,X3) = strong
%
% Model:  X1+X2+X3 = 1  (Co-monotonic)
%         Define: U1 = F1(X1), U2 = F2(X2), U3 = F3(X3)
%         Bivariate Association Conditions
%         Pi(U1,U2) < C(U1,U2) < M(U1,U2)
%         Pi(U1,U3) < C(U1,U3) < M(U1,U3)
%         Pi(U2,U3) < C(U2,U3) < M(U2,U3)
%
% Simulation Procedure:
%  1.) X1 ~ U(0,1), X2 ~ U(0,1)
%  2.) X3 = 1 - (X1 + X2)
%  3.) Compute U1, U2, U3, C(U1,U2), C(U1,U3), C(U2,U3),
%                         Pi(U1,U2),Pi(U1,U3),Pi(U2,U3),
%                          M(U1,U2), M(U1,U3), M(U2,U3)
%  4.) Throw out points which don't satisfy the constraints

clear;
clc;

u_delta = 0.005;

uu = 0:u_delta:1;
[X1,X2] = meshgrid(uu,uu);
X1 = X1(:); X2 = X2(:);
X3 = 1 - (X1 + X2);

U1 = pobs(X1); U2 = pobs(X2); U3 = pobs(X3);
C12 = empcopula(U1,U2); Pi12 = U1.*U2; M12 = min(U1,U2); W12 = max(U1+U2-1,0);
fprintf('Step 1 complete!\n');
C13 = empcopula(U1,U3); Pi13 = U1.*U3; M13 = min(U1,U3); W13 = max(U1+U3-1,0);
fprintf('Step 2 complete!\n');
C23 = empcopula(U2,U3); Pi23 = U2.*U3; M23 = min(U2,U3); W23 = max(U2+U3-1,0);
fprintf('Step 3 complete!\n');

select_flag = ( ((C12 >= Pi12) & (C12 < M12)) & ...
                ((C23 >= Pi23) & (C23 < M23)) & ...
                ((C13 >= Pi13) & (C23 < M23)) );% | ...
%               ( ((C12 <= Pi12) & (C12 > W12)) & ...
%                 ((C23 <= Pi23) & (C23 > W23)) & ...
%                 ((C13 <= Pi13) & (C23 > W23)) ) ;

x1 = X1(select_flag);
x2 = X2(select_flag);
x3 = X3(select_flag);