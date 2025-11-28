clc;
close all;

% E = load('Density(Exact).plt');
% E = load('Density_Exact_160.plt');

N = load('Density(Numerical).plt');
% N = load('Density_Numerical_160.plt');

figure;
plot(x, u_exact, 'bo', 'LineWidth', 1.3); hold on;
plot(N(:,1), N(:,2), 'r*', 'LineWidth', 1.3); hold on;

legend( 'Exact', 'Numerical');
% legend( 'Numerical');

grid on;
xlim([-1 1]);     % 전체 영역 보기
xlabel('x');
ylabel('Density');
title('Exact vs Numerical');
% title( 'Numerical');

err = 0 ; 
for i= 1: size(x)
    temp = 0; 
    err = err + (N(i,2)-u_exact(1,i)).^2 ; 
end
L2_err = sqrt(sum(err)) ; 
disp(L2_err)