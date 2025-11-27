clc;
close all; 

E = load('Density(Exact).plt');
% E = load('Density_Exact_160.plt');

N = load('Density(Numerical).plt');
% N = load('Density_Numerical_160.plt');

figure;
plot(E(:,1), E(:,2), 'b-', 'LineWidth', 1.3); hold on;
plot(N(:,1), N(:,2), 'r*', 'LineWidth', 1.3); hold on;

legend( 'Exact', 'Numerical');
% legend( 'Numerical');

grid on;
xlim([-1 1]);     % 전체 영역 보기 
xlabel('x');
ylabel('Density');
title('Exact vs Numerical');
% title( 'Numerical');
