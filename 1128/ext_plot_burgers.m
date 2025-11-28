% shock이 발생하기 전, 버거스 방정식에 대해서 exact soln 찾고,
% 그 그래프를 그려주는 코드

clear; 
clc; 
close all; 

% 1. 초기조건 / shock 시간
u0  = @(x) 1.0/4.0 + 1.0/2.0 *sin(pi*x);
% du0 = @(x) 1/2*pi*cos(x) ;

t = 0.3;  % shock 생기기 전

% 2. 공간 격자
xL = -1; xR = 1;
length  = abs(xR - xL) ; 

Nx = 80 ; % 이거는 구할 공간에서의 개수를 의미 ( exact soln은 20개만 나올 예정임) 
temp_x = zeros(1, Nx + 1) ; % 경계에서의 값
x = zeros(1,Nx) ; %셀 중심에서의 값 

for i = 1 : Nx+1 
    temp_x(i) = xL + length/Nx * (i-1) ;
end 

for i = 1: Nx 
    x(i)  = (temp_x(i)+temp_x(i+1))/2 ; 
end

u_exact = zeros(size(x));

% 3. 각 x에 대해 implicit eq 풀기
for i = 1:Nx
    xi_guess = x(i);  % 초기 guess는 그냥 x에서 시작 
    % x(i) : 궁금한 x에서의 위치, xi :  초기조건(크시)
    F = @(xi) xi + t*u0(xi) - x(i);
    xi_star = fzero(F, xi_guess);
    u_exact(i) = u0(xi_star);
end

save('exact_burgers',"x","u_exact");

% 이제 u_exact이 shock 이전의 정확한 해
plot(x, u_exact,'o');
xlabel('x'); ylabel('u(x,t)');
title(sprintf('Burgers solution at t = %.2f (before shock)', t));
