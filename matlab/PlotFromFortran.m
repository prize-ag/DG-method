clc; 
close all; 
% 포트란에서 가져온 데이터를 이미지화 시킨 것!


% data = load('result.dat8');  % 파일 읽기(sin pi x ) 
% data = load('result6_박사님코드박스.dat');  % 파일 읽기(box)

x = data(:,1);              % 첫 번째 열
u = data(:,2);              % 두 번째 열


% ext = sin(pi*x); 
ext = zeros(size(x)) ; 
ext(x >= -0.5 & x <= 0.5) = 1.0;

plot(x, u,'ko', 'LineWidth', 1.5);
hold on; 
plot(x,ext,'LineWidth', 2); 

xlabel('x');
ylabel('u(x)');
title('Numerical Solution from Fortran Output');
grid on;