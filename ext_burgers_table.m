% 한번에 만들려고 했는데 일단 하다가 멈춤.

% 1. 초기조건 / shock 시간
u0  = @(x) 1.0/4.0 + 1.0/2.0 *sin(pi*x);
du0 = @(x) 1/2*pi*cos(x) ;
% Tc  = -1 / (-1);  % min du0 = -1 이라서 Tc = 1

t = 0.3;  % shock 생기기 전

% 2. 공간 격자
xL = -1; xR = 1;
length  = abs(xR - xL) ;

% Nx = 20 ; % 이거는 구할 공간에서의 개수를 의미 ( exact soln은 20개만 나올 예정임)
temp_x = zeros(1, Nx + 1) ; % 경계에서의 값
x = zeros(1,Nx) ; %셀 중심에서의 값

Nx_list = [20 40 80 160 320] ;

u_exact = zeros(size(x),size(Nx_list));

for j = 1 : size(Nx_list)
    Nx  = Nx_list(j) ;

    for i = 1 : Nx+1
        temp_x(i) = xL + length/Nx * (i-1) ;
    end

    for i = 1: Nx
        x(i)  = (temp_x(i)+temp_x(i+1))/2 ;
    end


    % 3. 각 x에 대해 implicit eq 풀기
    for i = 1:Nx
        xi_guess = x(i);  % 초기 guess는 그냥 x에서 시작
        % x(i) : 궁금한 x에서의 위치, xi :  초기조건(크시)
        F = @(xi) xi + t*u0(xi) - x(i);
        xi_star = fzero(F, xi_guess);
        u_exact(i) = u0(xi_star);
    end
end
    % 이제 u_exact이 shock 이전의 정확한 해
    plot(x, u_exact,'o');
    xlabel('x'); ylabel('u(x,t)');
    title(sprintf('Burgers solution at t = %.2f (before shock)', t));
