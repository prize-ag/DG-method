% 코드에 대한 설명
% 버거스 방정식에 대한 exact soln을 구하는 코드가 들어있음
% 초기조건은 변경 가능
% 그리고 포트란으로 돌린 케이스들 가지고 와서, exact soln과 numerical soln을 비교
% l2 l inf 에러를 구할 수 있고, 그에 따른 order 구함 
% 그리고 테이블로 보기 좋게 만들어짐.



clear;
clc;
close all;

% 1. 초기조건 / shock 시간
u0  = @(x) 1.0/4.0 + 1.0/2.0 *sin(pi*x);
du0 = @(x) 1/2*pi*cos(x) ;

t = 0.3;  % shock 생기기 전

% 2. 공간 격자
xL = -1; xR = 1;
length_space  = abs(xR - xL) ;

imax_list = [20 40 80 160 320 640 1280] ;

L2_list    = zeros(size(imax_list)); % L2 에러 설정
Linf_list  = zeros(size(imax_list)); % L inf 에러 설정 
L2_err_ord = zeros(size(imax_list)); % L2 err order 
Linf_err_ord = NaN(size(imax_list)); % L inf err order 


for j = 1:length(imax_list)

    imax = imax_list(j) ;

    temp_x = zeros(1, imax + 1) ; % 경계에서의 값
    x = zeros(1,imax) ; %셀 중심에서의 값

    for i = 1 : imax+1
        temp_x(i) = xL + length_space/imax * (i-1) ;
    end

    for i = 1: imax
        x(i)  = (temp_x(i)+temp_x(i+1))/2 ;

    end

    u_exact = zeros(size(x));

    % 3. 각 x에 대해 implicit eq 풀기
    for i = 1:imax
        xi_guess = x(i);  % 초기 guess는 그냥 x에서 시작
        % x(i) : 궁금한 x에서의 위치, xi :  초기조건(크시)
        F = @(xi) xi + t*u0(xi) - x(i);
        xi_star = fzero(F, xi_guess);
        u_exact(i) = u0(xi_star);
    end

    % %%newton method － 잘못된 해를 찾아줌。。 ㅠ
    % for i = 1:imax
    %     xi_guess = x(i);    % initial guess
    % 
    %     F  = @(xi) xi + t*u0(xi) - x(i);
    %     Fp = @(xi) 1 + t*du0(xi);
    % 
    %     xi_star = newton_burgers(F, Fp, xi_guess, 20, 1e-14);
    % 
    %     u_exact(i) = u0(xi_star);
    % end





    %% save file and read file
    filename = sprintf('ext_result_%d.mat', imax);  
    save(filename,"x","u_exact");

    filename1 = sprintf('Density(Numerical_%d).plt', imax);
    N = load(filename1);   % 또는 dlmread


    %%plot

    % figure;
    % plot(x, u_exact, 'bo', 'LineWidth', 1.3); hold on;
    % plot(N(:,1), N(:,2), 'r*', 'LineWidth', 1.3); hold on;
    % 
    % legend( 'Exact', 'Numerical');
    % 
    % grid on;
    % xlim([-1 1]);     % 전체 영역 보기
    % xlabel('x');
    % ylabel('Density');
    % title('Exact vs Numerical');

    %%plot

    err1 = 0 ;
    err = 0;
    err_vec = N(:,2)' - u_exact;   %0 크기 맞는 오차 벡터

    L2_list(j)   = sqrt(sum(err_vec.^2)) ;
    Linf_list(j) = max(abs(err_vec));

    % fprintf("imax=%d → L2 = %.12e, L∞ = %.12e\n", ...
        % imax, L2_list(j), Linf_list(j));
end

% ---  수렴차수(order) 계산 ---
L2_err_ord(1) = NaN;  % 첫 값은 기준이라 order 없음
for j = 2:length(imax_list)
    L2_err_ord(j) = log2(L2_list(j-1) / L2_list(j));
    Linf_err_ord(j) = log2( Linf_list(j-1) / Linf_list(j) );
end

% Error 부분을 string(exponential)으로 저장
L2_str   = arrayfun(@(x) string(sprintf('%.4e', x)), L2_list, 'UniformOutput', true);
Linf_str = arrayfun(@(x) string(sprintf('%.4e', x)), Linf_list, 'UniformOutput', true);

% Order 계산
Order_L2   = NaN(size(imax_list));
Linf_err_ord = NaN(size(imax_list));

for j = 2:length(imax_list)
    Order_L2(j)   = log2(L2_list(j-1)   / L2_list(j));
    Linf_err_ord(j) = log2(Linf_list(j-1) / Linf_list(j));
end

% 테이블 생성
T = table(imax_list(:), L2_str(:), Linf_str(:), Order_L2(:), Linf_err_ord(:), ...
    'VariableNames', {'Cells','L2_Error','Linf_Error','Order_L2','Order_Linf'});



% T = table(imax_list(:), L2_list(:), Linf_list(:), ...
%           L2_err_ord(:), Linf_err_ord(:), ...
%     'VariableNames', {'Cells', 'L2_Error', 'Linf_Error', 'Order_L2', 'Linf_err_ord'});

disp(T);

