figure;

for lev = 1:5
    fname = sprintf('DG1D_GRIDLEVEL_%d.TXT', lev);
    data  = load(fname);
    
    x      = data(:, 1);
    u_ref  = data(:, 2);
    u_num  = data(:, 3);
    
    subplot(5, 1, lev);           
    plot(x, u_ref, 'k-', 'LineWidth', 1.0); hold on;
    plot(x, u_num, 'r-', 'LineWidth', 1.0);
    
    if lev == 5
        xlabel('x');
    end
    ylabel(sprintf('u', lev));
    
    title(sprintf('Grid Level %d', lev));
    grid on;
    
    if lev == 1
        legend('Reference', 'DG', 'Location', 'north');
    end
end
