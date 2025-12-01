%% INITIALIZE WORKSPACE
clear;
clc;
close all;

%% 1. INITIAL CONDITION / SHOCK TIME
u0 = @( x ) 1.0 / 4.0 + 1.0 / 2.0 * sin( pi * x );

% SHOCK-FREE TIME (PRE-SHOCK)
t = 0.5 / pi ^ 2; 

%% 2. LIST OF IMAX (NUMBER OF CELLS)
imax_list  = [ 20, 40, 80, 160, 320, 640, 1280 ];
num_levels = length( imax_list );

L2_list    = zeros( 1, num_levels );  % L2 ERROR STORAGE
Linf_list  = zeros( 1, num_levels );  % LINF ERROR STORAGE
Order_L2   = NaN( 1, num_levels );    % L2 CONVERGENCE ORDER
Order_Linf = NaN( 1, num_levels );    % LINF CONVERGENCE ORDER

%% 3. LOOP OVER GRID LEVELS (READ FORTRAN RESULTS AND COMPUTE ERRORS)
for j = 1 : num_levels

    imax = imax_list( j );

    %% 3-1. LOAD NUMERICAL SOLUTION FROM FORTRAN (TXT FILE)
    filename_num = sprintf( 'DG1D_NUMERICAL_%05d.TXT', imax );

    data = load( filename_num );    % COLUMN 1: X, COLUMN 2: NUMERICAL SOLUTION

    x_num  = data( :, 1 );
    u_num  = data( :, 2 );

    % ENSURE ROW VECTORS
    x      = x_num.';               
    u_num  = u_num.';              
    imax   = length( x );

    % UNIFORM GRID SIZE (ASSUMING UNIFORM CELLS)
    dx_val = x( 2 ) - x( 1 );       
    dx     = dx_val * ones( 1, imax );

    %% 3-2. COMPUTE EXACT SOLUTION AT THE SAME GRID POINTS
    u_exact = zeros( size( x ) );

    % SOLVE IMPLICIT EQUATION FOR BURGERS' EQUATION:
    %   xi + t * u0( xi ) = x( i )
    for i = 1 : imax
        xi_guess = x( i ); % INITIAL GUESS

        F       = @( xi ) xi + t * u0( xi ) - x( i );
        xi_star = fzero( F, xi_guess );

        u_exact( i ) = u0( xi_star );
    end

    %% 3-3. COMPUTE ERRORS (EXACT VS NUMERICAL)
    err_vec = u_num - u_exact;              

    % L2 ERROR (UNIFORM GRID: DX CONSTANT)
    L2_local = err_vec .^ 2;
    L2_local = L2_local * dx_val;

    L2_list( j )   = sqrt( sum( L2_local ) );
    Linf_list( j ) = max( abs( err_vec ) );

end

%% 4. COMPUTE CONVERGENCE ORDERS
% FIRST ENTRY HAS NO ORDER (REFERENCE LEVEL)
Order_L2( 1 )   = NaN;
Order_Linf( 1 ) = NaN;

for j = 2 : num_levels
    Order_L2( j )   = log2( L2_list( j - 1 )   / L2_list( j ) );
    Order_Linf( j ) = log2( Linf_list( j - 1 ) / Linf_list( j ) );
end

%% 5. FORMAT ERRORS AS STRINGS (EXPONENTIAL FORMAT, WITH 2 DECIMALS)
L2_str   = arrayfun( @( x ) string( sprintf( '%.2e', x ) ), L2_list,   'UniformOutput', true );
Linf_str = arrayfun( @( x ) string( sprintf( '%.2e', x ) ), Linf_list, 'UniformOutput', true );

%% 6. FORMAT ORDERS FOR DISPLAY
Order_L2_disp   = strings( num_levels, 1 );
Order_Linf_disp = strings( num_levels, 1 );

for j = 1 : num_levels

    if isnan( Order_L2( j ) )
        Order_L2_disp( j ) = " - ";
    else
        Order_L2_disp( j ) = sprintf( "%.2f", Order_L2( j ) );
    end

    if isnan( Order_Linf( j ) )
        Order_Linf_disp( j ) = " - ";
    else
        Order_Linf_disp( j ) = sprintf( "%.2f", Order_Linf( j ) );
    end

end

%% 7. BUILD SUMMARY TABLE
T = table( imax_list( : ), ...
           L2_str( : ), ...
           Linf_str( : ), ...
           Order_L2_disp, ...
           Order_Linf_disp, ...
           'VariableNames', { 'Cells', 'L^2 Error', 'L^Inf Error', 'L^2 Order', 'L^Inf Order' } );

disp( T );

%% 8. BUILD LATEX TABULAR STRING AND SAVE TO FILE
latex_lines = strings( num_levels + 4, 1 );
idx         = 1;

latex_lines( idx ) = "\begin{tabular}{rccccc}"; idx = idx + 1;
latex_lines( idx ) = "\hline";                   idx = idx + 1;
latex_lines( idx ) = "N & $L^2$ error & $L^\infty$ error & $L^2$ Order & $L^\infty$ Order\\"; 
idx = idx + 1;
latex_lines( idx ) = "\hline";                   idx = idx + 1;

for j = 1 : num_levels

    row_str = sprintf( "%4d & %s & %s & %s & %s \\", ...
                       imax_list( j ), ...
                       L2_str( j ), ...
                       Linf_str( j ), ...
                       Order_L2_disp( j ), ...
                       Order_Linf_disp( j ) );

    latex_lines( idx ) = row_str;
    idx = idx + 1;

end

latex_lines( idx ) = "\hline";       idx = idx + 1;
latex_lines( idx ) = "\end{tabular}";

latex_table = strjoin( latex_lines, newline );

fprintf( '\nLATEX TABLE:\n\n' );
disp( latex_table );

% OPTIONAL: WRITE TO FILE
tex_filename = 'DG1D_Burgers_error_table.tex';
fid = fopen( tex_filename, 'w' );
if fid ~= -1
    fprintf( fid, '%s\n', latex_table );
    fclose( fid );
    fprintf( '\nLATEX TABLE WRITTEN TO FILE: %s\n', tex_filename );
else
    fprintf( '\nWARNING: COULD NOT OPEN FILE %s FOR WRITING.\n', tex_filename );
end
