% solve the model and print a graph
function y_params = solveModel()
global rho_empty_1;
global rho_empty_2;
global rho_1_2;
global rho_2_1;
global r_empty_initial;
global r_1_initial;
global r_2_initial;
global r_1_2_initial;
global r_empty_i_1_initial;
global r_empty_i_2_initial;
global r_1_i_2_initial;
global r_2_i_1_initial;
global d_initial;
global days;
global N;
r_empty_initial = 1000 - 70;
r_1_initial = 10;
r_2_initial = 10;
r_1_2_initial = 10;
r_empty_i_1_initial = 10;
r_empty_i_2_initial = 10;
r_1_i_2_initial = 10;
r_2_i_1_initial = 10;
d_initial = 0;

% consts
HOURS_PER_DAY = 24;
days = 28;
initial_condition = [r_empty_initial; r_1_initial; r_2_initial; r_1_2_initial; r_empty_i_1_initial; r_empty_i_2_initial; r_1_i_2_initial; r_2_i_1_initial; d_initial;];
N = r_empty_initial + r_1_initial + r_2_initial + r_1_2_initial + r_empty_i_1_initial + r_empty_i_2_initial + r_1_i_2_initial + r_2_i_1_initial + d_initial;

% solve the ODEs
% assume each step is an hour - very important
x = lsode("mutationSirModel", initial_condition, (t = linspace(0, HOURS_PER_DAY * days, HOURS_PER_DAY * days)' ));

# normalize matrix
for i = 1:size(x)(1)
  for j = 1:size(x)(2)
    x(i, j) = x(i, j) / N;
  endfor
endfor

% colors consts
colors_letters = {"b", "g", "g", "g", "r", "r", "r", "r", "k"};

% plot the results
fig = figure();
hold on;
for i = 1:size(x)[2]
  plot(t, x(:,i), "-");
end
hold off;

% set plot info
xlim ([0, HOURS_PER_DAY * days * 1.05]);
ylim ([0, max(max(x)) * 1.05]);
title ("Mutation Model Dynamics");
xlabel ("Time [hours]");
ylabel ("Normalized population size");
h = legend ("R_{emptyset}", "R_{emptyset}I_1", "R_{emptyset}I_2", "R_1", "R_2", "R_1I_2", "R_2I_1", "R_{1,2}", "D");
set (h, "interpreter", "tex");

# save figure for next runs 
file_name = strcat(save_name ,".jpg");
print (fig, file_name, "-djpg");
printf("\nSaved Dynamics Figure to: '%s'\n", file_name);
# close figures to save memory on long runs
close all

endfunction


function xdot = mutationSirModel(x, t)
global N;

# use parameters
beta_empty_1 = 0.01;
beta_empty_2 = 0.01;
beta_1_2 = 0.1;
beta_2_1 = 0.1;
gama_empty_1 = 0.01;
gama_empty_2 = 0.01;
gama_1_2 = 0.01 ;
gama_2_1 = 0.01 ;
rho_empty_1 = 0.95;
rho_empty_2 = 0.95;
rho_1_2 = 0.75;
rho_2_1 = 0.75;

# model itself
xdot(1) = -1 * ( beta_empty_1 * x(2) + beta_empty_2 * x(3) ) * x(1) ;
xdot(2) = beta_empty_1 * x(2) * x(1) - ( gama_empty_1 * x(2) ) ;
xdot(3) = beta_empty_2 * x(3) * x(1) - ( gama_empty_2 * x(3) ) ;
xdot(4) = (gama_empty_1 * rho_empty_1 * x(2)) - ( beta_1_2 * x(6) * x(4) ) ;
xdot(5) = (gama_empty_2 * rho_empty_2 * x(3)) - ( beta_2_1 * x(7) * x(5) ) ;
xdot(6) = beta_1_2 * x(6) * x(4) - gama_1_2 * x(6) ;
xdot(7) = beta_2_1 * x(7) * x(5) - gama_2_1 * x(7) ;
xdot(8) = gama_2_1 * rho_2_1 * x(7) + gama_1_2 * rho_1_2 * x(6) ;
xdot(9) = gama_empty_1 * (1 - rho_empty_1) * x(2) + gama_2_1 * (1 - rho_2_1) * x(7) + gama_empty_2 *(1 - rho_empty_2) * x(3) + gama_1_2 * (1 - rho_1_2) * x(6) ;

endfunction 
