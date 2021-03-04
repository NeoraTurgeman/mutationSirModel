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

N = 10000;
r_1_initial = 1;
r_2_initial = 1;
r_1_2_initial = 0;
r_empty_i_1_initial = 0;
r_empty_i_2_initial = 0;
r_1_i_2_initial = 0;
r_2_i_1_initial = 0;
d_initial = 0;
r_empty_initial = N - (r_1_initial + r_2_initial + r_1_2_initial + r_empty_i_1_initial + r_empty_i_2_initial + r_1_i_2_initial + r_2_i_1_initial + d_initial);

% consts
HOURS_PER_DAY = 24;
days = 10;
initial_condition = [r_empty_initial; r_1_initial; r_2_initial; r_1_2_initial; r_empty_i_1_initial; r_empty_i_2_initial; r_1_i_2_initial; r_2_i_1_initial; d_initial;];

% solve the ODEs
% assume each step is an hour - very important
x = lsode("mutationSirModel", initial_condition, (t = linspace(0, HOURS_PER_DAY * days, HOURS_PER_DAY * days)' ));

# normalize matrix
for i = 1:size(x)(1)
  for j = 1:size(x)(2)
    x(i, j) = x(i, j) / N;
  endfor
endfor

% plot the results
fig = figure();

% set plot info

hold on;

plot(t, x(:,1), "-b");
plot(t, x(:,2), "-r");
plot(t, x(:,3), "-.r");
plot(t, x(:,4), "-g");
plot(t, x(:,5), "-.g");
plot(t, x(:,6), "--r");
plot(t, x(:,7), "--.r");
plot(t, x(:,8), "-g");
plot(t, x(:,9), "-k");

xlim([0, HOURS_PER_DAY * days * 1.05]);
ylim([0, max(max(x)) * 1.05]);
title("Mutation Model Dynamics");
xlabel("Time [hours]");
ylabel("Normalized population size");
legend({"R_{0}", "R_{0}I_1", "R_{0}I_2", "R_1", "R_2", "R_1I_2", "R_2I_1", "R_{1,2}", "D"});

hold off;

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
beta_empty_1 = 0.002;
beta_empty_2 = 0.001;
beta_1_2 = 0.001;
beta_2_1 = 0.002;
gama_empty_1 = 0.05;
gama_empty_2 = 0.03;
gama_1_2 = 0.000008 ;
gama_2_1 = 0.00001 ;
rho_empty_1 = 0.98;
rho_empty_2 = 0.98;
rho_1_2 = 0.90;
rho_2_1 = 0.95;

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
