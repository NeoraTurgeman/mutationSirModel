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

	# use parameters
	beta_empty_1 = 0.0002/check  ;
	beta_empty_2 = 0.00017/check ;
	beta_1_2 = 0.0017/check  ;
	beta_2_1 = 0.002/check ;
	gama_empty_1 = 0.19/check;
	gama_empty_2 = 0.099/check;
	gama_1_2 = 0.15/check ;
	gama_2_1 = 0.2/check ;
	rho_empty_1 = 0.98 ;
	rho_empty_2 = 0.98 ;
	rho_1_2 = 0.95 ;
	rho_2_1 = 0.95 ;

	N = 10000;
	r_1_initial = 1;
	r_2_initial = 1;
	r_1_2_initial = 1;
	r_empty_i_1_initial = 1;
	r_empty_i_2_initial = 1;
	r_1_i_2_initial = 1;
	r_2_i_1_initial = 1;
	d_initial = 0;
	r_empty_initial = N - (r_1_initial + r_2_initial + r_1_2_initial + r_empty_i_1_initial + r_empty_i_2_initial + r_1_i_2_initial + r_2_i_1_initial + d_initial);
	
	x = calcModel();
	printModel(x);
endfunction

% solve the model and calc sensitiviy analysis
function y_params = sensitivityAnalysis()
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
	r_1_2_initial = 1;
	r_empty_i_1_initial = 1;
	r_empty_i_2_initial = 1;
	r_1_i_2_initial = 1;
	r_2_i_1_initial = 1;
	d_initial = 0;
	r_empty_initial = N - (r_1_initial + r_2_initial + r_1_2_initial + r_empty_i_1_initial + r_empty_i_2_initial + r_1_i_2_initial + r_2_i_1_initial + d_initial);

	# use parameters
	beta_empty_1 = 0.0002/check  ;
	beta_empty_2 = 0.00017/check ;
	beta_1_2 = 0.0017/check  ;
	beta_2_1 = 0.002/check ;
	gama_empty_1 = 0.19/check;
	gama_empty_2 = 0.099/check;
	gama_1_2 = 0.15/check ;
	gama_2_1 = 0.2/check ;
	rho_empty_1 = 0.98 ;
	rho_empty_2 = 0.98 ;
	rho_1_2 = 0.95 ;
	rho_2_1 = 0.95 ;
	
	# play with that, it will be easy to control the graph
	start_val = 0;
	end_val = 0.0002;
	var_number_of_steps = 10;
	val_step_size = (end_val - start_val) / var_number_of_steps;
	
	sens = sensitivitySingleRun(start_val, val_step_size, end_val);
	printSensitivity(sens);
endfunction


function x = calcModel()
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
	
	% consts
	check = 4;
	HOURS_PER_DAY = 24;
	days = 40;
	initial_condition = [r_empty_initial; r_1_initial; r_2_initial; r_1_2_initial; r_empty_i_1_initial; r_empty_i_2_initial; r_1_i_2_initial; r_2_i_1_initial; d_initial;];

	% solve the ODEs
	% assume each step is an hour - very important
	x = lsode("mutationSirModel", initial_condition, (t = linspace(0, days * check , days * check )' ));

	# normalize matrix
	for i = 1:size(x)(1)
	  for j = 1:size(x)(2)
		x(i, j) = x(i, j) / N;
	  endfor
	endfor
endfunction


function x = printModel(x)
	% plot the results
	fig = figure();

	% set plot info

	hold on;

	plot(t/check, x(:,1), "-b");
	plot(t/check, x(:,2), "-r");
	plot(t/check, x(:,3), "-.r");
	plot(t/check, x(:,4), "-g");
	plot(t/check, x(:,5), "-.g");
	plot(t/check, x(:,6), "-m");
	plot(t/check, x(:,7), "-.m");
	plot(t/check, x(:,8), "-c");
	plot(t/check, x(:,9), "--k");

	xlim([0, days * 1.05]);
	ylim([0, max(max(x)) * 1.05]);
	title("Mutation Model Dynamics");
	xlabel("Time [Days]");
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


function sens = printSensitivity(sens, parm_name)
	fig = figure();

	% set plot info
	hold on;

	plot(sens(:, 1), x(:,2), "-k");

	xlim([0, max(sens(:, 1)) * 1.05]);
	ylim([0, max(x(:, 2)) * 1.05]);
	title(strcat("Mutation Model sensitivy for", parm_name));
	xlabel(strcat(strcat("Parameter ", parm_name), " value"));
	ylabel("NEORA - THE METRIC YOU GATHER");
	hold off;

	# save figure for next runs 
	file_name = strcat(save_name ,"_sensitiy.jpg");
	print (fig, file_name, "-djpg");
	printf("\nSaved Dynamics Figure to: '%s'\n", file_name);
	# close figures to save memory on long runs
	close all
endfunction


function xdot = mutationSirModel(x, t)
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
	check = 4;

	# model itself
	xdot(1) = -1 * ( beta_empty_1 * x(2) + beta_empty_2 * x(3) ) * x(1) ;
	xdot(2) = beta_empty_1 * x(2) * x(1) - ( gama_empty_1 * x(2) ) ;
	xdot(3) = beta_empty_2 * x(3) * x(1) - ( gama_empty_2 * x(3) ) ;
	xdot(4) = gama_empty_1 * rho_empty_1 * x(2) - beta_1_2 * x(4) * x(6) ;
	xdot(5) = gama_empty_2 * rho_empty_2 * x(3) - beta_2_1 * x(5) * x(7);
	xdot(6) = (beta_1_2 * x(6) * x(4) ) - (gama_1_2 * x(6) ) ;
	xdot(7) = (beta_2_1 * x(7) * x(5) ) - (gama_2_1 * x(7) ) ;
	xdot(8) = (gama_2_1 * rho_2_1 * x(7)) + (gama_1_2 * rho_1_2 * x(6)) ;
	xdot(9) = gama_empty_1 * (1 - rho_empty_1) * x(2) + gama_2_1 * (1 - rho_2_1) * x(7) + gama_empty_2 *(1 - rho_empty_2) * x(3) + gama_1_2 * (1 - rho_1_2) * x(6) ;

endfunction


function sens = sensitivitySingleRun(start_val, val_step_size, end_val)
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
	
	sens = [];
	index = 1;
	for val = start_val:end_val:val_step_size
		# NEORA - play here with the var you are exam
		rho_1_2 = val;
		x = calc_model();
		
		# find the target var
		sens(index, 1) = val;
		sense(index, 2) = targetVar(x);
    endfor 
	
endfunction


function answer = targetVar(x)
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
	
	# extract some info from the dynamics 
	# TODO - this is just an example, you probably want to use something else
	answer = max(x(:,2)) + max(x(:,3));
	
endfunction
