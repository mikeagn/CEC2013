% ******************************************************************************
% * Version: 1.0.1
% * Last modified on: 27 October, 2016 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: m_(DOT)_epitropakis_(AT)_lancaster_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************
%ATTENTION: USE THIS FUNCTION ONLY FOR STATISTICAL PURPOSES
%TODO: Accuracy of the optima positions might vary per function.
%The function returns a number_of_global_optima x number_of_dimension matrix
%each row contains a different optimum position
function [o] = get_copy_of_goptima(nfunc)
total_func_no = 20;

if nfunc > 10 & nfunc <= total_func_no
	load data/optima.mat; % saved the predefined optima, a 10*100 matrix;
	D = get_dimension(nfunc);
	o = o(:,1:D);
	return;
end

if nfunc == 1	    fname = 'data/F1_opt.dat';
elseif nfunc == 2	fname = 'data/F2_opt.dat';
elseif nfunc == 3	fname = 'data/F3_opt.dat';
elseif nfunc == 4	fname = 'data/F4_opt.dat';
elseif nfunc == 5	fname = 'data/F5_opt.dat';
elseif nfunc == 6	fname = 'data/F6_2D_opt.dat';
elseif nfunc == 7	fname = 'data/F6_3D_opt.dat';
elseif nfunc == 8	fname = 'data/F7_2D_opt.dat';
elseif nfunc == 9	fname = 'data/F7_3D_opt.dat';
elseif nfunc == 10	fname = 'data/F8_2D_opt.dat';
else
	fprintf('ERROR: Wrong function number: (%d).\n', nfunc);
	fprintf('       Please provide a function number in {1,2,...,%d}\n', total_func_no);
	fprintf('       For now function number == 1\n');
	fname = '';
end

o = load(fname); % saved the predefined optima
