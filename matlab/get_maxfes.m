% ******************************************************************************
% * Version: 1.0
% * Last modified on: 4 April, 2016 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

function [fit] = get_maxfes(nfunc)
Max_FEs = [50000*ones(1,5) 200000 200000 400000 400000 200000*ones(1,4) 400000*ones(1,7)]; 
fit = Max_FEs(nfunc);
