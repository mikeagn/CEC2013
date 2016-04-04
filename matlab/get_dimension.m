% ******************************************************************************
% * Version: 1.0
% * Last modified on: 4 April, 2016 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

function [fit] = get_dimension(nfunc)
Dims = [1 1 1 2 2 2 2 3 3 2 2 2 2 3 3 5 5 10 10 20]; % dimensionality of benchmark functions
fit = Dims(nfunc);
