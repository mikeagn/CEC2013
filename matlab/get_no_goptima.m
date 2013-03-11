% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2013 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

function [no] = get_no_goptima(nfunc)
nopt = [2 5 1 4 2 18 36 81 216 12 6 8 6 6 8 6 8 6 8 8];
no = nopt(nfunc);
