% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2013 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

function [rho] = get_rho(nfunc)
rho_ = [0.01*ones(1,4) 0.5 0.2 0.5 0.2 0.01*ones(1,11)];
rho = rho_(nfunc);
