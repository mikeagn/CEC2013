% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2013 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

function fit = niching_func(x,func_num)
% Benchmark Functions for CEC'2013 Special Session and Competition on 
%        Niching Methods for Multimodal Function Optimization
%
% INPUT:  x :	 	is a 1xD input vector for evaluation
%         func_num : 	denotes the number of the objective function which 
%         		is going to be used.
%
% OUTPUT: fit : 	The objective function value of the x input vector.
%
% This benchmark set includes the following 12 multimodal test functions:
%F1 : Five-Uneven-Peak Trap (1D)
%F2 : Equal Maxima (1D)
%F3 : Uneven Decreasing Maxima (1D)
%F4 : Himmelblau (2D)
%F5 : Six-Hump Camel Back (2D)
%F6 : Shubert (2D, 3D)
%F7 : Vincent (2D, 3D)
%F8 : Modified Rastrigin - All Global Optima (2D)
%F9 : Composition Function 1 (2D)
%F10 : Composition Function 2 (2D)
%F11 : Composition Function 3 (2D, 3D, 5D, 10D)
%F12 : Composition Function 4 (3D, 5D, 10D, 20D)
%
% For more information please refer to the Technical Report of the 
% Special Session/Competition
% at: http://goanna.cs.rmit.edu.au/~xiaodong/cec13-niching/
%
% This source code is based on the following two works:
% P. N. Suganthan, N. Hansen, J. J. Liang, K. Deb, Y. P. Chen, A. Auger, and S. Tiwari, 
% "Problem definitions and evaluation criteria for the CEC 2005 special 
% session on real-parameter optimization," Nanyang Technological University 
% and KanGAL Report #2005005, IIT Kanpur, India., Tech. Rep., 2005.
% and
% B.-Y. Qu and P. N. Suganthan, "Novel multimodal problems and differential evolution 
% with ensemble of restricted tournament selection," in Proceedings of the 
% IEEE Congress on Evolutionary Computation, CEC 2010. Barcelona, Spain, 2010, pp. 1–7.


persistent fname f_bias
total_func_no = 20;
MINMAN=1;      % Maximization

if func_num == 1	fname = str2func('five_uneven_peak_trap');
elseif func_num == 2	fname = str2func('equal_maxima');
elseif func_num == 3	fname = str2func('uneven_decreasing_maxima');
elseif func_num == 4	fname = str2func('himmelblau');
elseif func_num == 5	fname = str2func('six_hump_camel_back');
elseif func_num == 6	fname = str2func('shubert');
elseif func_num == 7	fname = str2func('vincent');
elseif func_num == 8	fname = str2func('shubert');
elseif func_num == 9	fname = str2func('vincent');
elseif func_num == 10	fname = str2func('modified_rastrigin_all');
elseif func_num == 11	fname = str2func('CF1');
elseif func_num == 12	fname = str2func('CF2');
elseif func_num == 13	fname = str2func('CF3');
elseif func_num == 14	fname = str2func('CF3');
elseif func_num == 15	fname = str2func('CF4');
elseif func_num == 16	fname = str2func('CF3');
elseif func_num == 17	fname = str2func('CF4');
elseif func_num == 18	fname = str2func('CF3');
elseif func_num == 19	fname = str2func('CF4');
elseif func_num == 20	fname = str2func('CF4');
else
	fprintf('ERROR: Wrong function number: (%d).\n', func_num);
	fprintf('       Please provide a function number in {1,2,...,%d}\n', total_func_no);
	fprintf('       For now function number == 1\n');
	fname = str2func('five_uneven_peak_trap');
end

f_bias = zeros(1,total_func_no);
fit = f_bias(func_num) + MINMAN*feval(fname,x);

%==============================================================================
% F1: Five-Uneven-Peak Trap
%==============================================================================
function fit = five_uneven_peak_trap(x)
% F1: Five-Uneven-Peak Trap
% Variable ranges: x in [0, 30]
% No. of global peaks: 2
% No. of local peaks:  3.

result = -1;
if (x >=0 & x < 2.5)
	result = 80*(2.5-x);
elseif (x >= 2.5 & x < 5.0)
	result = 64*(x-2.5);
elseif (x >= 5.0 & x < 7.5)
	result = 64*(7.5-x);
elseif (x >= 7.5 & x < 12.5)
	result = 28*(x-7.5);
elseif (x >= 12.5 & x < 17.5)
	result = 28*(17.5-x);
elseif (x >= 17.5 & x < 22.5)
	result = 32*(x-17.5);
elseif (x >= 22.5 & x < 27.5)
	result = 32*(27.5-x);
elseif (x >= 27.5 & x <= 30)
	result = 80*(x-27.5);
end
fit = result;

%==============================================================================
% F2: Equal Maxima
%==============================================================================
function fit = equal_maxima(x)
% F2: Equal Maxima
% Variable ranges: x in [0, 1]
% No. of global peaks: 5
% No. of local peaks:  0.

fit = sin( 5*pi*x ).^6;

%==============================================================================
% F3: Uneven Decreasing Maxima
%==============================================================================
function fit = uneven_decreasing_maxima(x)
% F3: Uneven Decreasing Maxima
% Variable ranges: x in [0, 1]
% No. of global peaks: 1
% No. of local peaks:  4.

fit = exp(-2*log(2)*((x-0.08)/0.854)^2)*sin(5*pi*(x^0.75-0.05))^6;

%==============================================================================
% F4: Himmelblau
%==============================================================================
function fit = himmelblau(x)
% F4: Himmelblau
% Variable ranges: x, y in [−6, 6]
% No. of global peaks: 4
% No. of local peaks:  0.

fit = 200 - (x(1).^2 + x(2) - 11).^2 - (x(1) + x(2).^2 - 7).^2;

%==============================================================================
% F5: Six-Hump Camel Back
%==============================================================================
function fit = six_hump_camel_back(x)
% F5: Six-Hump Camel Back
% Variable ranges: x in [−1.9, 1.9]; y in [−1.1, 1.1]
% No. of global peaks: 2
% No. of local peaks:  2.

fit = -((4-2.1*x(1).^2+(x(1).^4)/3).*x(1).^2+x(1).*x(2)+(4*x(2).^2-4).*x(2).^2);

%==============================================================================
% F6: Shubert
%==============================================================================
function fit = shubert(x)
% F6: Shubert
% Variable ranges: x_i in  [−10, 10]^n, i=1,2,...,n
% No. of global peaks: n*3^n
% No. of local peaks: many.

[tmp,D] = size(x);
result = 1;
j = [1:5];
for i=1:D
	result = result * sum(j.*cos((j+1).*x(i)+j));
end
fit = -result;

%==============================================================================
% F7: Vincent
%==============================================================================
function fit = vincent(x)
% F7: Vincent
% Variable range: x_i in [0.25, 10]^n, i=1,2,...,n
% No. of global optima: 6^n
% No. of local optima:  0.

[tmp,D] = size(x);
fit = sum( sin( 10*log(x) ) )/D;

%==============================================================================
% F8: Modified Rastrigin - All Global Optima
%==============================================================================
function fit = modified_rastrigin_all(x)
% Variable ranges: x_i in [0, 1]^n, i=1,2,...,n
% No. of global peaks: \prod_{i=1}^n k_i
% No. of local peaks:  0.
MMP = 0;
[tmp,D] = size(x);
if D == 2
	MMP = [3 4];
elseif D == 8
	MMP = [1 2 1 2 1 3 1 4];
elseif D == 16
	MMP = [1 1 1 2 1 1 1 2 1 1 1 3 1 1 1 4];
end

fit = -sum( 10 + 9*cos( 2*pi*MMP.*x ) );

%==============================================================================
%1. Composition Function 1, n=6
%==============================================================================
function fit = CF1(x)
global initial_flag
persistent func_num func o sigma lambda bias M

[ps,D] = size(x);
func_num = 6;
lb = -5; ub = 5;
if initial_flag==0
	load data/optima.mat % saved the predefined optima
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	initial_flag=1;
	func.f1 = str2func('FGriewank');
	func.f2 = str2func('FGriewank');
	func.f3 = str2func('FWeierstrass');
	func.f4 = str2func('FWeierstrass');
	func.f5 = str2func('FSphere');
	func.f6 = str2func('FSphere');
	bias = zeros(1,func_num);
	sigma = ones(1,func_num);
	lambda = [1; 1; 8; 8; 1/5; 1/5];
	lambda = repmat(lambda,1,D);
	for i = 1:func_num
		eval(['M.M' int2str(i) '= diag(ones(1,D));']);
	end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);

%==============================================================================
%2. Composition Function 2, n=8
%==============================================================================
function fit = CF2(x)
global initial_flag
persistent func_num func o sigma lambda bias M

[ps,D] = size(x);
func_num = 8;
lb = -5; ub = 5;
if initial_flag==0
	initial_flag=1;
	load data/optima.mat % saved the predefined optima
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	func.f1 = str2func('FRastrigin');
	func.f2 = str2func('FRastrigin');
	func.f3 = str2func('FWeierstrass');
	func.f4 = str2func('FWeierstrass');
	func.f5 = str2func('FGriewank');
	func.f6 = str2func('FGriewank');
	func.f7 = str2func('FSphere');
	func.f8 = str2func('FSphere');
	bias = zeros(1,func_num);
	sigma = ones(1,func_num);
	lambda = [1; 1; 10; 10; 1/10; 1/10; 1/7; 1/7];
	lambda = repmat(lambda,1,D);
	for i = 1:func_num
		eval(['M.M' int2str(i) '= diag(ones(1,D));']);
	end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);

%==============================================================================
%3. Composition Function 3, n=6
%==============================================================================
function fit = CF3(x)
global initial_flag
persistent func_num func o sigma lambda bias M

[ps,D] = size(x);
func_num = 6;
lb = -5; ub = 5;
if initial_flag==0
	initial_flag=1;
	load data/optima.mat % saved the predefined optima, a 10*100 matrix;
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	func.f1 = str2func('FEF8F2');
	func.f2 = str2func('FEF8F2');
	func.f3 = str2func('FWeierstrass');
	func.f4 = str2func('FWeierstrass');
	func.f5 = str2func('FGriewank');
	func.f6 = str2func('FGriewank');
	bias = zeros(1,func_num);
	sigma = [1,1,2,2,2,2];
	lambda = [1/4; 1/10; 2; 1; 2; 5];
	lambda = repmat(lambda,1,D);
	c = ones(1,func_num);
	if 	D==2,	load data/CF3_M_D2.mat,
	elseif 	D==3,	load data/CF3_M_D3.mat,
	elseif 	D==5,	load data/CF3_M_D5.mat,
	elseif 	D==10,	load data/CF3_M_D10.mat,
	elseif 	D==20,	load data/CF3_M_D20.mat,
	else 
		for i = 1:func_num
			%A = normrnd(0,1,D,D);
			%eval(['M.M' int2str(i) '= LocalGramSchmidt( A );']);
			eval(['M.M' int2str(i) '= RotMatrixCondition( D,c(i) );']);
		end
	end
end
fit = hybrid_composition_func(x,func_num,func,o,sigma,lambda,bias,M);

%==============================================================================
%4. Composition Function 4, n=8
%==============================================================================
function fit = CF4(x)
global initial_flag
persistent func_num func o sigma lambda bias M

[ps,D] = size(x);
func_num = 8;
lb = -5; ub = 5;
if initial_flag==0
	initial_flag=1;
	load data/optima.mat % saved the predefined optima, a 10*100 matrix;
	if length( o(1,:) ) >= D
		o = o(:,1:D);
	else
		o = lb + (ub - lb) * rand(func_num,D);
	end
	func.f1 = str2func('FRastrigin');
	func.f2 = str2func('FRastrigin');
	func.f3 = str2func('FEF8F2');
	func.f4 = str2func('FEF8F2');
	func.f5 = str2func('FWeierstrass');
	func.f6 = str2func('FWeierstrass');
	func.f7 = str2func('FGriewank');
	func.f8 = str2func('FGriewank');
	bias = zeros(1,func_num);
	sigma = [1,1,1,1,1,2,2,2];
	lambda = [ 4; 1; 4; 1; 1/10; 1/5; 1/10; 1/40];
	lambda = repmat(lambda,1,D);
	c = ones(1,func_num);
	if 	D==2,	load data/CF4_M_D2.mat,
	elseif 	D==3,	load data/CF4_M_D3.mat,
	elseif 	D==5,	load data/CF4_M_D5.mat,
	elseif 	D==10,	load data/CF4_M_D10.mat,
	elseif 	D==20,	load data/CF4_M_D20.mat,
	else 
		for i = 1:func_num
			%A = normrnd(0,1,D,D);
			%eval(['M.M' int2str(i) '= LocalGramSchmidt( A );']);
			eval(['M.M' int2str(i) '= RotMatrixCondition( D,c(i) );']);
		end
	end
end
fit = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M);

%==============================================================================
% Hybrid Composition General Framework
%==============================================================================
function res = hybrid_composition_func(x, func_num, func, o, sigma, lambda, bias, M)
[ps,D] = size(x);
for i = 1:func_num
	oo = repmat( o(i,:), ps, 1 );
	weight(:,i) = exp( -sum( (x-oo).^2, 2 )./2./( D * sigma(i)^2 ) );
end

[tmp,tmpid] = sort(weight,2);
for i = 1:ps
	weight(i,:) = (weight(i,:)==tmp(i,func_num)) .* weight(i,:) + (weight(i,:)~=tmp(i,func_num)) .* (weight(i,:).*(1-tmp(i,func_num).^10));
end
if sum(weight,2) == 0
	weight = weight + 1;
end
weight = weight ./ repmat( sum(weight,2), 1, func_num );
it = 0;
res = 0;
for i = 1:func_num
	oo = repmat(o(i,:),ps,1);
	eval(['f = feval(func.f' int2str(i) ',((x-oo)./repmat(lambda(i,:),ps,1))*M.M' int2str(i) ');']);
	x1 = 5*ones(1,D);
	eval(['f1 = feval(func.f' int2str(i) ',(x1./lambda(i,:))*M.M' int2str(i) ');']);
	fit1 = 2000 .* f ./ f1;
	res = res + weight(:,i) .* ( fit1 + bias(i) );
end
res = -res;

%==============================================================================
% Basic Functions
%==============================================================================

%------------------------------------------------------------------------------
% Sphere Function
%------------------------------------------------------------------------------
function f = FSphere(x)
%Please notice there is no use to rotate a sphere function, with rotation
%here just for a similar structure as other functions and easy programming
[ps,D] = size(x);
f = sum( x.^2, 2);

%------------------------------------------------------------------------------
% Griewank's Function
%------------------------------------------------------------------------------
function f = FGriewank(x)
[ps,D] = size(x);
f = 1;
for i = 1:D
	f = f.*cos( x(:,i)./sqrt(i) );
end
f = sum( x.^2, 2)./4000 - f + 1;

%------------------------------------------------------------------------------
% Rastrigin's Function
%------------------------------------------------------------------------------
function f = FRastrigin(x)
[ps,D] = size(x);
f = sum( x.^2-10.*cos( 2.*pi.*x )+10, 2 );

%------------------------------------------------------------------------------
% Weierstrass Function
%------------------------------------------------------------------------------
function f = FWeierstrass(x)
[ps,D] = size(x);
x = x+0.5;
a = 0.5;
b = 3;
kmax = 20;
c1(1:kmax+1) = a.^(0:kmax);
c2(1:kmax+1) = 2*pi*b.^(0:kmax);
f = 0;
c = -w(0.5,c1,c2);
for i = 1:D
	f = f + w( x(:,i)', c1, c2 );
end
f = f + c*D;

function y = w(x,c1,c2)
y = zeros(length(x),1);
for k = 1:length(x)
	y(k) = sum( c1.*cos( c2.*x(:,k) ) );
end

%------------------------------------------------------------------------------
% FEF8F2 Function
%------------------------------------------------------------------------------
function f = FEF8F2(x)
[ps,D] = size(x);
f = 0;
for i = 1:(D-1)
	f = f + F8F2( x(:,[i,i+1])+1 );     % (1,...,1) is minimum
end
f = f + F8F2( x(:,[D,1]) +1 );           % (1,...,1) is minimum

%------------------------------------------------------------------------------
% F8F2 Function
%------------------------------------------------------------------------------
function f = F8F2(x) 
f2 = 100.*(x(:,1).^2-x(:,2)).^2+(1-x(:,1)).^2; 
f = 1+f2.^2./4000-cos(f2); 

%------------------------------------------------------------------------------
% classical Gram Schmid  %TODO: why not use matlab's internal functions?
%------------------------------------------------------------------------------
function [q,r] = LocalGramSchmidt (A)
% computes the QR factorization of $A$ via
% classical Gram Schmid 

[n,m] = size(A); 
q = A;    
for j=1:m
    for i=1:j-1 
        r(i,j) = q(:,j)'*q(:,i);
    end
    for i=1:j-1   
      q(:,j) = q(:,j) -  r(i,j)*q(:,i);
    end
    t =  norm(q(:,j),2 ) ;
    q(:,j) = q(:,j) / t ;
    r(j,j) = t  ;
end 
 
%------------------------------------------------------------------------------
% Generates a D-dimensional rotation matrix with predifined Condition Number (c)
%------------------------------------------------------------------------------
function M = RotMatrixCondition(D,c)

% A random normal matrix
A = normrnd(0,1,D,D);

% P Orthogonal matrix
P = LocalGramSchmidt(A);

% A random normal matrix
A = normrnd(0,1,D,D);

% Q Orthogonal matrix
Q = LocalGramSchmidt(A);

% Make a Diagonal matrix D with prespecified Condition Number
u = rand(1,D);
D = c .^ ( (u-min(u))./(max(u)-min(u)) );
D = diag(D);

% M rotation matrix with Condition Number c
M = P * D * Q;

%==============================================================================
%============================= End of file ====================================
%==============================================================================
