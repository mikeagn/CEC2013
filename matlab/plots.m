% ******************************************************************************
% * Version: 1.0
% * Last modified on: 21 January, 2013 
% * Developers: Michael G. Epitropakis, Xiaodong Li.
% *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
% *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
% * ****************************************************************************

% F1: Five-Uneven-Peak Trap
figure(1);
x=[0:0.1:30];
y=zeros(size(x));
for i=1:length(x)
	y(i) = niching_func(x(i),1);
end
plot(x,y,'b-');
saveas(gcf,'figs/F1.eps','psc2')
saveas(gcf,'figs/F1.png','png')

% F2: Equal Maxima
figure(2);
x=[0:0.001:1];
y=zeros(size(x));
for i=1:length(x)
	y(i) = niching_func(x(i),2);
end
plot(x,y,'b-');
saveas(gcf,'figs/F2.eps','psc2')
saveas(gcf,'figs/F2.png','png')


% F3: Uneven Decreasing Maxima
figure(3);
x=[0:0.001:1];
y=zeros(size(x));
for i=1:length(x)
	y(i) = niching_func(x(i),3);
end
plot(x,y,'b-');
saveas(gcf,'figs/F3.eps','psc2')
saveas(gcf,'figs/F3.png','png')

% F4: Himmelblau
figure(4);
[X,Y] = meshgrid(-6:.05:6);
Z=zeros(size(X));
i=1;
for x=-6:0.05:6
	j=1;
	for y=-6:0.05:6
		Z(i,j) = niching_func([x y],4);
		j=j+1;
	end
	i=i+1;
end
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([-6 6 -6 6 -2000 200])
saveas(gcf,'figs/F4.eps','psc2')
saveas(gcf,'figs/F4.png','png')

% F5: Six-Hump Camel Back
figure(5);
[X,Y] = meshgrid(linspace(-1.9,1.9,200),linspace(-1.1,1.1,200));
Z=zeros(size(X));
i=1;
for x=linspace(-1.9,1.9,200)
	j=1;
	for y=linspace(-1.1,1.1,200)
		Z(i,j) = niching_func([x y],5);
		j=j+1;
	end
	i=i+1;
end
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([-1.9 1.9 -1.1 1.1 -6 2])
saveas(gcf,'figs/F5.eps','psc2')
saveas(gcf,'figs/F5.png','png')


% F6: Shubert
figure(6);
[X,Y] = meshgrid(-10:0.1:10);
Z=zeros(size(X));
i=1;
for x=-10:0.1:10
	j=1;
	for y=-10:0.1:10
		Z(i,j) = niching_func([x y],6);
		j=j+1;
	end
	i=i+1;
end
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
%axis([-10 10 -10 10 -300 200])
saveas(gcf,'figs/F6.eps','psc2')
saveas(gcf,'figs/F6.png','png')

% F7: Vincent
figure(7);
[X,Y] = meshgrid(0.25:0.05:10);
Z=zeros(size(X));
i=1;
for x=0.25:0.05:10
	j=1;
	for y=0.25:0.05:10
		Z(i,j) = niching_func([x y],7);
		j=j+1;
	end
	i=i+1;
end
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
%axis([0.25 10 0.25 10 -1 1])
saveas(gcf,'figs/F7.eps','psc2')
saveas(gcf,'figs/F7.png','png')

% F8: Modified Rastrigin - All Global Optima
figure(8);
[X,Y] = meshgrid(0:0.01:1);
Z=zeros(size(X));
i=1;
for x=0:0.01:1
	j=1;
	for y=0:0.01:1
		Z(i,j) = niching_func([x y],10);
		j=j+1;
	end
	i=i+1;
end
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([0 1 0 1 -40 10])
saveas(gcf,'figs/F9.eps','psc2')
saveas(gcf,'figs/F9.png','png')

% F9: Composition function 1
figure(9);
x=-5:0.1:5; y=x;
global initial_flag
initial_flag=0;
func_num = 11;

L = length(x);
f = zeros(L);

for i=1:L
	for j=1:L
		f(i,j) = niching_func([x(i),y(j)],func_num);
	end
end

surfc(x,y,f,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
grid off;
fname=sprintf('figs/CF_%d.eps',100+func_num)
saveas(gcf,fname,'psc2')
fname=sprintf('figs/CF_%d.png',100+func_num)
saveas(gcf,fname,'png')

% F10: Composition function 2
figure(10);
x=-5:0.1:5; y=x;
global initial_flag
initial_flag=0;
func_num = 12;

L = length(x);
f = zeros(L);

for i=1:L
	for j=1:L
		f(i,j) = niching_func([x(i),y(j)],func_num);
	end
end

surfc(x,y,f,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
grid off;
fname=sprintf('figs/CF_%d.eps',100+func_num)
saveas(gcf,fname,'psc2')
fname=sprintf('figs/CF_%d.png',100+func_num)
saveas(gcf,fname,'png')

% F11: Composition function 3
figure(11);
x=-5:0.1:5; y=x;
global initial_flag
initial_flag=0;
func_num = 13;

L = length(x);
f = zeros(L);

for i=1:L
	for j=1:L
		f(i,j) = niching_func([x(i),y(j)],func_num);
	end
end

surfc(x,y,f,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
grid off;
fname=sprintf('figs/CF_%d.eps',100+func_num)
saveas(gcf,fname,'psc2')
fname=sprintf('figs/CF_%d.png',100+func_num)
saveas(gcf,fname,'png')

% F12: Composition function 4
figure(12);
x=-5:0.1:5; y=x;
global initial_flag
initial_flag=0;
func_num = 14;

L = length(x);
f = zeros(L);

for i=1:L
	for j=1:L
		f(i,j) = niching_func([x(i),y(j)],func_num);
	end
end

surfc(x,y,f,'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
grid off;
fname=sprintf('figs/CF_%d.eps',100+func_num)
saveas(gcf,fname,'psc2')
fname=sprintf('figs/CF_%d.png',100+func_num)
saveas(gcf,fname,'png')

