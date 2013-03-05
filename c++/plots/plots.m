/******************************************************************************
 * Version: 1.0
 * Last modified on: 21 January, 2013 
 * Developers: Michael G. Epitropakis, Xiaodong Li.
 *      email: mge_(AT)_cs_(DOT)_stir_(DOT)_ac_(DOT)_uk 
 *           : xiaodong_(DOT)_li_(AT)_rmit_(DOT)_edu_(DOT)_au 
 * ***************************************************************************/

% F1: Five-Uneven-Peak Trap
figure(1);
data=load('f1.dat');
plot(data(:,1),data(:,2),'b-');
saveas(gcf,'figs/five_uneven_peak_trap.eps','psc2')
saveas(gcf,'figs/five_uneven_peak_trap.png','png')

% F2: Equal Maxima
figure(2);
data=load('f2.dat');
plot(data(:,1),data(:,2),'b-');
saveas(gcf,'figs/equal_maxima.eps','psc2')
saveas(gcf,'figs/equal_maxima.png','png')


% F3: Uneven Decreasing Maxima
figure(3);
data=load('f3.dat');
plot(data(:,1),data(:,2),'b-');
saveas(gcf,'figs/uneven_decreasing_maxima.eps','psc2')
saveas(gcf,'figs/uneven_decreasing_maxima.png','png')

% F4: Himmelblau
figure(4);
[X,Y] = meshgrid(-6:.05:6);
Z=load('f4.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([-6 6 -6 6 -2000 200])
saveas(gcf,'figs/himmelblau.eps','psc2')
saveas(gcf,'figs/himmelblau.png','png')

% F5: Six-Hump Camel Back
figure(5);
[X,Y] = meshgrid([-1.9:0.019:1.9],[-1.1:0.011:1.1]);
Z=load('f5.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([-1.9 1.9 -1.1 1.1 -25 5])
saveas(gcf,'figs/six_hump_camel_back.eps','psc2')
saveas(gcf,'figs/six_hump_camel_back.png','png')


% F6: Shubert
figure(6);
[X,Y] = meshgrid(-10:0.1:10);
Z=load('f6.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([-10 10 -10 10 -300 200])
saveas(gcf,'figs/shubert.eps','psc2')
saveas(gcf,'figs/shubert.png','png')

% F7: Vincent
figure(7);
[X,Y] = meshgrid(0.25:0.05:10);
Z=load('f7.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([0.25 10 0.25 10 -1 1])
saveas(gcf,'figs/vincent.eps','psc2')
saveas(gcf,'figs/vincent.png','png')

% F8: Modified Rastrigin - All Global Optima
figure(8);
[X,Y] = meshgrid(0:0.01:1);
Z=load('f8.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
axis([0 1 0 1 -40 10])
saveas(gcf,'figs/mod_rastrigin_all.eps','psc2')
saveas(gcf,'figs/mod_rastrigin_all.png','png')

% CF1: Composition Function 1
figure(9);
[X,Y] = meshgrid(-5:0.01:5);
Z=load('CF1.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
saveas(gcf,'figs/CF1.eps','psc2')
saveas(gcf,'figs/CF1.png','png')

% CF2: Composition Function 2
figure(10);
[X,Y] = meshgrid(-5:0.01:5);
Z=load('CF2.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
saveas(gcf,'figs/CF2.eps','psc2')
saveas(gcf,'figs/CF2.png','png')

% CF3: Composition Function 3
figure(11);
[X,Y] = meshgrid(-5:0.01:5);
Z=load('CF3.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
saveas(gcf,'figs/CF3.eps','psc2')
saveas(gcf,'figs/CF3.png','png')

% CF4: Composition Function 4
figure(12);
[X,Y] = meshgrid(-5:0.01:5);
Z=load('CF4.dat');
surfc(X,Y,Z, 'FaceColor','interp','FaceLighting','phong','EdgeColor','none');
saveas(gcf,'figs/CF4.eps','psc2')
saveas(gcf,'figs/CF4.png','png')
