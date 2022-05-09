clear; clc;

h_z = load('h_z.mat');
h_z = h_z.out;

h_z_2 = load('h_z_2.mat');
h_z_2 = h_z_2.out;

v = 9;

front_delta=figure('Position',[200 200 1400 700]);

[m,n,p] = size(h_z) ;
[X,Y] = meshgrid(linspace(0,1,p), linspace(0,3,n)) ;

subplot(1,2,1)
h1 = surf(X,Y,squeeze(h_z(1,:,:))) ;
xlabel('$y$', 'interpreter','latex');
ylabel('$x$', 'interpreter','latex');
zlabel('$H_z(x,y,t)$', 'interpreter','latex');
colormap(jet);
set(gca, 'FontSize', 25)
shading interp;
axis tight;
axis vis3d;
axis equal;
caxis([-0.25 0.25]);
zlim([-0.25 0.75])
xlim([-0.5 1.5])
zticks([])
yp = get(gca,'Ylim');
z1 = [ -0.1 0.5 0.5 -0.1];
y1 = [ yp(1) yp(1) yp(2) yp(2)];
x1 = ones(1,numel(y1))*0;
x2 = ones(1,numel(y1))*1;
p1 = patch(x1,y1,z1, 'b');
set(p1,'facealpha',0.2)
set(p1,'edgealpha',0.2)
p2 = patch(x2,y1,z1, 'b');
set(p2,'facealpha',0.2)
set(p2,'edgealpha',0.2)

subplot(1,2,2)
h2 = surf(X,Y,squeeze(h_z_2(1,:,:))) ;
xlabel('$y$', 'interpreter','latex');
ylabel('$x$', 'interpreter','latex');
zlabel('$H_z(x,y,t)$', 'interpreter','latex');
colormap(jet);
set(gca, 'FontSize', 25)
shading interp;
axis tight;
axis vis3d;
axis equal;
caxis([-0.25 0.25]);
zlim([-0.25 0.75])
xlim([-0.5 1.5])
zticks([])
yp = get(gca,'Ylim');
z1 = [ -0.1 0.5 0.5 -0.1];
y1 = [ yp(1) yp(1) yp(2) yp(2)];
x1 = ones(1,numel(y1))*0;
x2 = ones(1,numel(y1))*1;
p1 = patch(x1,y1,z1, 'b');
set(p1,'facealpha',0.2)
set(p1,'edgealpha',0.2)
p2 = patch(x2,y1,z1, 'b');
set(p2,'facealpha',0.2)
set(p2,'edgealpha',0.2)


for i = 1:v:m
    set(h1,'ZData',squeeze(h_z(i,:,:))) ;
    set(h2,'ZData',squeeze(h_z_2(i,:,:))) ;
    sgtitle({[],[strcat('t =',num2str(i),'\Deltat')]},'FontSize',27);
    
    drawnow
end

