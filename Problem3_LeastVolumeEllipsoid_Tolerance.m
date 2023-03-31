clc
clear all
% Define parameters
a = 3;
b = 4;
c = 5; 
point_on_surface=500;
point_inside=500;
rng(0)

% Generate points on the surface
u = 2*pi*rand(point_on_surface,1);
v = pi*rand(point_on_surface,1);
x = a * cos(u) .* sin(v);
y = b * sin(u) .* sin(v);
z = c * cos(v);

% Rotate the ellipsoid
P_surface=[x y z]';
R=Rot([1;1;1],45);
P_surface=R*P_surface;

% Generate points inside the ellipsoid
u = 2 * pi * rand(point_inside,1);
v = pi * rand(point_inside,1);
x = [a * rand(point_inside,1) .* cos(u) .* sin(v)];
y = [b * rand(point_inside,1) .* sin(u) .* sin(v)];
z = [c * rand(point_inside,1) .* cos(v)];

% Rotate the ellipsoid
P_inside=[x y z]';
R=Rot([1;1;1],45);
P_inside=R*P_inside;

% MOSEK test
ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='mosek';
acc_mosek_tight=[];
time_mosek_tight=[];
iter_mosek_tight=[];
for i=1:30
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 10^(-i);
    M = sdpvar(3,3);
    F = [M >= 0];
    for j=1:350
        F = [F, [eye(3) M*P_surface(:,j);(M*P_surface(:,j))' 1]>=0];
    end
    output=optimize(F,-log(det(M)),ops);
    M = value(M);
    acc_mosek_tight=[acc_mosek_tight;i -log10(abs(4/3*pi/sqrt(det(M*M))-4/3*pi*a*b*c))];
    iter_mosek_tight=[iter_mosek_tight;i output.solveroutput.res.info.MSK_IINF_INTPNT_ITER];
    time_mosek_tight=[time_mosek_tight;i output.solvertime];
    M = sdpvar(3,3);
end

% SeDuMi test
M = sdpvar(3,3);
ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='sedumi';
acc_sedumi_tight=[];
time_sedumi_tight=[];
iter_sedumi_tight=[];
for i=1:30
    ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 10^(-i);
    M = sdpvar(3,3);
    F = [M >= 0];
    for j=1:350
        F = [F, [eye(3) M*P_surface(:,j);(M*P_surface(:,j))' 1]>=0];
    end
    output=optimize(F,-log(det(M)),ops);
    M = value(M);
    acc_sedumi_tight=[acc_sedumi_tight;i -log10(abs(4/3*pi/sqrt(det(M*M))-4/3*pi*a*b*c))];
    iter_sedumi_tight=[iter_sedumi_tight;i output.solveroutput.info.iter];
    time_sedumi_tight=[time_sedumi_tight;i output.solvertime];
    M = sdpvar(3,3);
end

subplot(2,2,1);
plot(time_mosek_tight(:,1),time_mosek_tight(:,2),'DisplayName','MOSEK');
hold on
plot(time_sedumi_tight(:,1),time_sedumi_tight(:,2),'DisplayName','SeDuMi');
hold on
title('Cost Time')
 
subplot(2,2,2);
plot(iter_mosek_tight(:,1),iter_mosek_tight(:,2),'DisplayName','MOSEK');
hold on
plot(iter_sedumi_tight(:,1),iter_sedumi_tight(:,2),'DisplayName','SeDuMi');
hold on
title('Iteration')

subplot(2,2,[3 4]);
plot(acc_mosek_tight(:,1),acc_mosek_tight(:,2),'DisplayName','MOSEK');
hold on
plot(acc_sedumi_tight(:,1),acc_sedumi_tight(:,2),'DisplayName','SeDuMi');
hold on
title('Accuracy')

function T=Rot(k,theta)
    k=k/norm(k);
    T=[k(1)^2*(1-cosd(theta))+cosd(theta) k(1)*k(2)*(1-cosd(theta))-k(3)*sind(theta) k(1)*k(3)*(1-cosd(theta))+k(2)*sind(theta);
       k(1)*k(2)*(1-cosd(theta))+k(3)*sind(theta) k(2)^2*(1-cosd(theta))+cosd(theta) k(2)*k(3)*(1-cosd(theta))-k(1)*sind(theta);
       k(1)*k(3)*(1-cosd(theta))-k(2)*sind(theta) k(2)*k(3)*(1-cosd(theta))+k(1)*sind(theta) k(3)^2*(1-cosd(theta))+cosd(theta)];
end