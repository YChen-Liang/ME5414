clear all
clc
format short e
A=[1/3 0;0 1/7];
ellipsoid_num=400;
Q=zeros(2,2,ellipsoid_num);
c=zeros(2,ellipsoid_num);
angle=360/ellipsoid_num;
for i=0:ellipsoid_num-1
    R=[cosd(i*angle) -sind(i*angle);sind(i*angle) cosd(i*angle)];
    Q(:,:,i+1)=R'*A*R;
    c(:,i+1)=R'*[-1;0];
end

ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='mosek';
ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-5;
acc_mosek=[];
time_mosek=[];
iter_mosek=[];
for i=1:400
    t = sdpvar(1,1);
    x = sdpvar(2,1);
    F = [[t*eye(2) x;x' 1]>=0];
    for j=1:i
        F = [F,[inv(Q(:,:,j)) x-c(:,j);(x-c(:,j))' 1]>=0];
    end
    output=optimize(F,t,ops);
    x=value(x);
    acc_mosek=[acc_mosek;i -log10(norm(x))];
    time_mosek=[time_mosek;i output.solvertime];
    iter_mosek=[iter_mosek;i output.solveroutput.res.info.MSK_IINF_INTPNT_ITER];
end

ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='sedumi';
ops.sedumi.eps=1.0e-5;
acc_sedumi=[];
time_sedumi=[];
iter_sedumi=[];
for i=1:400
    t = sdpvar(1,1);
    x = sdpvar(2,1);
    F = [[t*eye(2) x;x' 1]>=0];
    for j=1:i
        F = [F,[inv(Q(:,:,j)) x-c(:,j);(x-c(:,j))' 1]>=0];
    end
    output=optimize(F,t,ops);
    x=value(x);
    acc_sedumi=[acc_sedumi;i -log10(norm(x))];
    time_sedumi=[time_sedumi;i output.solvertime];
    iter_sedumi=[iter_sedumi;i output.solveroutput.info.iter];
end

subplot(2,2,1);
plot(time_mosek(:,1),time_mosek(:,2),'DisplayName','MOSEK');
hold on
plot(time_sedumi(:,1),time_sedumi(:,2),'DisplayName','SeDuMI');
title('Cost Time')

subplot(2,2,2);
plot(iter_mosek(:,1),iter_mosek(:,2),'DisplayName','MOSEK');
hold on
plot(iter_sedumi(:,1),iter_sedumi(:,2),'DisplayName','SeDuMI');
title('Iteration')

subplot(2,2,[3 4]);
plot(acc_mosek(:,1),acc_mosek(:,2),'DisplayName','MOSEK');
hold on
plot(acc_sedumi(:,1),acc_sedumi(:,2),'DisplayName','SeDuMI');
title('Accuracy')