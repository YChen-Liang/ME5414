rng(0);
A = [-2 0;0 -3];
B = eye(2);
Ci = 10*rand(2,800);

ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='mosek';
ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1.0e-5;
iter_mosek = [];
time_mosek = [];
for i=1:3:400
    P = sdpvar(2,2);
    Lmix = sdpvar(2*i,2);
    Cmix = Ci(:,1:2*i);
    F = [P>=eye(2)];
    F = [F [P A*P+B*Cmix*Lmix;(A*P+B*Cmix*Lmix)' P]>=eye(4)];
    output=optimize(F,0,ops);
    iter_mosek = [iter_mosek;i output.solveroutput.res.info.MSK_IINF_INTPNT_ITER];
    time_mosek = [time_mosek;i output.solvertime];
end

ops=sdpsettings();
ops.savesolveroutput=1;
ops.solver='sedumi';
ops.sedumi.eps=1.0e-5;
iter_sedumi = [];
time_sedumi = [];
for i=1:3:400
    P = sdpvar(2,2);
    Lmix = sdpvar(2*i,2);
    Cmix = Ci(:,1:2*i);
    F = [P>=eye(2)];
    F = [F [P A*P+B*Cmix*Lmix;(A*P+B*Cmix*Lmix)' P]>=eye(4)];
    output=optimize(F,0,ops);
    iter_sedumi = [iter_sedumi;i output.solveroutput.info.iter];
    time_sedumi = [time_sedumi;i output.solvertime];
end

subplot(1,2,1);
plot(time_mosek(:,1),time_mosek(:,2),'DisplayName','MOSEK');
hold on
plot(time_sedumi(:,1),time_sedumi(:,2),'DisplayName','SeDuMi');
title('Cost Time')
 
subplot(1,2,2);
plot(iter_mosek(:,1),iter_mosek(:,2),'DisplayName','MOSEK');
hold on
plot(iter_sedumi(:,1),iter_sedumi(:,2),'DisplayName','SeDuMi');
title('Iteration')