syms mt N mb Jr Jg Ng Bt Bb Bs Kt Kb rb Ks
M =[mt+N*mb N*mb*rb 0 0;
    N*mb*rb N*mb*rb^2 0 0;
    0 0 Jr 0;
    0 0 0 Jg*Ng^2];

Ce = [Bt 0 0 0;
    0 N*Bb*rb^2 0 0;
    0 0 Bs -Bs;
    0 0 -Bs Bs];
K = [Kt 0 0;
    0 N*Kb*rb^2 0;
    0 0 Ks;
    0 0 -Ks];

Q = [1 0 0 ;
    rb 0 0 ;
     0 1 0;
     0 0 -Ng];
L = [eye(3), [0;0;-1]];
A1 = [zeros(3), L;...
    -M\K, -M\Ce];
B1=[zeros(3); M\Q];