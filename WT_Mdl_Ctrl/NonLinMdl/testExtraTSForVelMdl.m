%% Create cont first order model
TS = 0.01;
tau = 0.5;
atest_c = - 1/tau;
btest_c = - atest_c;
ctest_c = 1;
dtest_c = 0;
sstest_c = ss(atest_c,btest_c,ctest_c,dtest_c);

%% Create velocity second order model
sstestVel_c = ss([0,ctest_c; 0, atest_c], [0;btest_c], [1,0],0);
sstestVel_d = c2d(sstestVel_c,TS);
sstestVelTS_c = ss([0,ctest_c; 0, atest_c], 1/TS * [0;btest_c], [1,0],0);

%% Create test input
tTest = (0:TS:10)';
utest = ones(size(tTest));
nSamples = length(tTest);
mSample = round(nSamples/2);
utest(mSample:end) = 2;
udotTest = [0;diff(utest)];

utest1 = utest;
%utest1(mSample:end) = utest1(mSample:end) + rand(size(utest1(mSample:end)));
utest1(mSample:end) = utest1(mSample:end) + rand(size(utest1(mSample:end)));
udotTest1 = [0;diff(utest1)];

%% Run lsim
figure;
subplot(1,3,1)
lsim(sstest_c,utest,tTest,1)
subplot(1,3,2)
y = lsim(sstestVel_d,udotTest,tTest,[1;0]);
plot(tTest,y); %,tTestudotTest,'Color',[0.5,0.5,0.5]);
title('Test Same B','FontSize',12);
ylabel('Time (s)')
subplot(1,3,3)
lsim(sstestVelTS_c,udotTest,tTest,[1;0]);

figure;
subplot(1,3,1)
lsim(sstest_c,utest1,tTest,1)
subplot(1,3,2)
y = lsim(sstestVel_c,udotTest1,tTest,[1;0]);
plot(tTest,y); %,tTestudotTest,'Color',[0.5,0.5,0.5]);
title('Test Same B','FontSize',12);
ylabel('Time (s)')
subplot(1,3,3)
lsim(sstestVelTS_c,udotTest1,tTest,[1;0]);

% figure;
% subplot(1,3,1)
% lsim(sstest_c,utest1,tTest,1)
% subplot(1,3,2)
% y = lsim(sstestVel_c,udotTest1,tTest,[1;1]);
% plot(tTest,y); %,tTestudotTest,'Color',[0.5,0.5,0.5]);
% title('Test Same B','FontSize',12);
% ylabel('Time (s)')
% subplot(1,3,3)
% lsim(sstestVelTS_c,udotTest1,tTest,[1;1]);

%% Test switch from input to state
DT = 0.0080;
TS = 0.01;
tau = TS;
atest_c = - 1/tau;
btest_c = - atest_c;
ctest_c = 1;
dtest_c = 0;
sstest_c = ss(atest_c,btest_c,ctest_c,dtest_c);
sstest_c_d = c2d(sstest_c, DT);

Atest_c = [0, 1/DT*ctest_c; 0,atest_c]; 
Btest_c = [0; btest_c];
Ctest_c = [1,0]; 
sstest2_c = ss(Atest_c,Btest_c,Ctest_c, 0);
sstest2_c_d = c2d(sstest2_c, DT);

sstest3_c_d = ss([1, sstest_c_d.A; 0,sstest_c_d.A],[sstest_c_d.B; sstest_c_d.B], Ctest_c, 0,DT);
%sstest3_c_d = ss([1, ctest_c; 0,sstest_c_d.A],[0; sstest_c_d.B], Ctest_c, 0,DT);

aExplicit = exp(atest_c*DT);
bExplicit = 1/atest_c * (exp(atest_c*DT) -1) *btest_c;
fprintf('A discrete %2.3f,  exp(atest_c*DT) %2.3f  A(end,end,1) %2.3f \n', sstest_c_d.A, aExplicit, A(end,end,1))
fprintf('B discrete %2.3f,  1/atest(exp(atest_c*DT)-1) %2.3f  B(end,end,1) %2.3f \n', sstest_c_d.B, bExplicit, B(end,end,1)*DT)

load('OutDataTestOL_LPV_windChange', 'xDot06','xkLPV','xk06', 'u','uDot','A','B','wind');
Tref = [5087; repmat(6516,20,1)];

xTest = nan(size(Tref,1),2);
xTest(1,:) = xk06(1,[2,6]);

xTest2 = nan(size(Tref,1),2);
xTest2(1,:) = xk06(1,[2,6]);

xTest3 = nan(size(Tref,1),2);
xTest3(1,:) = xk06(1,[2,6]);

xTestDiff = nan(size(Tref,1),1);
xTestDiff(1) = xk06(1,6);
xTestDiff(2) = xk06(2,6);

TrefDif = [ 0; Tref(2:end) - Tref(1:end-1)];

xTest6States = nan(size(Tref,1),6);
xTest6States(1,:) = xkLPV(1,:);


xTest6StatesLPV = nan(size(Tref,1),6);
xTest6StatesLPV(1,:) = xk06(1,:);


figure;
plot(xk06(1:10,4)')
hold on;
plot(xkLPV(1:10,4)')


for idx = 1: length(Tref)
    xTest(idx+1,1) = sstest_c_d.A *xTest(idx,1) +  sstest_c_d.B * Tref(idx,1);
    xTest(idx+1,2) = xTest(idx +1,1) - xTest(idx,1);
    xTestDiff(idx+1) = sstest_c_d.A * xTestDiff(idx,1) +  sstest_c_d.B * TrefDif(idx);
    
    xTest2(idx+1,:) = sstest2_c_d.A * xTest2(idx,:)' +  sstest2_c_d.B * TrefDif(idx);
    
    xTest3(idx+1,:) = sstest3_c_d.A * xTest3(idx,:)' +  sstest_c_d.B * TrefDif(idx);
    
    xTest6States(idx+1,:) = A(:,:,idx) * xTest6States(idx,:)' +   B(:,:,idx) * uDot(idx,:)'; 
    xTest6StatesLPV(idx+1,:) = A(:,:,idx) * xTest6StatesLPV(idx,:)' +   B(:,:,idx) * uDot(idx,:)'; 
end

tau = 0.1;
Atest_c = [atest_c, -btest_c; 0,-1/tau];
Btest_c = [0; 1/tau];
Ctest_c = [1,0];
sstest2_c = ss(Atest_c,Btest_c,Ctest_c,dtest_c);
figure;
step(sstest_c); hold on; step(sstest2_c )


diffT_out = xk06(2:end,2) -  xk06(1:end-1,2);
diffT_state = xk06(2:end,6);
max(abs(diffT_out - diffT_state))
diffOmegaG_out = xk06(2:end,1) -  xk06(1:end-1,1);
diffOmegaG_state = xk06(2:end,5)*97;
max(abs(diffOmegaG_out - diffOmegaG_state))

figure; 
subplot(2,1,1) 
plot(diffT_out); hold on; plot(xk06(2:end,6),'--')
subplot(2,1,2) 
plot(diffOmegaG_out); hold on; plot(diffOmegaG_state,'--')

diffT_out = xkLPV(2:end,2) -  xkLPV(1:end-1,2);
diffT_state = xkLPV(2:end,6);
max(abs(diffT_out - diffT_state))
diffOmegaG_out = xkLPV(2:end,1) -  xkLPV(1:end-1,1);
diffOmegaG_state = xkLPV(2:end,5)*97;
max(abs(diffOmegaG_out - diffOmegaG_state))

figure; 
subplot(2,1,1) 
plot(diffT_out); hold on; plot( xkLPV(2:end,6),'--'); hold on; plot(xk06(2:end,6),'k:');
subplot(2,1,2) 
plot(diffOmegaG_out); hold on;  plot(97*xkLPV(2:end,5),'--'); hold on; plot(97*xk06(2:end,5),'k:');

%% Get the matrices for the non-linear system

% Aeval = [[0,1,-1,0] ;
%     [-Ks/Jr, -Bs /Jr, Bs/Jr,0];
%     Ks/Jgr, Bs/Jgr, -Bs/Jgr, - Ng/Jgr; %0.37*Ng/Jgr;
%     zeros(1,3), -1/tau_g];


vecOmegaR = squeeze(A(4,3:6,:))';
ThetaS = q(:,1) - q(:,2);
omegaR = qdot(:,1);
omegaGR = qdot(:,2);
T_G = xk06(:,2);

statesNonLin = [ThetaS,omegaR,omegaGR,T_G]';
deltaStatesNonLin = xk06(:,3:6)';

deltaOmegaR = xk06(:,4);

lenSim = length(xk06);
omegaR_nonLin_test = NaN(lenSim,1);
omegaR_nonLin_test(1) = omegaR(1);

% deltaOmegaR_LPV_test = NaN(lenSim,1);
% deltaOmegaR_LPV_test(1) = deltaOmegaR(1);

deltaStates_LPV_test = NaN(size(deltaStatesNonLin));
deltaStates_LPV_test(:,1) = deltaStatesNonLin(:,1);

C_testR = [0,1,0,0];
ny = 2;

for idx = 1: (lenSim-1) 
    [Aeval,Beval,AwithTG,BwithTG,dTrdomega,Tr,TrDiff] = ABqL_WT_KMC(wind(idx),Rotor_Lambda,Rotor_Pitch,Rotor_cQ);
    
    ssC_testR = ss(AwithTG,BwithTG,C_testR, zeros(1,2));
    ssD_testR = c2d(ssC_testR,DT);
    
    Ak = ssD_testR.A;
    Bk = ssD_testR.B;
    
%     A = [eye(ny) diag([97,1])* Ak(3:4,:);
%     zeros(size(Ak,1),ny), Ak];
    %B = [diag([97,1])* Bk(3:4,:); Bk];
    
    % vecOmegaR(idx,:) * statesNonLin(:,idx) + ssD_testR.B(2,:)* [T_r(idx); T_G(idx)]
    omegaR_nonLin_test(idx+1) = ssD_testR.A(2,:) * statesNonLin(:,idx) + ssD_testR.B(2,:)* [T_r(idx); T_G(idx)];
    
    deltaStates_LPV_test(:,idx+1) = ssD_testR.A * deltaStates_LPV_test(:,idx) + ssD_testR.B * uDot(idx,:)' ;
    deltaStates_LPV_test(2,idx+1) = deltaStates_LPV_test(2,idx+1) + dTrdomega * ssD_testR.B(2,1) * deltaStates_LPV_test(2,idx);
    
end

%[omegaR_nonLin_test(2:end) - omegaR_nonLin_test(1:end-1), omegaR(2:end) - omegaR(1:end-1), deltaOmegaR(2:end)]
[deltaStates_LPV_test(2,:)', deltaOmegaR]

deltaOmegaR_NonLin = xk06(:,4);
deltaOmegaR_LPV = xkLPV(:,4);

uNonLin = [T_r, u(:,2)];
plot(T_r); hold on; plot(u(:,2) *97,'--');




endMat = 6;
startMat = 1;
omegaRdotVec_kplus1 = xk06(startMat+1:endMat,4);
ThetaOmegaROmegaGrTg_k =  xk06(startMat:endMat-1,3:6);
vecNonLin = ThetaOmegaROmegaGrTg_k \ (omegaRdotVec_kplus1);

omegaRdotLPV_kplus1 = xkLPV(startMat+1:endMat,4);
ThetaOmegaROmegaGrTgLPV_k =  xkLPV(1:endMat-1,3:6);
vecLin = ThetaOmegaROmegaGrTggLPV_k \ (omegaRdotLPV_kplus1);

% vec sys matrix: Ks/Jgr, Bs/Jgr, -Bs/Jgr, Ng/Jgr;
rowOmegaR = A(4,:,endMat);

fprintf('Ks/Jr:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaR(3), vecNonLin(1), vecLin(1))
fprintf('Bs:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaR(4), vecNonLin(2), vecLin(2))
fprintf('-Bs:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaR(5), vecNonLin(3), vecLin(3))
fprintf('Ng:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n', rowOmegaR(6), vecNonLin(4), vecLin(4))


omegaRdotVec_kplus1 = xkLPV(startMat+1: endMat,5);
ThetaOmegaROmegaGrTg_k =  xk06(startMat:endMat-1,3:6);
vecNonLin = ThetaOmegaROmegaGrTg_k \ (omegaRdotVec_kplus1);

omegaRdotLPV_kplus1 = xkLPV(2:endMat,5);
ThetaOmegaROmegaGrTg_k =  xkLPV(1:endMat-1,3:6);
vecLin = ThetaOmegaROmegaGrTg_k \ (omegaRdotLPV_kplus1);

% vec sys matrix: Ks/Jgr, Bs/Jgr, -Bs/Jgr, Ng/Jgr;
rowOmegaG = A(5,:,endMat);

fprintf('Ks/Jr:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaG(3), vecNonLin(1), vecLin(1))
fprintf('Bs:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaG(4), vecNonLin(2), vecLin(2))
fprintf('-Bs:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n',  rowOmegaG(5), vecNonLin(3), vecLin(3))
fprintf('Ng:A: %2.3e, NonLin  %2.3e, Lin %2.3e \n', rowOmegaG(6), vecNonLin(4), vecLin(4))




xDiffRel = nan(length(xk06) -1 ,6);
xDiffAbs = xDiffRel;

for idx = 0:  (length(xk06) -2)
 xDiffRel(idx+1,:) = [(A(:,:,1+idx)* xk06(1+idx,:)' + B(:,:,1+idx) * uDot(1+idx,:)')' - xk06(2+idx,:)] ./ (xk06(2+idx,:)+eps);
 xDiffAbs(idx+1,:) = [(A(:,:,1+idx)* xk06(1+idx,:)' + B(:,:,1+idx) * uDot(1+idx,:)')' - xk06(2+idx,:)];
end

mat_kplus1 = xk06(2:end,3:6);
matrixNonLin = ThetaOmegaROmegaGrTg_k \ mat_kplus1;





