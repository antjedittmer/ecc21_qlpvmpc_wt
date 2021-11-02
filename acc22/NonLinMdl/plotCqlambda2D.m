function plotCqlambda2D(Lin_points,wecs,beta,lambda,Cq)
% plotCqlambda2D plots lambda^(-2) C_Q as a functiof lambda for different
% beta
%
% To estimate the aerodynamic Torque Tr, use the equation
%      Tr = 0.5 *rho *pi*R^3 *C_Q *V^2
%         = 0.5 *rho *pi*R^3 *C_Q *(omega^2*R^2/lambda^2)
% =>  2* Tr /(rho *pi * R^5 omega^2) = lambda^(-2) C_Q

% Get C_Q,beta as a function of TSR lambda at all relevant pitch angle
[betaLUT,idxLUT] = unique(Lin_points.Pitch); %Select all pitch angle in eq. points
lambdaOP = Lin_points.RSpeed * wecs.rb ./ Lin_points.V; % Calculate TSR in OP
lambdaOPLUT = lambdaOP(idxLUT);
idxLambda1 = 1; %min(lambdaOP) 
idxLambda = (lambda >= idxLambda1 & lambda <= max(lambdaOP)); %index TSR

lambdaLUT = lambda(idxLambda); %select all relevant TSR with index TSR
CqPitch = interp2(lambda,beta,Cq,lambdaLUT,betaLUT);
lambdaInvMat = repmat((lambdaLUT .^ (-2))', size(CqPitch,1),1);
lambdaCqMat = lambdaInvMat .* CqPitch;

lambdaCqFun = nan(length(betaLUT), length(lambdaLUT)); %init lambda^-2 C_Q,beta
CqPitch = nan(length(betaLUT),sum(idxLambda)); % initialize C_Q,beta
pitchCt = beta;
for idx = 1: length(betaLUT)
    diffPitch = abs(pitchCt-betaLUT(idx)); % difference to LUT pitch vector
    idxNearestPitch =  min(diffPitch) == diffPitch; % index of nearest pitch
    CqPitch(idx,:) = Cq(idxNearestPitch,idxLambda); % select C_Q,beta
    lambdaCqFun(idx,:) = (lambdaLUT .^ (-2)) .* CqPitch(idx,:)';
end

%% Plot tip speed ratio lambda vs. function lambda^-2 C_Q,beta
cl = colormap('lines');
numCl = numel(unique(cl(:,1)));
numPitch = length(betaLUT);
noDiffLineStyles = ceil(numPitch/numCl);
diffLineStyles = {'-','--','-.'};

figure;
for idx = 1: noDiffLineStyles
    nvecAll = (1:numCl) + (idx -1)*numCl;
    nvec = nvecAll(nvecAll<= numPitch);
    plot(lambdaLUT,lambdaCqFun(nvec,:),'LineStyle',diffLineStyles{idx});
    hold on;
end

% Check monotonocity
diffLambdaCqFun = diff(lambdaCqFun);

xlabel('Tip Speed Ratio \lambda [-]');
ylabel('\lambda^{-2} C_{Q,\beta}(\lambda)[-]');

% for operating points
vecV = [4,11,25];
mCell = {'o','d','s'};
lambdaOPPlot = nan(size(vecV ));
idxSpacing = 1;

for idx = 1: length(vecV)%: idxSpacing: length(lambdaOPLUT)
    lamdaVecPlot = lambdaOP(Lin_points.V == vecV(idx));
    betaVecPlot = Lin_points.Pitch(Lin_points.V == vecV(idx));
    
    idxVec = abs(lambdaLUT - lamdaVecPlot);
    idxOP = idxVec  == min(idxVec);
    
    idxVecBeta = abs(betaLUT - betaVecPlot);
    idxB = idxVecBeta == min(idxVecBeta);
    
    lambdaOPPlot(idx) = lambdaLUT(idxOP);
    
    plot(lambdaLUT(idxOP),lambdaCqFun(idxB,idxOP),mCell{idx},'MarkerFaceColor',cl(idxB,:),'MarkerEdgeColor','k');
end


cellLeg = regexp(sprintf('Pitch: %2.2f [deg],',...
    betaLUT(1:idxSpacing:end)*180/pi),',','split');
temp = [vecV', lambdaOPPlot']';

cellLeg1 = regexp(sprintf('OP %d m/s: %2.2f,',...
     temp(:)),',','split');

X = [1, min(lambdaOP), min(lambdaOP),1];
Y = [0,0, max(lambdaCqMat(:)),  max(lambdaCqMat(:))];
fill(X,Y,'k','FaceAlpha',0.2,'LineStyle','none');
axis tight; grid on
legend([cellLeg(1:end-1),cellLeg1(1:end-1)]);


