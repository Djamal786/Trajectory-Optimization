function [J] = trajectoryCriterion_Djamal(nbDOF, paramVectorSize, Q0, t_f, t_i, n_f, nbTrajSamples, caractPuls, trajectoryParameters)

W = zeros(nbDOF*nbTrajSamples, paramVectorSize);
timeSamples = linspace(t_i, t_f, nbTrajSamples);
Q = zeros(7,nbTrajSamples); 
Qp = zeros(7,nbTrajSamples);
Qpp = zeros(7,nbTrajSamples);
n = 10*(1+rand(7,1));% * 100;  

%% Trajectory with sinus, no Fourier 
%phi = linspace(0,2*pi,nbTrajSamples);
% for i=1:nbTrajSamples
%     Q = sin_traj(phi(i),n,nbTrajSamples,1);
%     Qp = sin_traj(phi(i),n,nbTrajSamples,2);
%     Qpp = sin_traj(phi(i),n,nbTrajSamples,3);
%     W(robot.nbDOF*(i-1)+1:robot.nbDOF*i,:) =  Y_kukalwr_nom_red_mex(Q,Qp,Qpp); %feval(sprintf('Regressor_Y_%s', robot.name),Q, Qp, Qpp, robot.numericalParameters.Geometry, robot.numericalParameters.Gravity); 
% end

%% Trajectory random, no Fourier 
% for i=1:nbTrajSamples
%     Q = rand(7,1);%-0.5)*pi;
%     Qp = rand(7,1);%-0.5)*4*pi;
%     Qpp = rand(7,1);%-0.5)*8*pi;
%     W(robot.nbDOF*(i-1)+1:robot.nbDOF*i,:) =  Y_kukalwr_nom_red_mex(Q,Qp,Qpp); %feval(sprintf('Regressor_Y_%s', robot.name),Q, Qp, Qpp, robot.numericalParameters.Geometry, robot.numericalParameters.Gravity); 
% end

%% Trajectory with Fourier (Toolbox)

for i=1:nbTrajSamples
    augmentedState = trajectoryGeneratorFourier_djamal(timeSamples(i), caractPuls, n_f, nbDOF, Q0, trajectoryParameters); % augmentedState = [Qpp; Qp; Q], [21,1]
    Qpp = augmentedState(1:nbDOF);
    Qp = augmentedState(nbDOF+1:2*nbDOF);
    Q = augmentedState(2*nbDOF+1:end);
    W(nbDOF*(i-1)+1:nbDOF*i,:) =  Y_kukalwr_nom_red_mex(Q,Qp,Qpp); %feval(sprintf('Regressor_Y_%s', robot.name),Q, Qp, Qpp, robot.numericalParameters.Geometry, robot.numericalParameters.Gravity); 
end

%% Compute the optimization cost:
k1 = 1;
k2 = 100;
S = svd(W);
sig_min = min(S);
C = cond(W'*W);
J = k1*C+k2*1/sig_min;


end

