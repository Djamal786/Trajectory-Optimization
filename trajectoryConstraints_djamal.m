function [C, Ceq] = trajectoryConstraints_djamal(nbDOF, name, paramVectorSize, Q0, Geometry, t_f, t_i, n_f, nbTrajSamples, caractPuls, trajectoryParameters, augmentedState_max, augmentedState_min)

augmentedJointState = zeros(3*nbDOF*nbTrajSamples,1);
cartesianState1 = zeros(3*nbTrajSamples,1);
cartesianState2 = zeros(3*nbTrajSamples,1);
timeSamples = linspace(t_i, t_f, nbTrajSamples);
d = [0.42;0.4]; %DH-Parameters of the robot 
for i=1:nbTrajSamples 
    augmentedState = trajectoryGeneratorFourier_djamal(timeSamples(i), caractPuls, n_f, nbDOF, Q0, trajectoryParameters);  % augmentedState = [Qpp; Qp; Q];
    augmentedJointState(3*nbDOF*(i-1)+1:3*nbDOF*i,1) = augmentedState;
    Q = augmentedState(2*nbDOF+1:end);
    HT1 = T_full_7(Q,d); 
    cartesianState1(i,1) = -HT1(3,4)+0.3; 
 end

C = [augmentedJointState-repmat(augmentedState_max-0.1,nbTrajSamples,1);repmat(augmentedState_min+0.1,nbTrajSamples,1)-augmentedJointState; cartesianState1]; % C <= 0
% C = [augmentedJointState-repmat(augmentedState_max-0.1,nbTrajSamples,1);repmat(augmentedState_min+0.1,nbTrajSamples,1)-augmentedJointState]; % C <= 0
Ceq = [];

end
