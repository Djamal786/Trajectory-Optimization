     %% Generating an optimized trajectory  
    
    nbDOF = 7; %Degrees of freedom of the robot
    name = "iiwa 14"; %name of the robot 
    nbTrajSamples = 100; %Trajectory Samples 
    t_f = 10; %duration of the trajectory, end time [s] 
    t_i = 0; %start time [s]
    caractFreq = 1/(t_f-t_i); %frequency 
    caractPuls = 2*pi*caractFreq;
    n_f = 20; %Fourier parameter
    nbVars = 2*n_f*nbDOF; %280
    trajectoryParameters = zeros(nbVars,1); %Important for optimization with fmincon [280,2]
    paramVectorSize = 55;
    nbCtrlSamples = 10000;
    Geometry =   [0, 0 , 0, pi/2;  %DH-Parameters of the robot [q,d,a,alpha]
                  0, 0, 0, -pi/2
                  0, 0.42, 0, -pi/2
                  0, 0, 0, pi/2
                  0, 0.4, 0, pi/2
                  0, 0, 0, -pi/2
                  0, 0, 0, 0];
            
    Q_L = [80,60,80,-40,-170,-120,-175]' * pi/180; %Constraints of the joints
    Q_U = [120,40,120,40,170,120,175]' * pi/180; 
    Qp_L = -0.5*[85,85,100,75,130,135,135]' * pi/180;
    Qp_U = 0.5*[85,85,100,75,130,135,135]'* pi/180;
    Qpp_U = [10,12,11,10,13,23,22]';
    Qpp_L = -[12,16,12,11,12,21,22]';
    Tau_L = 100*[5,60,23,49,3,9,0.3]' * (-1);
    Tau_U = 100*[5,60,18,47,6,10,0.4]';
    
    augmentedState_max = [Qpp_U; Qp_U; Q_U];
    augmentedState_min = [Qpp_L; Qp_L; Q_L];
    Q0 = [90;50;100;0;0;0;0]*pi/180; %inital condition of the joint angles [rad]
    
       
    %% Setting constraints on initial position and velocity:   
    
    Aeq = zeros(3*nbDOF, nbVars);
    beq = zeros(3*nbDOF,1);
        
    for i=1:nbDOF
        for j=1:n_f
            Aeq(i, 2*n_f*(i-1)+n_f+j)=1/(caractPuls*j); % -> set the initial position to the desired position qi0
            beq(i)=0;
            Aeq(nbDOF+i, 2*n_f*(i-1)+j)=1; % -> set the initial velocities to 0
            Aeq(2*nbDOF+i, 2*n_f*(i-1)+n_f+j)=caractPuls*j; % -> set the initial acceleration to 0
        end
    end
    %% Trajectory optimization 
    
    max_traj_it = 10 %iteration parameter for fmincon
    fval_final = trajectoryCriterion_djamal(nbDOF, paramVectorSize, Q0, t_f, t_i, n_f, nbTrajSamples, caractPuls, trajectoryParameters);
    iteration = 0;

    alg = 'fmin'; % 'fmin' or 'ga'
    
    fun = @(trajectoryParameters_tmp)trajectoryCriterion_djamal(nbDOF, paramVectorSize, Q0, t_f, t_i, n_f, nbTrajSamples, caractPuls, trajectoryParameters_tmp);
    nonlcon = @(trajectoryParameters_tmp) trajectoryConstraints_djamal(nbDOF, name, paramVectorSize, Q0, Geometry, t_f, t_i, n_f, nbTrajSamples, caractPuls, trajectoryParameters_tmp, augmentedState_max, augmentedState_min);
                                           
    [C0,Ceq0]=trajectoryConstraints_djamal(nbDOF, name, paramVectorSize, Q0, Geometry, t_f, t_i, n_f, nbTrajSamples, caractPuls, 0.001*(-2 + 4*rand(nbVars, 1)), augmentedState_max, augmentedState_min)
    
    if ~all(C0 <= 0) % Trajectory generation constraints are violated by the initial point...
        %error('Trajectory generation constraints are violated by the initial point...');
    end
    
    while iteration<max_traj_it && fval_final > 5
        if strcmp(alg,'fmin') % 'interior-point' (default), 'trust-region-reflective', sqp', 'sqp-legacy', 'active-set'. 
            %options = optimoptions('fmincon','Algorithm' ,'interior-point', 'Display','iter', 'MaxFunctionEvaluations', 10000*robot.nbDOF, 'UseParallel', true);
            %leboutet status: options = optimoptions('fmincon','Algorithm' ,'interior-point', 'Display','iter', 'MaxFunctionEvaluations', 10, 'UseParallel', true);
            options = optimoptions('fmincon','Algorithm' ,'interior-point', 'Display','iter','MaxFunctionEvaluations',1e4, 'MaxIterations',1e3, 'UseParallel', true);
            [trajectoryParameters_optim,fval,exitFlag] = fmincon(fun, 0.1*(-2 + 4*rand(nbVars, 1)), [], [], Aeq, beq, [], [], nonlcon, options);
        elseif strcmp(alg,'ga')
            population = 240;
            %leboutet status: options = optimoptions('ga','MaxStallGenerations',20,'MutationFcn',@mutationadaptfeasible, 'FunctionTolerance',1e-10,'Display','iter','MaxGenerations',50*robot.nbDOF,'PlotFcn',{@gaplotbestindiv,@gaplotdistance,@gaplotrange}, 'PopulationSize', population,'UseParallel', true, 'EliteCount', floor(population/10), 'CrossoverFraction', 0.9);
            %leboutet status: [trajectoryParameters_optim,fval] = ga(fun, nbVars, [], [], Aeq, beq, [], [], nonlcon, options);
            [trajectoryParameters_optim,fval] = ga(fun, nbVars, [], [], Aeq, beq, [], [], nonlcon);
        else
            error("Trajectory generation: unknown algorithm");
        end
        
        if fval < fval_final
            trajectoryParameters = trajectoryParameters_optim;
            fval_final = fval;
        end
        
        iteration = iteration + 1;
    end
    
    fprintf('Optimization Complete!');
    fprintf('Lowest value of the cost functon: %d\n', fval_final);
    
    % Trajectory data generation:
    
    fprintf('Generating Trajectory Data...\n');
    trajectoryData = generateTrajectoryData(t_i, t_f, nbCtrlSamples, Q0, caractPuls, n_f, nbDOF, trajectoryParameters);
    fprintf('Trajectory Data Generation Completed!\n');
 

    %% Test the cost function 
%    J = trajectoryConstraints()

    %% Generation of trajectory data, extra function 
    function [trajectoryData] = generateTrajectoryData(t_i, t_f, nbCtrlSamples, Q0, caractPuls, n_f, nbDOF, trajectoryParameters)

    t = linspace(t_i, t_f, nbCtrlSamples);
    trajectoryData.t = t;
    trajectoryData.Q=zeros(nbDOF,nbCtrlSamples);
    trajectoryData.Qp=zeros(nbDOF,nbCtrlSamples);
    trajectoryData.Qpp=zeros(nbDOF,nbCtrlSamples);
    
    for i = 1:nbCtrlSamples
        augmentedState=trajectoryGeneratorFourier_djamal(t(i), caractPuls, n_f, nbDOF, Q0, trajectoryParameters);
        trajectoryData.Qpp(:,i) = augmentedState(1:nbDOF);
        trajectoryData.Qp(:,i) = augmentedState(nbDOF+1:2*nbDOF);
        trajectoryData.Q(:,i) = augmentedState(2*nbDOF+1:end);
    end
    
    end

    








