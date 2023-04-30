function [augmentedState] = trajectoryGeneratorFourier_djamal(t, caractPuls, n_f, nbDOF, Q0, trajectoryParameters)

% Authors: Julien Roux, Quentin Leboutet, Alexandre Janot, Gordon Cheng
%
% Generate a desired joint space trajectory by parametrizing Fourier series
Q = Q0;
Qp = zeros(nbDOF, 1);
Qpp = zeros(nbDOF, 1);
%A = [20,10,20,40,170,120,175] * 10* pi/180; 
if (nbDOF ~= numel(Q))
    error('The Fourier trajectory generation tool detected inconsistancy in nbDOF !')
end

for i=1:nbDOF
    for j=1:n_f
        Aij = trajectoryParameters(2*n_f*(i-1)+j);
        Bij = trajectoryParameters(2*n_f*(i-1)+n_f+j);
        
        Cojt = cos(caractPuls*j*t);
        Sojt = sin(caractPuls*j*t);
        
        Q(i) = Q(i) + Aij/(caractPuls*j)*Sojt - Bij/(caractPuls*j)*Cojt;
        Qp(i) = Qp(i) + Aij*Cojt + Bij*Sojt;
        Qpp(i) = Qpp(i) + Bij*j*Cojt - Aij*j*Sojt;
    end
    Qpp(i) = caractPuls*Qpp(i);

end

augmentedState = [Qpp; Qp; Q];


end
