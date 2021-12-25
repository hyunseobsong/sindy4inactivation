function Xi = sparsifyDynamics(Theta,dXdt,lambda)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz

%{
INPUT
Theta: The library of monomial terms created from polynomial combinations.
dXdt: The column containing the observed D values
lambda: The minimum threshold that coefficients must reach for inclusion in
the model equation

OUTPUT
Xi: an array containing the coefficients of the model equation.
%}

% compute Sparse regression: sequential least squares
Xi = Theta\dXdt;  % initial guess: Least-squares

%******
ndXdt=size(dXdt,2);
% lambda is our sparsification knob.
for k=1:10
    idxSmall = (abs(Xi)<lambda);   % find small coefficients
    Xi(idxSmall)=0;                % and threshold
    for ind = 1:ndXdt                   % n is state dimension
        idxBig = ~idxSmall(:,ind);
        % Regress dynamics onto remaining terms to find sparse Xi
        Xi(idxBig,ind) = Theta(:,idxBig)\dXdt(:,ind); 
    end
end

