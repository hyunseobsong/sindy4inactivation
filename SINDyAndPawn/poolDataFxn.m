function [xout,Theta] = poolDataFxn(xin,polyorder, name, useReciprocal)
%{
Creates the library of candidate terms.

INPUT
xin: The columns containing the dataset's input variable values
polyorder: The maximum order of polynomials that the library will hold
name: The dataset name. Used to determine how to use the symbolic variables 
and build the library
useReciprocal: If useReciprocal is true, then the reciprocal of temperature, 
or 1 divided by temperature, will be used instead of temperature itself 
when pooling the data. Depending on the files we're working with, the 
symbolic variable corresponding to temperature may change. For example, in 
the Cerf dataset, x3 corresponds to temperature, but in Villa-Rojas, x2 
corresponds to temperature.

OUTPUT
xout: Array of symbolic variables representing candidates of the model
equation's coefficients.
Theta: The library of monomial terms created from polynomial combinations.
%}
xout = poolsym(polyorder, name, useReciprocal);          % Symbolic Theta column
xout_eqn = matlabFunction(xout);    % Theta column functions

m = size(xin,1);    % no. of data points

% Building library row-by-row

Theta = zeros(m,length(xout));

for i = 1:m
    % Supply input (must match nVars)
	if name == "Cerf" || name == "Gil"
		Theta(i,:) = xout_eqn(xin(i,1),xin(i,2),xin(i,3)); %Cerf, Gil
	elseif name == "Gaillard" || name == "Villa-Rojas"
		Theta(i,:) = xout_eqn(xin(i,1),xin(i,2));  %Gaillard_1998a, Villa-Rojas
	else
		Theta(i,:) = xout_eqn(xin(i,1),xin(i,2),xin(i,3),xin(i,4)); %Juneja
    end
    
end

end

%%-------------------------------------------------------------------------

function xout = poolsym(polyorder, name, useReciprocal)
% Input variables (must match nVars)

syms xout x1 x2 x3 x4 %Juneja

x = [];

if name == "Cerf"
	if useReciprocal
		%disp("Cerf x = [x1 x2 1/x3]")
		x = [x1 x2 1/x3];
	else
		%disp("Cerf x = [x1 x2 x3]")
		x = [x1 x2 x3];
	end
elseif name == "Gil"
	if useReciprocal
		%disp("Gil x = [1/x1 x2 x3]")
		x = [1/x1 x2 x3];
	else
		%disp("Gil x = [x1 x2 x3]")
		x = [x1 x2 x3];
	end
elseif name == "Gaillard" || name == "Villa-Rojas" 
	if useReciprocal
		%disp("Gaillard or Villa-Rojas x = [x1 1/x2]")
		x = [x1 1/x2];
	else
		%disp("Gaillard or Villa-Rojas x = [x1 x2]")
		x = [x1 x2];
	end
elseif name == "Juneja"
	if useReciprocal
		%disp("Juneja x = [1/x1 x2 x3 x4]")
		x = [1/x1 x2 x3 x4];
	else
		%disp("Juneja x = [x1 x2 x3 x4]")
		x = [x1 x2 x3 x4];
	end
end

%The library is created by multiplying the terms amongst themselves.

ind = 1;                    % Poly order 0
xout(ind) = 1;
ind = ind+1;

for i = 1:length(x)         % Poly order 1
    xout(ind) = x(i);
    ind = ind+1;
end

if polyorder >= 2           % Poly order 2
    for i = 1:length(x)
        for j = i:length(x)
            xout(ind) = x(i)*x(j);
            ind = ind+1;
        end
    end
end

if polyorder >=3            % Poly order 3
    for i = 1:length(x)
        for j = i:length(x)
            for k = i:length(x)
                xout(ind) = x(i)*x(j)*x(k);
                ind = ind+1;
            end
        end
    end
end

if polyorder >=4           % Poly order 4
    for i = 1:length(x)
        for j = i:length(x)
            for k = i:length(x)
                for l = i:length(x)
                    xout(ind) = x(i)*x(j)*x(k)*x(l);
                    ind = ind+1;
                end
            end
        end
    end
end

for i = 1:length(xout)      % Remove possible duplicates
    for j = i+1:length(xout)
        if xout(i)==xout(j)
            xout(j)=0;
        end
    end
end
xout(xout == 0) = [];

end
