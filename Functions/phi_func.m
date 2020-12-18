function f = phi_func(y, lure)

% Extract parameters
W1 = lure.W1;
W2 = lure.W2;
b1 = lure.b1;
b2 = lure.b2;

% Evaluate nonlinearity
f = zeros(size(y));
for i = 1:length(y)
    f(i) =  W2.'*tanh(W1*y(i)+b1) + b2;
end

end