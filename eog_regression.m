function [corrected_rw,b,X] = eog_regression(raw,eog)
% specific to two bipolar channels
Y = raw';
X = [ones(size(Y)), eog'];
X_inv = pinv(X);
b = X_inv*Y;
corrected_rw = (Y-X*b)';