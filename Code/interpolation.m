%%%%%%%%%%%%%%
%Torque curve%
%%%%%%%%%%%%%%

% Generate Torque Curve using polynomial fitting

% knots
x = [0 260 628 1099 1214 1570];  %[rad/s]
v = [30 500 610 700 450 280]; %[Nm]

%% Define the query points, xq, and interpolate.
x_foo = linspace(0,1570,1570);
p = polyfit(x,v,2);
torque = polyval(p,x_foo);

figure
plot(x_foo,torque)