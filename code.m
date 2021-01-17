%% Constants and System Parameters
hbar = 1.0545718*10^-34;
energy = 1000*1.6*10^-22;
lambda = 10^-3;
omega = energy*lambda/hbar;

%% Solving the coupled differential equation
tbegin = 0;
tend = 100*6.62607004*10^(-34)/(lambda*energy); %100 cycles of oscillations
steps = 1000;
initial = [1; 0];

tspan = linspace(tbegin, tend, steps);
[t, B] = ode45(@solver, tspan, initial);

%% Plotting the result
B = B';
ground = B(1, :);
first = B(2, :);


prob_ground = conj(ground).*ground;
prob_first = conj(first).*first

figure();
subplot(2,1,1);
plot(t,prob_ground);
xlabel('time (s)');
ylabel('Ground State Population');
subplot(2,1,2);
plot(t,prob_first);
xlabel('time (s)');
ylabel('First Excited State Population');

%% Defining the equation
function B_dot = solver(t, B)

hbar = 1.0545718*10^-34;
energy = 1000*1.6*10^-22;
lambda = 10^-3;
omega = energy*lambda/hbar;

B_dot = zeros(2, 1);

% Comment out any one of the two pair of coupled differential equation to
% obtain the solution

% Differential equation after substitution and simplification
B_dot(1) = 3i*energy*(1-2*lambda*sin(omega*t))*B(1)/(hbar) - 4*lambda*omega*B(2)*cos(omega*t)/(3*(1 + lambda*sin(omega*t)));
B_dot(2) = -3i*energy*(1-2*lambda*sin(omega*t))*B(2)/(hbar) + 4*lambda*omega*B(1)*cos(omega*t)/(3*(1 + lambda*sin(omega*t)));

% Original differential equation 
% B_dot(1) = 4*lambda*omega*B(2)*cos(omega*t)*exp(-3i*energy*(t - 2*lambda*(1-cos(omega*t))/omega)/hbar)/(3*(1 + lambda*sin(omega*t)));
% B_dot(2) = -4*lambda*omega*B(1)*cos(omega*t)*exp(3i*energy*(t - 2*lambda*(1-cos(omega*t))/omega)/hbar)/(3*(1 + lambda*sin(omega*t)));
end