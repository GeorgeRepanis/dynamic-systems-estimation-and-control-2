% Ορισμός παραμέτρων
m = 0.75;
L = 1.25;
c = 0.15;
g = 9.81;
A0 = 4;
omega = 2;

% Συνάρτηση εξισώσεων κατάστασης
A = [0 1; -g/L -c/(m*L^2)];
B = [0; 1/(m*L^2)];

% Χρόνος προσομοίωσης και βήμα ολοκλήρωσης
tspan = [0 20];
dt = 9e-4; % Βήμα ολοκλήρωσης μικρότερο από 10^(-3)
t = 0:dt:tspan(2);

% Είσοδος u(t) = A0*sin(omega*t)
u = A0 * sin(omega * t);

% Ορισμός της συνάρτησης κατάστασης
odefun = @(t, x) A*x + B*A0*sin(omega*t);

% Αρχικές συνθήκες
x0 = [0; 0];

% Επίλυση με ODE45
options = odeset('RelTol',1e-6,'AbsTol',1e-6); % Αύξηση ακρίβειας
[t, x] = ode45(odefun, tspan, x0, options);

% Γραφικές Παραστάσεις
figure;
subplot(2,1,1);
plot(t, x(:,1), 'b', 'LineWidth', 1.5);
xlabel('Χρόνος (s)'); ylabel('Θέση q (rad)'); grid on;
title('Απόκριση θέσης q(t)');

subplot(2,1,2);
plot(t, x(:,2), 'r', 'LineWidth', 1.5);
xlabel('Χρόνος (s)'); ylabel('Ταχύτητα dq/dt (rad/s)'); grid on;
title('Απόκριση ταχύτητας dq/dt');

