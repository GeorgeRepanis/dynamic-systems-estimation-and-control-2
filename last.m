clc; clear;

% Πραγματικές παράμετροι
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
omega = 2;
Tfinal = 20;
Ts = 0.1;

% Δημιουργία "πραγματικών" σημάτων με πολύ μικρό Ts
Ts_fine = 0.001;
t_fine = 0:Ts_fine:Tfinal;
A0_ref = 4;
u_fine = A0_ref * sin(omega * t_fine);

% Πραγματική προσομοίωση
q = zeros(size(t_fine));
dq = zeros(size(t_fine));
for k = 1:length(t_fine)-1
    ddq = (1 / (m * L^2)) * (u_fine(k) - c * dq(k) - m * g * L * q(k));
    dq(k+1) = dq(k) + Ts_fine * ddq;
    q(k+1) = q(k) + Ts_fine * dq(k);
end

% Πραγματικές παράμετροι
theta_true = [m*L^2; c; m*g*L];

% Τιμές A0 προς δοκιμή
A0_values = [1, 2, 4, 6, 8];
errors = zeros(size(A0_values));

% Επαναληπτικός υπολογισμός εκτιμήσεων για κάθε A0
for i = 1:length(A0_values)
    A0 = A0_values(i);
    t = 0:Ts:Tfinal;
    u = A0 * sin(omega * t);

    % Δειγματοληψία
    q_s = interp1(t_fine, q, t);
    dq_s = interp1(t_fine, dq, t);
    ddq_s = diff(dq_s) / Ts;

    q_trim = q_s(1:end-1);
    dq_trim = dq_s(1:end-1);
    u_trim = u(1:end-1);

    Phi = [ddq_s', dq_trim', q_trim'];
    theta_hat = (Phi' * Phi) \ (Phi' * u_trim');
    errors(i) = norm(theta_true - theta_hat);  % Σφάλμα ||θ - θ̂||
end

% ΔΙΑΓΡΑΜΜΑ
figure;
plot(A0_values, errors, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('A_0 (πλάτος εισόδου)');
ylabel('Σφάλμα εκτίμησης ||θ - θ̂||');
title('Επίδραση του πλάτους εισόδου στην εκτίμηση παραμέτρων');
grid on;
