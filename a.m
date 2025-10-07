clc; clear;

% ΠΡΑΓΜΑΤΙΚΟ ΣΥΣΤΗΜΑ
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
A0 = 4; omega = 2;
Tfinal = 20;
u_func = @(t) A0 * sin(omega * t);
Ts_fine = 0.001;
t_fine = 0:Ts_fine:Tfinal;
u_fine = u_func(t_fine);

% Πραγματική προσομοίωση με Euler
q = zeros(size(t_fine));
dq = zeros(size(t_fine));
for k = 1:length(t_fine)-1
    ddq = (1 / (m * L^2)) * (u_fine(k) - c * dq(k) - m * g * L * q(k));
    dq(k+1) = dq(k) + Ts_fine * ddq;
    q(k+1) = q(k) + Ts_fine * dq(k);
end

% Πλήρης μέτρηση x(t), u(t)
Ts_a = 0.1;
t_a = 0:Ts_a:Tfinal;
q_a = interp1(t_fine, q, t_a);
dq_a = interp1(t_fine, dq, t_a);
u_a = interp1(t_fine, u_fine, t_a);

ddq_a = diff(dq_a)/Ts_a;
q_a = q_a(1:end-1);
dq_a = dq_a(1:end-1);
u_a = u_a(1:end-1);

% Least Squares
Phi_a = [ddq_a', dq_a', q_a'];
theta_hat = (Phi_a' * Phi_a) \ (Phi_a' * u_a');
theta1 = theta_hat(1); theta2 = theta_hat(2); theta3 = theta_hat(3);

% Ανακατασκευή q̂(t)
q_hat_a = zeros(size(q_a));
dq_hat_a = zeros(size(dq_a));
for k = 1:length(q_hat_a)-1
    ddq_hat = (1 / theta1) * (u_a(k) - theta2 * dq_hat_a(k) - theta3 * q_hat_a(k));
    dq_hat_a(k+1) = dq_hat_a(k) + Ts_a * ddq_hat;
    q_hat_a(k+1) = q_hat_a(k) + Ts_a * dq_hat_a(k);
end
eq_a = q_a - q_hat_a;

% Διαγράμματα
figure('Name','Πλήρης μέτρηση x(t), u(t)');
subplot(3,1,1); plot(t_a(1:end-1), q_a, 'b', 'LineWidth', 1.5);
ylabel('q(t)'); title('Πραγματική γωνία q(t)'); grid on;
subplot(3,1,2); plot(t_a(1:end-1), q_hat_a, 'r', 'LineWidth', 1.5);
ylabel('q̂(t)'); title('Εκτιμημένη γωνία q̂(t)'); grid on;
subplot(3,1,3); plot(t_a(1:end-1), eq_a, 'k', 'LineWidth', 1.5);
ylabel('e_q(t)'); xlabel('t [sec]');
title('Σφάλμα e_q(t) = q(t) - q̂(t)'); grid on;
