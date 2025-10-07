clc; clear;

% ΠΡΑΓΜΑΤΙΚΟ ΣΥΣΤΗΜΑ (ίδιο όπως πριν)
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
A0 = 4; omega = 2;
Tfinal = 20;
u_func = @(t) A0 * sin(omega * t);
Ts_fine = 0.001;
t_fine = 0:Ts_fine:Tfinal;
u_fine = u_func(t_fine);

q = zeros(size(t_fine));
dq = zeros(size(t_fine));
for k = 1:length(t_fine)-1
    ddq = (1 / (m * L^2)) * (u_fine(k) - c * dq(k) - m * g * L * q(k));
    dq(k+1) = dq(k) + Ts_fine * ddq;
    q(k+1) = q(k) + Ts_fine * dq(k);
end

% Μερική μέτρηση q(t), u(t)
Ts_b = 0.2;
t_b = 0:Ts_b:Tfinal;
q_b = interp1(t_fine, q, t_b);
u_b = interp1(t_fine, u_fine, t_b);

dq_b = diff(q_b)/Ts_b;
q_b = q_b(1:end-1);
u_b = u_b(1:end-1);

Phi_b = [q_b', u_b'];
theta_b = (Phi_b' * Phi_b) \ (Phi_b' * dq_b');
a_hat = theta_b(1);
b_hat = theta_b(2);

q_hat_b = zeros(size(q_b));
for k = 1:length(q_hat_b)-1
    dq_est = a_hat * q_hat_b(k) + b_hat * u_b(k);
    q_hat_b(k+1) = q_hat_b(k) + Ts_b * dq_est;
end
eq_b = q_b - q_hat_b;

% Διαγράμματα
figure('Name','Μερική μέτρηση q(t), u(t)');
subplot(3,1,1); plot(t_b(1:end-1), q_b, 'b', 'LineWidth', 1.5);
ylabel('q(t)'); title('Πραγματική γωνία q(t)'); grid on;
subplot(3,1,2); plot(t_b(1:end-1), q_hat_b, 'r', 'LineWidth', 1.5);
ylabel('q̂(t)'); title('Εκτιμημένη γωνία q̂(t)'); grid on;
subplot(3,1,3); plot(t_b(1:end-1), eq_b, 'k', 'LineWidth', 1.5);
ylabel('e_q(t)'); xlabel('t [sec]');
title('Σφάλμα e_q(t) = q(t) - q̂(t)'); grid on;
