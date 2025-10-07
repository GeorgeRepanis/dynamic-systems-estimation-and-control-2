clc; clear;

% Πραγματικές παράμετροι
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
A0 = 4; omega = 2;
Tfinal = 20;
Ts_fine = 0.001;
t_fine = 0:Ts_fine:Tfinal;
u_fine = A0 * sin(omega * t_fine);

% Προσομοίωση "πραγματικού" συστήματος
q = zeros(size(t_fine));
dq = zeros(size(t_fine));
for k = 1:length(t_fine)-1
    ddq = (1 / (m * L^2)) * (u_fine(k) - c * dq(k) - m * g * L * q(k));
    dq(k+1) = dq(k) + Ts_fine * ddq;
    q(k+1) = q(k) + Ts_fine * dq(k);
end

% Δειγματοληψία
Ts = 0.1;
t = 0:Ts:Tfinal;
q_s = interp1(t_fine, q, t);
dq_s = interp1(t_fine, dq, t);
u_s = A0 * sin(omega * t);

ddq_s = diff(dq_s) / Ts;
q_trim = q_s(1:end-1); dq_trim = dq_s(1:end-1); u_trim = u_s(1:end-1);

% Εκτίμηση χωρίς θόρυβο
Phi_clean = [ddq_s', dq_trim', q_trim'];
theta_clean = (Phi_clean' * Phi_clean) \ (Phi_clean' * u_trim');

% Εκτίμηση με θόρυβο
noise_level = 0.02;
q_noisy = q_trim + noise_level * randn(size(q_trim));
dq_noisy = dq_trim + noise_level * randn(size(dq_trim));
ddq_noisy = ddq_s + noise_level * randn(size(ddq_s));
Phi_noisy = [ddq_noisy', dq_noisy', q_noisy'];
theta_noisy = (Phi_noisy' * Phi_noisy) \ (Phi_noisy' * u_trim');

% Ανακατασκευή q̂(t) χωρίς θόρυβο
q_hat_clean = zeros(size(q_trim));
dq_hat_clean = zeros(size(dq_trim));
for k = 1:length(q_hat_clean)-1
    ddq_est = (1 / theta_clean(1)) * (u_trim(k) - theta_clean(2)*dq_hat_clean(k) - theta_clean(3)*q_hat_clean(k));
    dq_hat_clean(k+1) = dq_hat_clean(k) + Ts * ddq_est;
    q_hat_clean(k+1) = q_hat_clean(k) + Ts * dq_hat_clean(k);
end

% Ανακατασκευή q̂(t) με θόρυβο
q_hat_noisy = zeros(size(q_trim));
dq_hat_noisy = zeros(size(dq_trim));
for k = 1:length(q_hat_noisy)-1
    ddq_est = (1 / theta_noisy(1)) * (u_trim(k) - theta_noisy(2)*dq_hat_noisy(k) - theta_noisy(3)*q_hat_noisy(k));
    dq_hat_noisy(k+1) = dq_hat_noisy(k) + Ts * ddq_est;
    q_hat_noisy(k+1) = q_hat_noisy(k) + Ts * dq_hat_noisy(k);
end

% Σφάλματα
eq_clean = q_trim - q_hat_clean;
eq_noisy = q_trim - q_hat_noisy;

% ΓΡΑΦΙΚΗ ΣΥΓΚΡΙΣΗ
figure('Name','Θέμα 3(α) - Σύγκριση Εκτίμησης με και χωρίς θόρυβο');

subplot(3,1,1);
plot(t(1:end-1), q_trim, 'b', 'LineWidth', 1.5); hold on;
plot(t(1:end-1), q_hat_clean, 'g--', 'LineWidth', 1.5);
plot(t(1:end-1), q_hat_noisy, 'r:', 'LineWidth', 1.5);
legend('q(t)', 'q̂(t) καθαρό', 'q̂(t) με θόρυβο');
title('Σύγκριση πραγματικής και εκτιμημένης γωνίας');
ylabel('q(t)');
grid on;

subplot(3,1,2);
plot(t(1:end-1), eq_clean, 'k', 'LineWidth', 1.5);
title('Σφάλμα χωρίς θόρυβο: e_q(t) = q(t) - q̂(t)');
ylabel('e_q(t)');
grid on;

subplot(3,1,3);
plot(t(1:end-1), eq_noisy, 'm', 'LineWidth', 1.5);
title('Σφάλμα με θόρυβο: e_q(t) = q(t) - q̂_noisy(t)');
ylabel('e_q(t)');
xlabel('t [sec]');
grid on;

