clc; clear;

% ΠΡΑΓΜΑΤΙΚΕΣ ΠΑΡΑΜΕΤΡΟΙ
m = 0.75; L = 1.25; c = 0.15; g = 9.81;
A0 = 4; omega = 2;
Tfinal = 20;
Ts_fine = 0.001;
t_fine = 0:Ts_fine:Tfinal;
u_fine = A0 * sin(omega * t_fine);

% ΠΡΑΓΜΑΤΙΚΗ ΠΡΟΣΟΜΟΙΩΣΗ
q = zeros(size(t_fine));
dq = zeros(size(t_fine));
for k = 1:length(t_fine)-1
    ddq = (1 / (m * L^2)) * (u_fine(k) - c * dq(k) - m * g * L * q(k));
    dq(k+1) = dq(k) + Ts_fine * ddq;
    q(k+1) = q(k) + Ts_fine * dq(k);
end

% ΔΕΙΓΜΑΤΟΛΗΨΙΑ
Ts = 0.1;
t = 0:Ts:Tfinal;
q_s = interp1(t_fine, q, t);
dq_s = interp1(t_fine, dq, t);
u_s = A0 * sin(omega * t);

% ΥΠΟΛΟΓΙΣΜΟΣ ΠΑΡΑΓΩΓΩΝ
ddq_s = diff(dq_s) / Ts;
q_s = q_s(1:end-1);
dq_s = dq_s(1:end-1);
u_s = u_s(1:end-1);

% ΠΡΟΣΘΗΚΗ ΛΕΥΚΟΥ ΓΚΑΟΥΣΙΑΝΟΥ ΘΟΡΥΒΟΥ
noise_level = 0.02; % Θόρυβος ~ 2%
q_noisy = q_s + noise_level * randn(size(q_s));
dq_noisy = dq_s + noise_level * randn(size(dq_s));
ddq_noisy = ddq_s + noise_level * randn(size(ddq_s));

% ΠΙΝΑΚΑΣ ΠΑΡΑΤΗΡΗΣΙΜΟΤΗΤΑΣ ΚΑΙ LEAST SQUARES
Phi_noisy = [ddq_noisy', dq_noisy', q_noisy'];
theta_hat_noisy = (Phi_noisy' * Phi_noisy) \ (Phi_noisy' * u_s');

% ΠΡΑΓΜΑΤΙΚΕΣ ΠΑΡΑΜΕΤΡΟΙ ΓΙΑ ΣΥΓΚΡΙΣΗ
theta_true = [m*L^2; c; m*g*L];

% ΕΜΦΑΝΙΣΗ ΑΠΟΤΕΛΕΣΜΑΤΩΝ
disp('--- Εκτιμώμενες Παράμετροι με Θόρυβο ---');
disp(theta_hat_noisy);

disp('--- Πραγματικές Παράμετροι ---');
disp(theta_true);

disp('--- Απόλυτο Σφάλμα ---');
disp(abs(theta_hat_noisy - theta_true));
