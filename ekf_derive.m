pkg load symbolic

% EKF state variables
syms q0 q1 q2 q3 bp bq br
X = [q0; q1; q2; q3; bp; bq; br];

% f(x)
syms f(X) p q r
f_q = [1/2 * [-q1 -q2 -q3; q0 -q3 q2; q3 q0 -q1; -q2 q1 q0]] * [p-bp; q-bq; r-br];
f(X) = [f_q; 0; 0; 0];

% h(X)
syms h(X) g
h11 = 2*g*(q1*q3 - q0*q2);
h21 = 2*g*(q2*q3 + q0*q1);
h31 = g*(q0^2-q1^2-q2^2+q3^2);
h41 = atan((2*(q1*q2+q0*q3))/(q0^2-q1^2-q2^2+q3^2));
h(X) = [h11; h21; h31; h41];

% F
F = jacobian(f(X));

% H
H = jacobian(h(X));
