clear;

syms alpha beta gamma phi theta psi real
syms dalpha dbeta dgamma dphi dtheta dpsi real
syms ddalpha ddbeta ddgamma ddphi ddtheta ddpsi real
syms m1 m2 m3 real
syms I1x I1y I1z real
syms I2x I2y I2z real
syms I3x I3y I3z real

syms k1 k2 k3 real % "spring" constants
syms u1 u2 u3 real % inputs

q = [psi; theta; phi];
dq = [dpsi; dtheta; dphi];
ddq = [ddpsi; ddtheta; ddphi];

R1 = rotz(psi);
R2 = roty(theta);
R3 = rotx(phi);

R = R1 * R2 * R3;
R00 = eye(3);
R01 = R1;
R02 = R1*R2;
R03 = R;

I1 = diag([I1x, I1y, I1z]);
I2 = diag([I2x, I2y, I2z]);
I3 = diag([I3x, I3y, I3z]);

J1 = R00(:,3);
J2 = R01(:,2);
J3 = R02(:,1);

J1_ = [J1, zeros(3, 2)];
J2_ = [J1, J2, zeros(3,1)];
J3_ = [J1, J2, J3];

D1 = J1_' * R01 * I1 * R01' * J1_;
D2 = J2_' * R02 * I2 * R02' * J2_;
D3 = J3_' * R03 * I3 * R03' * J3_;
D = simplify(D1 + D2 + D3, 100);

C = christoffel(D, q);

% Assume center of rotation is center of mass.


tau1 = k1 * sin(u1 - q(1));
tau2 = k2 * sin(u2 - q(2));
tau3 = k3 * sin(u3 - q(3));

tau = [tau1; tau2; tau3];

deq = simplify(D*ddq + C * dq, 100) == tau;

x = [q; dq];

f = simplify([dq; D^-1*(-C*dq+tau)], 100);

f_q = simplify(D^-1 * (-C * dq + tau), 100);

matlabFunction(f, 'File', 'get_f');
matlabFunction(f_q, 'File', 'get_f_q');
save('matrices', 'D', 'C', 'tau')


