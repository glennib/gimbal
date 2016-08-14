clear;

syms alpha beta gamma phi theta psi real
syms dalpha dbeta dgamma dphi dtheta dpsi real
syms ddalpha ddbeta ddgamma ddphi ddtheta ddpsi real
syms tau1 tau2 tau3 real
syms m1 m2 m3 real
syms I1x I1y I1z real
syms I2x I2y I2z real
syms I3x I3y I3z real

q = [psi; theta; phi];
dq = [dpsi; dtheta; dphi];

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
J2 = R01(:,3);
J3 = R02(:,3);

J1_ = [J1, zeros(3, 2)];
J2_ = [J1, J2, zeros(3,1)];
J3_ = [J1, J2, J3];

K1 = 0.5 * dq' * J1_' * R1 * I1 * R1' * J1_ * dq;
K2 = 0.5 * dq' * J2_' * R2 * I2 * R2' * J2_ * dq;
K3 = 0.5 * dq' * J3_' * R3 * I3 * R3' * J3_ * dq;

K = K1 + K2 + K3;
K = simplify(K);

% Assume center of rotation is center of mass.

P = 0;

dKddq = [
    diff(K, dq(1));
    diff(K, dq(2));
    diff(K, dq(3))
    ];

dKdq = [
    diff(K, q(1));
    diff(K, q(2));
    diff(K, q(3));
    ];

cg = ([I3z * ddpsi * cos(phi) ^ 2 + I3z * ddtheta * cos(phi) ^ 2 - I3y * cos(phi) ^ 2 * ddpsi + I2z * ddpsi * cos(theta) ^ 2 - I2x * cos(theta) ^ 2 * ddpsi + I3y * ddphi * cos(theta) - I3y * cos(phi) ^ 2 * ddtheta + I2z * ddtheta * cos(theta) ^ 2 - I2x * cos(theta) ^ 2 * ddtheta + I2x * ddpsi + I3y * ddpsi + I1z * ddpsi + I2x * ddtheta + I3y * ddtheta + 0.2e1 * I3y * dphi ^ 2 * cos(phi) ^ 2 * sin(psi) * sin(theta) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) ^ 2 * sin(psi) * sin(theta) + 0.2e1 * I3y * cos(phi) * cos(theta) * dphi ^ 2 * sin(phi) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) * cos(theta) * sin(phi) + 0.2e1 * I3y * cos(phi) * dtheta * dphi * sin(phi) - 0.2e1 * I3z * dtheta * cos(phi) * dphi * sin(phi) + 0.2e1 * I3y * cos(phi) * dpsi * dphi * sin(phi) - 0.2e1 * I3z * dpsi * cos(phi) * dphi * sin(phi) + I3y * cos(phi) ^ 2 * sin(theta) * dphi * dtheta - I3y * cos(phi) ^ 2 * cos(theta) * ddphi + I3z * ddphi * cos(phi) ^ 2 * cos(theta) + cos(psi) * sin(theta) * (-I3y * cos(phi) ^ 2 * sin(psi) * sin(theta) * dphi + I3z * cos(phi) ^ 2 * sin(psi) * sin(theta) * dphi - I3y * cos(phi) * cos(theta) * sin(phi) * dphi + I3z * cos(phi) * cos(theta) * sin(phi) * dphi + I3x * sin(psi) * sin(theta) * dphi - I3y * cos(phi) * sin(phi) * dtheta - I3y * cos(phi) * sin(phi) * dpsi + I3z * cos(phi) * sin(phi) * dtheta + I3z * cos(phi) * sin(phi) * dpsi - I3z * sin(psi) * sin(theta) * dphi) * dphi - I3y * dphi ^ 2 * sin(psi) * sin(theta) + I3z * dphi ^ 2 * sin(psi) * sin(theta) - I3y * dphi * dtheta * sin(theta) - I3z * ddphi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * ddphi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * dphi * dtheta * cos(phi) * cos(theta) * sin(phi) * sin(psi) - I3z * dphi * dtheta * cos(phi) * cos(theta) * sin(phi) * sin(psi) + I3y * dphi * dpsi * cos(phi) * cos(psi) * sin(phi) * sin(theta) - I3z * dphi * dpsi * cos(phi) * cos(psi) * sin(phi) * sin(theta) - I3z * dphi * dtheta * cos(phi) ^ 2 * sin(theta) + 0.2e1 * I2x * dpsi * dtheta * cos(theta) * sin(theta) - 0.2e1 * I2z * dpsi * dtheta * cos(theta) * sin(theta) + 0.2e1 * I2x * dtheta ^ 2 * cos(theta) * sin(theta) - 0.2e1 * I2z * dtheta ^ 2 * cos(theta) * sin(theta) I3z * ddpsi * cos(phi) ^ 2 + I3z * ddtheta * cos(phi) ^ 2 - I3y * cos(phi) ^ 2 * ddpsi + I2z * ddpsi * cos(theta) ^ 2 - I2x * cos(theta) ^ 2 * ddpsi + I3y * ddphi * cos(theta) - I3y * cos(phi) ^ 2 * ddtheta + I2z * ddtheta * cos(theta) ^ 2 - I2x * cos(theta) ^ 2 * ddtheta + I2x * ddpsi + I3y * ddpsi + I2x * ddtheta + I3y * ddtheta + 0.2e1 * I3y * dphi ^ 2 * cos(phi) ^ 2 * sin(psi) * sin(theta) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) ^ 2 * sin(psi) * sin(theta) + 0.2e1 * I3y * cos(phi) * cos(theta) * dphi ^ 2 * sin(phi) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) * cos(theta) * sin(phi) + 0.2e1 * I3y * cos(phi) * dtheta * dphi * sin(phi) - 0.2e1 * I3z * dtheta * cos(phi) * dphi * sin(phi) + 0.2e1 * I3y * cos(phi) * dpsi * dphi * sin(phi) - 0.2e1 * I3z * dpsi * cos(phi) * dphi * sin(phi) - I3y * cos(phi) ^ 2 * sin(theta) * dphi * dpsi - 0.2e1 * I3y * dphi ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(theta) + I3z * cos(psi) ^ 2 * cos(theta) * sin(theta) * dphi ^ 2 + I3y * dphi ^ 2 * cos(phi) * sin(phi) * sin(psi) - I3z * dphi ^ 2 * cos(phi) * sin(phi) * sin(psi) - I3y * cos(phi) ^ 2 * cos(theta) * ddphi + I3z * ddphi * cos(phi) ^ 2 * cos(theta) + I3y * dphi ^ 2 * cos(theta) * sin(theta) - I3z * dphi ^ 2 * cos(theta) * sin(theta) - I3y * dphi ^ 2 * sin(psi) * sin(theta) + I3z * dphi ^ 2 * sin(psi) * sin(theta) + I3y * dphi * dpsi * sin(theta) - I3z * ddphi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * ddphi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * cos(psi) ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(theta) * dphi ^ 2 - I3z * cos(psi) ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(theta) * dphi ^ 2 - I3y * dphi * dpsi * cos(phi) * cos(theta) * sin(phi) * sin(psi) + I3z * dphi * dpsi * cos(phi) * cos(theta) * sin(phi) * sin(psi) + I3y * dphi * dpsi * cos(phi) * cos(psi) * sin(phi) * sin(theta) - I3z * dphi * dpsi * cos(phi) * cos(psi) * sin(phi) * sin(theta) + I3z * dphi * dpsi * cos(phi) ^ 2 * sin(theta) + 0.2e1 * I3z * dphi ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(theta) - I3x * dphi ^ 2 * cos(psi) ^ 2 * cos(theta) * sin(theta) - I2x * dpsi ^ 2 * cos(theta) * sin(theta) + I2z * dpsi ^ 2 * cos(theta) * sin(theta) + I2x * dtheta ^ 2 * cos(theta) * sin(theta) - I2z * dtheta ^ 2 * cos(theta) * sin(theta) - 0.2e1 * I3y * dphi ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) * sin(psi) + 0.2e1 * I3z * dphi ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) * sin(psi) I3x * ddphi * cos(psi) ^ 2 + I3y * ddpsi * cos(theta) + I3y * ddphi * cos(phi) ^ 2 - I3y * dtheta ^ 2 * sin(theta) - I3z * cos(phi) ^ 2 * ddphi - I3z * cos(psi) ^ 2 * ddphi - I3z * cos(theta) ^ 2 * ddphi + I3y * ddphi * cos(theta) ^ 2 + I3y * ddtheta * cos(theta) + I3z * ddphi - 0.2e1 * I3z * dphi * cos(phi) * cos(theta) * sin(phi) * dpsi * cos(psi) * sin(theta) + 0.2e1 * I3y * dphi * cos(phi) * cos(theta) * sin(phi) * dpsi * cos(psi) * sin(theta) + I3y * cos(phi) ^ 2 * cos(theta) ^ 2 * cos(psi) ^ 2 * ddphi - I3z * cos(phi) ^ 2 * cos(theta) ^ 2 * cos(psi) ^ 2 * ddphi - I3z * dpsi * cos(phi) ^ 2 * dtheta * sin(theta) + I3y * cos(phi) ^ 2 * dtheta * sin(theta) * dpsi - 0.2e1 * I3x * dphi * cos(psi) * dpsi * sin(psi) + 0.2e1 * I3z * cos(psi) * dphi * dpsi * sin(psi) + 0.2e1 * I3z * cos(theta) * dphi * dtheta * sin(theta) - 0.2e1 * I3y * dphi * cos(theta) * dtheta * sin(theta) - I3y * dphi ^ 2 * cos(theta) * sin(psi) * sin(theta) + I3z * dphi ^ 2 * cos(theta) * sin(psi) * sin(theta) + I3y * cos(psi) ^ 2 * cos(phi) * sin(phi) * dphi ^ 2 - I3z * cos(psi) ^ 2 * cos(phi) * sin(phi) * dphi ^ 2 + 0.2e1 * I3y * ddphi * cos(phi) * cos(theta) * sin(phi) * sin(psi) * sin(theta) - 0.2e1 * I3z * ddphi * cos(phi) * cos(theta) * sin(phi) * sin(psi) * sin(theta) + I3y * dpsi * cos(phi) * sin(phi) * sin(psi) * dtheta * cos(theta) - I3z * dpsi * cos(phi) * sin(phi) * sin(psi) * dtheta * cos(theta) + I3y * dtheta * cos(phi) * sin(phi) * dpsi * cos(psi) * sin(theta) - I3z * dtheta * cos(phi) * sin(phi) * dpsi * cos(psi) * sin(theta) + 0.4e1 * I3y * dphi * cos(phi) * cos(theta) ^ 2 * sin(phi) * sin(psi) * dtheta - 0.4e1 * I3z * dphi * cos(phi) * cos(theta) ^ 2 * sin(phi) * sin(psi) * dtheta - 0.2e1 * I3y * cos(phi) ^ 2 * cos(theta) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) + 0.2e1 * I3z * cos(phi) ^ 2 * cos(theta) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) - 0.2e1 * I3y * cos(phi) ^ 2 * cos(theta) * cos(psi) ^ 2 * dphi * dtheta * sin(theta) + 0.2e1 * I3z * cos(phi) ^ 2 * cos(theta) * cos(psi) ^ 2 * dphi * dtheta * sin(theta) - I3x * cos(theta) ^ 2 * cos(psi) ^ 2 * ddphi + I3z * cos(theta) ^ 2 * cos(psi) ^ 2 * ddphi - I3y * dpsi * dtheta * sin(theta) - I3y * cos(phi) ^ 2 * cos(psi) ^ 2 * ddphi + I3z * cos(phi) ^ 2 * cos(psi) ^ 2 * ddphi - I3y * cos(phi) ^ 2 * cos(theta) * ddpsi + I3z * ddpsi * cos(phi) ^ 2 * cos(theta) + I3z * ddtheta * cos(phi) ^ 2 * cos(theta) - I3y * cos(phi) ^ 2 * cos(theta) * ddtheta + 0.2e1 * I3z * ddphi * cos(phi) ^ 2 * cos(theta) ^ 2 - 0.2e1 * I3y * cos(phi) ^ 2 * cos(theta) ^ 2 * ddphi - I3z * dtheta ^ 2 * cos(phi) ^ 2 * sin(theta) + I3y * cos(phi) ^ 2 * dtheta ^ 2 * sin(theta) - I3y * dphi ^ 2 * cos(phi) * sin(phi) + I3z * dphi ^ 2 * cos(phi) * sin(phi) - 0.4e1 * I3z * dphi * cos(phi) ^ 2 * cos(theta) * dtheta * sin(theta) + 0.4e1 * I3y * cos(phi) ^ 2 * cos(theta) * dphi * dtheta * sin(theta) + 0.2e1 * I3y * cos(phi) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) - 0.2e1 * I3z * cos(phi) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) - 0.2e1 * I3y * dphi * cos(phi) * dtheta * sin(phi) * sin(psi) + 0.2e1 * I3z * dphi * cos(phi) * dtheta * sin(phi) * sin(psi) + I3y * ddpsi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * ddtheta * cos(phi) * sin(phi) * sin(psi) * sin(theta) - I3z * ddtheta * cos(phi) * sin(phi) * sin(psi) * sin(theta) + I3y * dpsi ^ 2 * cos(phi) * sin(phi) * cos(psi) * sin(theta) - I3z * dpsi ^ 2 * cos(phi) * sin(phi) * cos(psi) * sin(theta) + I3y * dtheta ^ 2 * cos(phi) * sin(phi) * sin(psi) * cos(theta) - I3z * dtheta ^ 2 * cos(phi) * sin(phi) * sin(psi) * cos(theta) - I3z * ddpsi * cos(phi) * sin(phi) * sin(psi) * sin(theta) + 0.2e1 * I3x * cos(theta) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) - 0.2e1 * I3z * cos(theta) ^ 2 * cos(psi) * dphi * dpsi * sin(psi) + 0.2e1 * I3x * cos(theta) * cos(psi) ^ 2 * dphi * dtheta * sin(theta) - 0.2e1 * I3z * cos(theta) * cos(psi) ^ 2 * dphi * dtheta * sin(theta) - I3y * cos(psi) ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) * dphi ^ 2 + I3z * cos(psi) ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) * dphi ^ 2 - 0.2e1 * I3y * dpsi * dtheta * cos(phi) * sin(phi) + 0.2e1 * I3z * dpsi * dtheta * cos(phi) * sin(phi) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) + 0.2e1 * I3y * dphi ^ 2 * cos(phi) * cos(theta) ^ 2 * sin(phi) - I3y * dpsi ^ 2 * cos(phi) * sin(phi) + I3z * dpsi ^ 2 * cos(phi) * sin(phi) - I3y * dtheta ^ 2 * cos(phi) * sin(phi) + I3z * dtheta ^ 2 * cos(phi) * sin(phi) + 0.2e1 * I3y * dphi ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(psi) * sin(theta) - 0.2e1 * I3z * dphi ^ 2 * cos(phi) ^ 2 * cos(theta) * sin(psi) * sin(theta)]) == ([tau1 tau2 tau3])
