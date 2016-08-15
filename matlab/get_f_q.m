function f_q = get_f_q(I2x,I3x,I2y,I3y,I1z,I2z,I3z,dpsi,dtheta,k1,k2,k3,phi,psi,theta,u1,u2,u3)
%GET_F_Q
%    F_Q = GET_F_Q(I2X,I3X,I2Y,I3Y,I1Z,I2Z,I3Z,DPSI,DTHETA,K1,K2,K3,PHI,PSI,THETA,U1,U2,U3)

%    This function was generated by the Symbolic Math Toolbox version 7.0.
%    15-Aug-2016 13:52:00

t2 = cos(phi);
t3 = t2.^2;
t4 = cos(theta);
t5 = t4.^2;
t6 = psi-u1;
t7 = sin(t6);
t8 = I3z.^2;
t9 = sin(theta);
t10 = sin(u3);
t11 = I3y.^2;
t12 = sin(phi);
t13 = cos(u3);
t14 = sin(u2);
t15 = cos(u2);
t16 = I2x.*I2y;
t17 = I2x.*I3z;
t18 = I2y.*I1z;
t19 = I1z.*I3z;
t20 = I2x.*I3y.*t3;
t21 = I3y.*I1z.*t3;
t22 = I2y.*I3y.*t5;
t23 = I2y.*I2z.*t5;
t24 = I3y.*I3z.*t5;
t25 = I2z.*I3z.*t5;
t26 = I2x.*I3z.*t3.*t5;
t27 = I2y.*I3z.*t3.*t5;
t28 = I3y.*I2z.*t3.*t5;
t33 = I2x.*I3z.*t3;
t34 = I1z.*I3z.*t3;
t35 = I2x.*I2y.*t5;
t36 = I2x.*I3z.*t5;
t37 = I2x.*I3y.*t3.*t5;
t38 = I2y.*I3y.*t3.*t5;
t39 = I2z.*I3z.*t3.*t5;
t29 = t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28-t33-t34-t35-t36-t37-t38-t39;
t30 = 1.0./t29;
t31 = phi.*2.0;
t32 = sin(t31);
t40 = theta.*2.0;
t41 = sin(t40);
t42 = I3x.^2;
f_q = [t30.*(I2y.*k1.*t7.*2.0+I3z.*k1.*t7.*2.0+dtheta.*t4.*t8-dtheta.*t3.*t4.*t8+dtheta.*t3.*t4.*t11-I3x.*I2y.*dtheta.*t4-I2y.*I3y.*dtheta.*t4-I3x.*I3z.*dtheta.*t4+I2y.*I3z.*dtheta.*t4-I3y.*I3z.*dtheta.*t4+I3y.*k1.*t3.*t7.*2.0-I3z.*k1.*t3.*t7.*2.0-I3x.*I3y.*dtheta.*t3.*t4+I2y.*I3y.*dtheta.*t3.*t4.*2.0+I3x.*I3z.*dtheta.*t3.*t4-I2y.*I3z.*dtheta.*t3.*t4.*2.0-I2y.*k3.*t2.*t9.*t10.*2.0+I2y.*k3.*t9.*t12.*t13.*2.0-I3z.*k3.*t2.*t9.*t10.*2.0+I3z.*k3.*t9.*t12.*t13.*2.0-dpsi.*t2.*t5.*t8.*t12+dpsi.*t2.*t5.*t11.*t12-I3x.*I3y.*dpsi.*t2.*t5.*t12+I2y.*I3y.*dpsi.*t2.*t5.*t12.*2.0+I3x.*I3z.*dpsi.*t2.*t5.*t12-I2y.*I3z.*dpsi.*t2.*t5.*t12.*2.0-I3y.*k3.*t2.*t3.*t9.*t10.*2.0+I3y.*k2.*t2.*t5.*t12.*t14.*2.0+I3y.*k3.*t3.*t9.*t12.*t13.*2.0+I3z.*k3.*t2.*t3.*t9.*t10.*2.0-I3z.*k2.*t2.*t5.*t12.*t14.*2.0-I3z.*k3.*t3.*t9.*t12.*t13.*2.0-I3y.*k2.*t2.*t4.*t9.*t12.*t15.*2.0+I3z.*k2.*t2.*t4.*t9.*t12.*t15.*2.0).*(-1.0./2.0);t30.*(-dpsi.*t4.*t5.*t11+I2x.*I3x.*dpsi.*t4-I2x.*I3y.*dpsi.*t4+I3x.*I1z.*dpsi.*t4+I2x.*I3z.*dpsi.*t4-I3y.*I1z.*dpsi.*t4+I1z.*I3z.*dpsi.*t4-I2x.*I3y.*dtheta.*t32+I2x.*I3z.*dtheta.*t32-I3y.*I1z.*dtheta.*t32+I1z.*I3z.*dtheta.*t32-I2x.*k2.*t4.*t14.*2.0+I2x.*k2.*t9.*t15.*2.0-I1z.*k2.*t4.*t14.*2.0+I1z.*k2.*t9.*t15.*2.0-I2x.*I3x.*dpsi.*t4.*t5+I2x.*I3y.*dpsi.*t3.*t4.*2.0+I2x.*I3y.*dpsi.*t4.*t5+I3x.*I3y.*dpsi.*t4.*t5-I2x.*I3z.*dpsi.*t3.*t4.*2.0-I2x.*I3z.*dpsi.*t4.*t5+I3x.*I2z.*dpsi.*t4.*t5+I3y.*I1z.*dpsi.*t3.*t4.*2.0-I3y.*I2z.*dpsi.*t4.*t5+I3y.*I3z.*dpsi.*t4.*t5-I1z.*I3z.*dpsi.*t3.*t4.*2.0+I2z.*I3z.*dpsi.*t4.*t5+I2x.*k2.*t4.*t5.*t14.*2.0-I2x.*k2.*t5.*t9.*t15.*2.0-I3y.*k2.*t4.*t5.*t14.*2.0+I3y.*k2.*t5.*t9.*t15.*2.0-I2z.*k2.*t4.*t5.*t14.*2.0+I2z.*k2.*t5.*t9.*t15.*2.0-dpsi.*t3.*t4.*t5.*t8+dpsi.*t3.*t4.*t5.*t11+dtheta.*t2.*t5.*t8.*t12-dtheta.*t2.*t5.*t11.*t12-I2x.*I3y.*dpsi.*t3.*t4.*t5.*2.0-I3x.*I3y.*dpsi.*t3.*t4.*t5+I2x.*I3z.*dpsi.*t3.*t4.*t5.*2.0+I3x.*I3z.*dpsi.*t3.*t4.*t5+I3y.*I2z.*dpsi.*t3.*t4.*t5.*2.0-I2z.*I3z.*dpsi.*t3.*t4.*t5.*2.0+I2x.*I3y.*dtheta.*t2.*t5.*t12.*2.0+I3x.*I3y.*dtheta.*t2.*t5.*t12-I2x.*I3z.*dtheta.*t2.*t5.*t12.*2.0-I3x.*I3z.*dtheta.*t2.*t5.*t12-I3y.*I2z.*dtheta.*t2.*t5.*t12.*2.0+I2z.*I3z.*dtheta.*t2.*t5.*t12.*2.0-I3y.*k1.*t2.*t4.*t7.*t12.*2.0+I3y.*k2.*t3.*t4.*t5.*t14.*2.0-I3y.*k3.*t2.*t4.*t9.*t13.*2.0-I3y.*k2.*t3.*t5.*t9.*t15.*2.0+I3z.*k1.*t2.*t4.*t7.*t12.*2.0-I3z.*k2.*t3.*t4.*t5.*t14.*2.0+I3z.*k3.*t2.*t4.*t9.*t13.*2.0+I3z.*k2.*t3.*t5.*t9.*t15.*2.0+I3y.*k3.*t2.*t3.*t4.*t9.*t13.*2.0+I3y.*k3.*t3.*t4.*t9.*t10.*t12.*2.0-I3z.*k3.*t2.*t3.*t4.*t9.*t13.*2.0-I3z.*k3.*t3.*t4.*t9.*t10.*t12.*2.0).*(-1.0./2.0);(t30.*(I3x.*dtheta.*t8.*t41.*(1.0./2.0)-I2y.*dtheta.*t41.*t42.*(1.0./2.0)-I3z.*dtheta.*t41.*t42.*(1.0./2.0)-I3x.*I2y.*I3y.*dtheta.*t41.*(1.0./2.0)+I3x.*I2y.*I3z.*dtheta.*t41.*(1.0./2.0)-I3x.*I3y.*I3z.*dtheta.*t41.*(1.0./2.0)-I2x.*I2y.*k3.*t2.*t10.*2.0-I3x.*I2y.*k3.*t2.*t10.*2.0+I3x.*I2y.*k1.*t7.*t9.*2.0+I2x.*I2y.*k3.*t12.*t13.*2.0+I3x.*I2y.*k3.*t12.*t13.*2.0-I2x.*I3z.*k3.*t2.*t10.*2.0-I3x.*I3z.*k3.*t2.*t10.*2.0+I3x.*I3z.*k1.*t7.*t9.*2.0+I2x.*I3z.*k3.*t12.*t13.*2.0+I3x.*I3z.*k3.*t12.*t13.*2.0-I2y.*I1z.*k3.*t2.*t10.*2.0+I2y.*I1z.*k3.*t12.*t13.*2.0-I1z.*I3z.*k3.*t2.*t10.*2.0+I1z.*I3z.*k3.*t12.*t13.*2.0-I2x.*I3y.*k3.*t2.*t3.*t10.*2.0+I2x.*I2y.*k3.*t2.*t5.*t10.*2.0-I3x.*I3y.*k3.*t2.*t3.*t10.*2.0+I3x.*I2y.*k3.*t2.*t5.*t10.*2.0+I3x.*I3y.*k1.*t3.*t7.*t9.*2.0+I2x.*I3y.*k3.*t3.*t12.*t13.*2.0-I2x.*I2y.*k3.*t5.*t12.*t13.*2.0+I3x.*I3y.*k3.*t3.*t12.*t13.*2.0-I3x.*I2y.*k3.*t5.*t12.*t13.*2.0+I2x.*I3z.*k3.*t2.*t3.*t10.*2.0+I3x.*I3z.*k3.*t2.*t3.*t10.*2.0+I2x.*I3z.*k3.*t2.*t5.*t10.*2.0-I2y.*I3y.*k3.*t2.*t5.*t10.*2.0-I3x.*I3z.*k1.*t3.*t7.*t9.*2.0+I3x.*I3z.*k3.*t2.*t5.*t10.*2.0-I2x.*I3z.*k3.*t3.*t12.*t13.*2.0-I3x.*I3z.*k3.*t3.*t12.*t13.*2.0-I2x.*I3z.*k3.*t5.*t12.*t13.*2.0+I2y.*I3y.*k3.*t5.*t12.*t13.*2.0-I3x.*I3z.*k3.*t5.*t12.*t13.*2.0-I3y.*I1z.*k3.*t2.*t3.*t10.*2.0-I2y.*I2z.*k3.*t2.*t5.*t10.*2.0-I3y.*I3z.*k3.*t2.*t5.*t10.*2.0+I3y.*I1z.*k3.*t3.*t12.*t13.*2.0+I2y.*I2z.*k3.*t5.*t12.*t13.*2.0+I3y.*I3z.*k3.*t5.*t12.*t13.*2.0+I1z.*I3z.*k3.*t2.*t3.*t10.*2.0-I2z.*I3z.*k3.*t2.*t5.*t10.*2.0-I1z.*I3z.*k3.*t3.*t12.*t13.*2.0+I2z.*I3z.*k3.*t5.*t12.*t13.*2.0-I3x.*dtheta.*t3.*t4.*t8.*t9+I3x.*dtheta.*t3.*t4.*t9.*t11-I3y.*dtheta.*t3.*t4.*t9.*t42+I3z.*dtheta.*t3.*t4.*t9.*t42+I2x.*I3y.*k3.*t2.*t3.*t5.*t10.*2.0+I3x.*I3y.*k3.*t2.*t3.*t5.*t10.*2.0-I2x.*I3y.*k3.*t3.*t5.*t12.*t13.*2.0-I3x.*I3y.*k2.*t2.*t4.*t12.*t15.*2.0-I3x.*I3y.*k3.*t3.*t5.*t12.*t13.*2.0-I2x.*I3z.*k3.*t2.*t3.*t5.*t10.*2.0+I2y.*I3y.*k3.*t2.*t3.*t5.*t10.*2.0-I3x.*I3z.*k3.*t2.*t3.*t5.*t10.*2.0+I2x.*I3z.*k3.*t3.*t5.*t12.*t13.*2.0+I3x.*I3z.*k2.*t2.*t4.*t12.*t15.*2.0-I2y.*I3y.*k3.*t3.*t5.*t12.*t13.*2.0+I3x.*I3z.*k3.*t3.*t5.*t12.*t13.*2.0-I2y.*I3z.*k3.*t2.*t3.*t5.*t10.*2.0-I3y.*I2z.*k3.*t2.*t3.*t5.*t10.*2.0+I2y.*I3z.*k3.*t3.*t5.*t12.*t13.*2.0+I3y.*I2z.*k3.*t3.*t5.*t12.*t13.*2.0+I2z.*I3z.*k3.*t2.*t3.*t5.*t10.*2.0-I2z.*I3z.*k3.*t3.*t5.*t12.*t13.*2.0-I3x.*dpsi.*t2.*t5.*t8.*t9.*t12+I3x.*dpsi.*t2.*t5.*t9.*t11.*t12-I3y.*dpsi.*t2.*t5.*t9.*t12.*t42+I3z.*dpsi.*t2.*t5.*t9.*t12.*t42+I3x.*I2y.*I3y.*dtheta.*t3.*t4.*t9.*2.0-I3x.*I2y.*I3z.*dtheta.*t3.*t4.*t9.*2.0+I3x.*I2y.*I3y.*dpsi.*t2.*t5.*t9.*t12.*2.0-I3x.*I2y.*I3z.*dpsi.*t2.*t5.*t9.*t12.*2.0+I3x.*I3y.*k2.*t2.*t4.*t5.*t12.*t15.*2.0+I3x.*I3y.*k2.*t2.*t5.*t9.*t12.*t14.*2.0-I3x.*I3z.*k2.*t2.*t4.*t5.*t12.*t15.*2.0-I3x.*I3z.*k2.*t2.*t5.*t9.*t12.*t14.*2.0).*(-1.0./2.0))./I3x];