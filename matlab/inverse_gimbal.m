function [ q ] = inverse_gimbal( R )
%INVERSE_GIMBAL Summary of this function goes here
%   Detailed explanation goes here
    theta = 0;
    phi = 0;
    psi = 0;

    USE_POSITIVE_PITCH = true;

    if ( ~(is_zero(r(1,3))) || ~(is_zero(r(2,3))))
        if (USE_POSITIVE_PITCH)
            theta = atan2(r(3,3), sqrt(1-r(3,3)^2));
            phi = atan2(r(1,3), r(2,3));
            psi = atan2(-r(3,1),r(3,2));
        else
            theta = atan2(r(3,3), -sqrt(1-r(3,3)^2));
            phi = atan2(-r(1,3), -r(2,3));
            psi = atan2(r(3,1), -r(3,2));
        end
    else
        if (r(3,3) > 0)
            sum_phi_psi = atan2(r(1,1),r(2,1));
            phi = 0; % Some convention?
            psi = sum_phi_psi - phi;
            theta = 0;
        else
            difference_phi_psi = atan2(-r(1,1), -r(1,2));
            phi = 0; % Some convention
            psi = phi - difference_phi_psi;
            theta = pi;
        end
    end
    
    q = [phi; theta; psi];
end

function [ zero ] = is_zero(number)

    DELTA = 1e-5;

    if (abs(number) <= DELTA)
        zero = true;
    else
        zero = false;
    end
end