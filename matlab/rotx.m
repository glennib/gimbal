function [ R ] = rotx( theta )
%ROTX Summary of this function goes here
%   Detailed explanation goes here

R = [
    1 0 0;
    0 cos(theta) -sin(theta);
    0 sin(theta) cos(theta)
    ];

end

