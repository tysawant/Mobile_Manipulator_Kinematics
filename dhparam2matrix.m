%% Question 3(a)
function T = dhparam2matrix(theta, d, a, alpha)
    Rotz = [cos(theta) -sin(theta) 0 0; sin(theta) cos(theta) 0 0; 0 0 1 0; 0 0 0 1];
    Transz = [1 0 0 0; 0 1 0 0; 0 0 1 d; 0 0 0 1];
    Transx = [1 0 0 a; 0 1 0 0; 0 0 1 0; 0 0 0 1];
    Rotx = [1 0 0 0; 0 cos(alpha) -sin(alpha) 0; 0 sin(alpha) cos(alpha) 0; 0 0 0 1];
    T = Rotz*Transz*Transx*Rotx;
end