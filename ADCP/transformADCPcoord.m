function T = transformADCPcoord(heading, pitch, roll)
% Function that creates transformation matrix to go from
% instrument-referenced to (magnetic) Earth-referenced
% coordinate system: 
    hh = pi * (heading - 90) / 180;
    pp = pi * pitch / 180;
    rr = pi * roll / 180;
 
    % Make heading matrix
    H = [ cos(hh) sin(hh) 0 ; ...
         -sin(hh) cos(hh) 0 ; ...
            0       0     1 ];
 
    % Make tilt matrix
    P = [cos(pp) -sin(pp)*sin(rr) -cos(rr)*sin(pp) ; ...
           0          cos(rr)          -sin(rr)    ; ...
         sin(pp)  sin(rr)*cos(pp)  cos(pp)*cos(rr) ];
 
    % Make resulting transformation matrix
    T = H*P;
 
end