function R = R_E321(Euler_angles)

% Euler angles 3-2-1 sequence (yaw, pitch and roll angles), respective sines and cosines
yaw = Euler_angles(1); sy = sin(yaw); cy = cos(yaw); 
pitch = Euler_angles(2); cp = cos(pitch); sp = sin(pitch);
roll = Euler_angles(3); cr = cos(roll); sr = sin(roll);

% Rotation tensor that brings the reference basis to the final basis
R = [cp*cy, cy*sr*sp - cr*sy, sr*sy + cr*cy*sp
     cp*sy, cr*cy + sr*sp*sy, cr*sp*sy - cy*sr
       -sp,            cp*sr,            cr*cp];

end