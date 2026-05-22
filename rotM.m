function SO = rotM(s)
% Rotation Matrix (1-2-3) or (x-y-z)
% RotM = @(x)[cos(x(2))*cos(x(3)) cos(x(1))*sin(x(3))+sin(x(1))*sin(x(2))*cos(x(3)) sin(x(1))*sin(x(3))-cos(x(1))*sin(x(2))*cos(x(3));
%         -cos(x(2))*sin(x(3)) cos(x(1))*cos(x(3))-sin(x(1))*sin(x(2))*sin(x(3)) cos(x(1))*sin(x(2))*sin(x(3))+sin(x(1))*cos(x(3));
%         sin(x(2)) -sin(x(1))*cos(x(2)) cos(x(1))*cos(x(2))];

    roll = s(1);
    pitch = s(2);
    yaw = s(3);

    c_roll = cos(roll);
    c_pitch = cos(pitch);
    c_yaw = cos(yaw);
    s_roll = sin(roll);
    s_pitch = sin(pitch);
    s_yaw = sin(yaw);
    t_roll = tan(roll);
    t_pitch = tan(pitch);
    t_yaw = tan(yaw);
   
    SO = [c_pitch*c_yaw, s_roll*s_pitch*c_yaw-c_roll*s_yaw, c_roll*s_pitch*c_yaw+s_roll*s_yaw;
          c_pitch*s_yaw, s_roll*s_pitch*s_yaw+c_roll*c_yaw, c_roll*s_pitch*s_yaw-s_pitch*c_yaw;
          -s_pitch, s_roll*c_pitch, c_roll*c_pitch]; % Rotation Matrix

    SO = rotMz(yaw)*rotMy(pitch)*rotMx(roll);

    % SO = [c_pitch*c_roll, s_yaw*s_pitch*c_roll-c_yaw*s_roll, c_yaw*s_pitch*c_roll+s_yaw*s_roll;
    %       c_pitch*s_roll, s_yaw*s_pitch*s_roll+c_yaw*c_roll, c_yaw*s_pitch*s_roll-s_yaw*c_roll;
    %       -s_pitch, s_yaw*c_pitch, c_yaw*c_pitch];
end