function [robot_out] = modify_robot(rob_in, params, num_bodies)
    robot_out = rob_in;

    % other
    robot_out.DataFormat = 'col';
    robot_out.Gravity = [0; 0; -9.81];

    % check initialized
    % if init_in
    %     return;
    % end

    for i = 1:num_bodies
        % Mass
        robot_out.Bodies{i}.Mass = params(i,1);

        % COM
        if size(params(i, :), 2) >= 4
            robot_out.Bodies{i}.CenterOfMass = params(i,2:4);
        end

        % Inertia
        if size(params(i, :), 2) >= 10
            robot_out.Bodies{i}.Inertia = ...
                parallel_axis(params(i,5:10), params(i, 1), params(i, 2:4));
        end

    end
end
