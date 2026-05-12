function tot_wave = chatGPT_animate_traveling_wave(x, y, spatial_freq, speed, direction, duration, fps)
%ANIMATE_TRAVELING_WAVE_SCATTER Animate 2D scatter plot of a traveling wave.
%
%   animate_traveling_wave_scatter(x, y, spatial_freq, speed, direction, duration, fps)
%   shows a 2D color-coded scatter plot of a sinusoidal traveling wave evolving over time.
%
%   Inputs:
%     x, y           - column vectors of coordinates (same length)
%     spatial_freq   - spatial frequency (cycles per unit distance)
%     speed          - wave speed (units per second)
%     direction      - 2-element vector [dx, dy] (does not need to be unit)
%     duration       - total animation time (in seconds)
%     fps            - frames per second

    % Normalize direction
    direction = direction / norm(direction);

    % Project positions onto wave direction
    pos_proj = x * direction(1) + y * direction(2);

    % Wave parameters
    k = 2 * pi * spatial_freq;
    omega = 2 * pi * spatial_freq * speed;

    % Time points
    t_vals = linspace(0, duration, round(duration * fps));

    tot_wave = [];

    % Create figure and initial plot
    % figure;
    wave = sin(k * pos_proj - omega * t_vals(1));
    tot_wave = [tot_wave; wave'];
    % h = scatter(x, y, 50, wave, 'filled');
    % axis equal;
    % 
    % cmap = blue_cyclic_colormap(256);
    % colormap(cmap)
    % 
    % colorbar;
    % clim([-1, 1]);
    % title(sprintf('Traveling Wave at t = %.2f s', t_vals(1)));
    % xlabel('x');
    % ylabel('y');

    % Animation loop
    for t = t_vals
        wave = sin(k * pos_proj - omega * t);
        tot_wave = [tot_wave; wave'];
        % set(h, 'CData', wave);
        % title(sprintf('Traveling Wave at t = %.2f s', t));
        % drawnow;
    end
end