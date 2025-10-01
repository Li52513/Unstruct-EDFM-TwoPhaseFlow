function plot_one_step(folder, step)
% 绘制单个时间步的压力/温度“云图”（规则网格插值）
% folder: 输出目录（例如 out_txt_CO2）
% step  : 步号（整数，如 1000）

[fname, data] = load_step_file(folder, step);

% 散点 -> 规则网格（便于成“云图”）
nx = 200; ny = 200;
[xg, yg, Pgrid, Tgrid] = interpolate_to_grid(data.cx, data.cy, data.p, data.T, nx, ny);

figure('Name', sprintf('Step %d (t=%.6g)', data.step, data.time), 'Color', 'w');

subplot(1,2,1);
imagesc(xg(1,:), yg(:,1), Pgrid); set(gca, 'YDir','normal'); axis image tight;
title(sprintf('Pressure (Pa) @ step %d', data.step));
xlabel('x'); ylabel('y'); colorbar;
hold on; plot(data.cx, data.cy, 'k.', 'MarkerSize', 4); % 可视化单元中心

subplot(1,2,2);
imagesc(xg(1,:), yg(:,1), Tgrid); set(gca, 'YDir','normal'); axis image tight;
title(sprintf('Temperature (K) @ step %d', data.step));
xlabel('x'); ylabel('y'); colorbar;
hold on; plot(data.cx, data.cy, 'k.', 'MarkerSize', 4);

sgtitle(sprintf('t = %.6g s', data.time));
disp(['Plotted: ' fname]);
end