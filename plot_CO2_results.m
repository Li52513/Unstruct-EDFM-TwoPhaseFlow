function plot_CO2_results()
% 入口（示例）
% 1) 单步压力/温度云图：
%    plot_one_step('D:\Yongwei\博士生涯\100-Research\110Code\111-2D_EDFM_FVM_CO2PlumingSystem\B-Code\2D-Unstr-Quadrilateral-EDFM\out_txt_CO2', 1000);
%
% 2) 批量浏览（每隔若干步）：
%    browse_many_steps('...\out_txt_CO2', 1, 1000, 50);   % 从 1 到 1000，每 50 步画一次
%
% 3) 生成动画 GIF（压力、温度各一份）：
%    make_gif('...\out_txt_CO2', 'pressure.gif', 1, 1000, 10, 'p');  % 每 10 步一帧
%    make_gif('...\out_txt_CO2', 'temperature.gif', 1, 1000, 10, 'T');

disp('示例函数已加载：plot_one_step, browse_many_steps, make_gif');
end


% ============== 单步云图 ==============
function plot_one_step(folder, step)
% folder: 输出目录（out_txt_CO2）
% step:   步号（整数），例如 1000
[fname, data] = load_step_file(folder, step);

% 散点插值到规则网格（更平滑的“云图”）
nx = 200; ny = 200;
[xg, yg, Pgrid, Tgrid] = interpolate_to_grid(data.cx, data.cy, data.p, data.T, nx, ny);

figure('Name', sprintf('Step %d (t=%.6g)', data.step, data.time), 'Color', 'w');

% 压力
subplot(1,2,1);
imagesc(xg(1,:), yg(:,1), Pgrid); 
set(gca, 'YDir','normal'); axis image tight;
title(sprintf('Pressure (Pa) @ step %d', data.step));
xlabel('x'); ylabel('y'); colorbar;

% 温度
subplot(1,2,2);
imagesc(xg(1,:), yg(:,1), Tgrid); 
set(gca, 'YDir','normal'); axis image tight;
title(sprintf('Temperature (K) @ step %d', data.step));
xlabel('x'); ylabel('y'); colorbar;

% 叠加单元中心（可选）
hold on; subplot(1,2,1); hold on; plot(data.cx, data.cy, 'k.', 'MarkerSize', 4);
subplot(1,2,2); hold on; plot(data.cx, data.cy, 'k.', 'MarkerSize', 4);

sgtitle(sprintf('t = %.6g s', data.time));
disp(['Plotted: ' fname]);
end


% ============== 批量浏览若干步 ==============
function browse_many_steps(folder, stepStart, stepEnd, stepStride)
if nargin < 4, stepStride = 10; end
for s = stepStart:stepStride:stepEnd
    try
        plot_one_step(folder, s);
        drawnow;
    catch ME
        warning('Failed step %d: %s', s, ME.message);
    end
end
end


% ============== 生成 GIF 动画（p 或 T） ==============
function make_gif(folder, outGif, stepStart, stepEnd, stepStride, whichField)
% whichField: 'p' 或 'T'
if nargin < 6, whichField = 'p'; end
nx = 200; ny = 200;

[fname0, data0] = load_step_file(folder, stepStart);
[xg, yg, grid0, ~] = interpolate_to_grid(data0.cx, data0.cy, data0.p, data0.T, nx, ny);

% 统一颜色范围（按全程 min/max 也行；这里按第一帧设定）
cax = [min(grid0(:)), max(grid0(:))];

figure('Color','w');
for s = stepStart:stepStride:stepEnd
    try
        [~, data] = load_step_file(folder, s);
        [xg, yg, Pgrid, Tgrid] = interpolate_to_grid(data.cx, data.cy, data.p, data.T, nx, ny);
        switch lower(whichField)
            case 'p', grid = Pgrid; tit = sprintf('Pressure (Pa) @ step %d', data.step);
            otherwise, grid = Tgrid; tit = sprintf('Temperature (K) @ step %d', data.step);
        end
        imagesc(xg(1,:), yg(:,1), grid);
        set(gca,'YDir','normal'); axis image tight; colorbar;
        title({tit, sprintf('t = %.6g s', data.time)});
        xlabel('x'); ylabel('y');
        if all(isfinite(cax)), caxis(cax); end
        drawnow;

        % 写入 GIF
        frame = getframe(gcf);
        [im, map] = frame2im(frame);
        [A, cmap] = rgb2ind(im, 256);
        if s == stepStart
            imwrite(A, cmap, outGif, 'gif', 'LoopCount', inf, 'DelayTime', 0.1);
        else
            imwrite(A, cmap, outGif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
        end
    catch ME
        warning('Failed step %d: %s', s, ME.message);
    end
end
disp(['GIF saved: ' outGif]);
end


% ============== 工具：读单步文件 ==============
function [fname, data] = load_step_file(folder, step)
% 文件名兼容：step_******.txt 或 pT_step_******.csv
cand = dir(fullfile(folder, sprintf('*%06d*.txt', step)));
if isempty(cand)
    cand = dir(fullfile(folder, sprintf('*%06d*.csv', step)));
end
if isempty(cand)
    error('Cannot find step file for step=%d in %s', step, folder);
end
fname = fullfile(folder, cand(1).name);

% 使用 textscan，忽略以 # 开头的两行注释
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open %s', fname); end
C = textscan(fid, '%f %f %f %f %f %f %f', 'CommentStyle', '#', 'MultipleDelimsAsOne', true);
fclose(fid);

% 打包
data.cell_id = C{1};
data.cx      = C{2};
data.cy      = C{3};
data.cz      = C{4};
data.vol     = C{5};
data.p       = C{6};
data.T       = C{7};

% 从首行提取元信息（可选）
data.step = step;
data.time = parse_time_from_header(fname);
end


% ============== 工具：散点插值到规则网格 ==============
function [xg, yg, Pgrid, Tgrid] = interpolate_to_grid(x, y, p, T, nx, ny)
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
[xg, yg] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));

Fp = scatteredInterpolant(x, y, p, 'natural', 'none');   % 外延设为 NaN
Ft = scatteredInterpolant(x, y, T, 'natural', 'none');
Pgrid = Fp(xg, yg);
Tgrid = Ft(xg, yg);
end


% ============== 工具：从文件头解析时间（可选） ==============
function t = parse_time_from_header(fname)
% 若首行形如：# step 1000  time 9.99999e+00
% 解析 time；失败则返回 NaN
fid = fopen(fname, 'r'); t = NaN;
if fid < 0, return; end
line = fgetl(fid); fclose(fid);
if ~ischar(line), return; end
tok = regexp(line, 'time\s+([0-9eE\+\-\.]+)', 'tokens', 'once');
if ~isempty(tok), t = str2double(tok{1}); end
end
