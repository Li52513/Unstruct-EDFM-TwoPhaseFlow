function make_gif(folder, outGif, stepStart, stepEnd, stepStride, whichField)
% 生成 GIF 动画
% whichField: 'p' 或 'T'（默认 'p'）
if nargin < 6, whichField = 'p'; end
nx = 200; ny = 200;

[~, data0] = load_step_file(folder, stepStart);
[xg, yg, P0, T0] = interpolate_to_grid(data0.cx, data0.cy, data0.p, data0.T, nx, ny);
grid0 = strcmpi(whichField,'p') * P0 + strcmpi(whichField,'T') * T0;
cax = [min(grid0(:)), max(grid0(:))];  % 用第一帧定色标（也可扫全程 min/max）

figure('Color','w');
for s = stepStart:stepStride:stepEnd
    try
        [~, data] = load_step_file(folder, s);
        [xg, yg, Pgrid, Tgrid] = interpolate_to_grid(data.cx, data.cy, data.p, data.T, nx, ny);
        if strcmpi(whichField,'p')
            grid = Pgrid; tit = sprintf('Pressure (Pa) @ step %d', data.step);
        else
            grid = Tgrid; tit = sprintf('Temperature (K) @ step %d', data.step);
        end

        imagesc(xg(1,:), yg(:,1), grid); set(gca,'YDir','normal'); axis image tight; colorbar;
        if all(isfinite(cax)), caxis(cax); end
        title({tit, sprintf('t = %.6g s', data.time)});
        xlabel('x'); ylabel('y'); drawnow;

        % 写 GIF
        frame = getframe(gcf); [im, ~] = frame2im(frame); [A, cmap] = rgb2ind(im, 256);
        if s == stepStart
            imwrite(A, cmap, outGif, 'gif', 'LoopCount', inf, 'DelayTime', 0.12);
        else
            imwrite(A, cmap, outGif, 'gif', 'WriteMode', 'append', 'DelayTime', 0.12);
        end
    catch ME
        warning('Failed step %d: %s', s, ME.message);
    end
end
disp(['GIF saved: ' outGif]);
end
