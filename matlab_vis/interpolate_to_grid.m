function [xg, yg, Pgrid, Tgrid] = interpolate_to_grid(x, y, p, T, nx, ny)
% 把单元中心散点 (x,y) 的 p/T 插值到规则网格（nx×ny）
xmin = min(x); xmax = max(x);
ymin = min(y); ymax = max(y);
[xg, yg] = meshgrid(linspace(xmin, xmax, nx), linspace(ymin, ymax, ny));

Fp = scatteredInterpolant(x, y, p, 'natural', 'none');   % 外延 NaN
Ft = scatteredInterpolant(x, y, T, 'natural', 'none');
Pgrid = Fp(xg, yg);
Tgrid = Ft(xg, yg);
end
