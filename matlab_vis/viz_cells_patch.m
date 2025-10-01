function viz_cells_patch(cmd, varargin)
% 方案C：Voronoi + 域裁剪 + polyshape 可视化
% 命令：
%   list        : 列出现有步号（自动识别 step_* 或 pT_step_*；txt / csv 都行）
%   plot_step   : 画某一步（压力+温度）
%   make_gif    : 批量做 GIF（自动发现现有步号，跳过缺失项）
%
% 例：
%   viz_cells_patch list 'D:\...\out_txt_CO2'
%   viz_cells_patch plot_step 'D:\...\out_txt_CO2' 1000
%   viz_cells_patch make_gif  'D:\...\out_txt_CO2' auto auto 'pt_CO2.gif'
%   viz_cells_patch make_gif  'D:\...\out_txt_CO2' 1 1000 'pt_CO2.gif'

switch lower(cmd)
    case 'list'
        outFolder = varargin{1};
        S = scan_step_files(outFolder);
        if isempty(S.steps), fprintf('No step files found under: %s\n', outFolder); return; end
        fprintf('Found %d files. First=%d  Last=%d\n', numel(S.steps), S.steps(1), S.steps(end));
        if numel(S.steps) <= 30, disp(S.steps(:).'); else, disp(S.steps(1:15).'); disp(' ...'); disp(S.steps(end-14:end).'); end
        fprintf('Sample file: %s\n', S.files{1});

    case 'plot_step'
        outFolder = varargin{1};
        step      = double(varargin{2});
        S = scan_step_files(outFolder);
        [ok, fn] = resolve_step_file(S, step);
        if ~ok, error('Step %d not found. Run: viz_cells_patch list ''%s''', step, outFolder); end
        data  = read_step_file(fn);
        cells = make_voronoi_cells(data);
        pmin = data.p_min; pmax = data.p_max; Tmin = data.T_min; Tmax = data.T_max;

        figure('Color','w','Position',[100 100 980 420]);
        subplot(1,2,1);
        title(sprintf('Pressure (Pa) @ step %d, t=%.2f s', step, data.t), 'FontWeight','bold');
        draw_cells(cells, data.p, [pmin pmax]); xlabel('x'); ylabel('y'); axis image; box on; cb=colorbar; cb.Label.String='Pa';

        subplot(1,2,2);
        title(sprintf('Temperature (K) @ step %d, t=%.2f s', step, data.t), 'FontWeight','bold');
        draw_cells(cells, data.T, [Tmin Tmax]); xlabel('x'); ylabel('y'); axis image; box on; cb=colorbar; cb.Label.String='K';

    case 'make_gif'
        outFolder = varargin{1};
        s0_in     = varargin{2};
        s1_in     = varargin{3};
        gifPath   = varargin{4};

        S = scan_step_files(outFolder);
        if isempty(S.steps), error('No step files found in %s', outFolder); end
        if (ischar(s0_in) && strcmpi(s0_in,'auto')) || (isstring(s0_in) && strcmpi(s0_in,'auto')), s0=S.steps(1); else, s0=double(s0_in); end
        if (ischar(s1_in) && strcmpi(s1_in,'auto')) || (isstring(s1_in) && strcmpi(s1_in,'auto')), s1=S.steps(end); else, s1=double(s1_in); end

        % 仅保留区间内且确实存在的步
        mask = (S.steps >= s0) & (S.steps <= s1);
        steps = S.steps(mask); files = S.files(mask);
        if isempty(steps), error('No existing steps in [%d,%d]. Try "list".', s0, s1); end

        % 用第一帧确定网格 & 色标
        d0    = read_step_file(files{1});
        cells = make_voronoi_cells(d0);
        pmin  = d0.p_min; pmax = d0.p_max; Tmin = d0.T_min; Tmax = d0.T_max;

        fig = figure('Color','w','Position',[100 100 980 420]);
        for i=1:numel(steps)
            di = read_step_file(files{i});
            clf(fig);
            subplot(1,2,1);
            title(sprintf('Pressure (Pa) @ step %d, t=%.2f s', steps(i), di.t), 'FontWeight','bold');
            draw_cells(cells, di.p, [pmin pmax]); xlabel('x'); ylabel('y'); axis image; box on; cb=colorbar; cb.Label.String='Pa';

            subplot(1,2,2);
            title(sprintf('Temperature (K) @ step %d, t=%.2f s', steps(i), di.t), 'FontWeight','bold');
            draw_cells(cells, di.T, [Tmin Tmax]); xlabel('x'); ylabel('y'); axis image; box on; cb=colorbar; cb.Label.String='K';

            drawnow;
            [im, map] = rgb2ind(frame2im(getframe(fig)), 256, 'nodither');
            if i==1, imwrite(im,map,gifPath,'gif','LoopCount',Inf,'DelayTime',0.05);
            else,    imwrite(im,map,gifPath,'gif','WriteMode','append','DelayTime',0.05);
            end
        end
        fprintf('[GIF] wrote %s  (frames: %d from [%d..%d], missing steps auto-skipped)\n', gifPath, numel(steps), s0, s1);

    otherwise
        error('Unknown cmd. Use: list | plot_step | make_gif');
end
end

%==================== 扫描目录，找出现有步与文件 ====================%
function S = scan_step_files(outFolder)
d = [ dir(fullfile(outFolder,'step_*.txt')); ...
      dir(fullfile(outFolder,'pT_step_*.txt')); ...
      dir(fullfile(outFolder,'step_*.csv')); ...
      dir(fullfile(outFolder,'pT_step_*.csv')) ];

% 如果还没找到，兜底搜 *.txt / *.csv，并尝试解析文件名里的步号
if isempty(d)
    d = [ dir(fullfile(outFolder,'*.txt')); dir(fullfile(outFolder,'*.csv')) ];
end

steps = []; files = {};
for i=1:numel(d)
    name = d(i).name;
    m = regexp(name, '(?:^|_)step_(\d+)\.(txt|csv)$', 'tokens', 'once');
    if isempty(m)
        m = regexp(name, '(?:^|_)pT_step_(\d+)\.(txt|csv)$', 'tokens', 'once');
    end
    if isempty(m)
        % 最后兜底：抓取文件名里最后一个数字串
        m2 = regexp(name, '(\d+)(?=\.(txt|csv)$)', 'tokens', 'once');
        if ~isempty(m2), m = {m2{1} 'txt'}; end
    end
    if ~isempty(m)
        steps(end+1,1) = str2double(m{1}); %#ok<AGROW>
        files{end+1,1} = fullfile(outFolder, name); %#ok<AGROW>
    end
end

[steps, idx] = sort(steps);
files = files(idx);
S.steps = steps;
S.files = files;
end

%==================== 按步号解析文件路径（若存在） ====================%
function [ok, fn] = resolve_step_file(S, step)
idx = find(S.steps==step, 1, 'first');
ok = ~isempty(idx);
fn = '';
if ok, fn = S.files{idx}; end
end

%==================== 读取一步文件（txt/csv 均可） ====================%
function data = read_step_file(fn)
% 文件第一行可能含有 "# step xxx  time xxx"；主体是 7 列（空格或逗号分隔）
opts = detectImportOptions(fn,'FileType','text','CommentStyle','#');
T = readtable(fn, opts, 'ReadVariableNames', false);
if size(T,2) < 7
    error('Unexpected columns in %s (need >=7: id,cx,cy,cz,vol,p,T).', fn);
end
data.id = T{:,1}; data.cx = T{:,2}; data.cy = T{:,3};
data.cz = T{:,4}; data.vol= T{:,5}; data.p  = T{:,6}; data.T  = T{:,7};
data.t  = try_read_time_from_header(fn);
data.p_min = min(data.p); data.p_max = max(data.p);
data.T_min = min(data.T); data.T_max = max(data.T);
end

function t = try_read_time_from_header(fn)
t = NaN; fid=fopen(fn,'r'); if fid<0, return; end
l1=fgetl(fid); fclose(fid);
tok = regexp(l1, 'time\s+([Ee0-9\.\+\-]+)', 'tokens', 'once');
if ~isempty(tok), t = str2double(tok{1}); end
end

%==================== Voronoi 单元构建 ====================%
function cells = make_voronoi_cells(data)
x = data.cx(:); y = data.cy(:);
pad = 1e-9; xmin=min(x)-pad; xmax=max(x)+pad; ymin=min(y)-pad; ymax=max(y)+pad;
dom  = polyshape([xmin xmax xmax xmin],[ymin ymin ymax ymax]);

[V,C] = voronoin([x y]);
n = numel(C); P = repmat(polyshape,n,1);
for i=1:n
    ci = C{i}; ci = ci(~isinf(ci));
    if isempty(ci), continue; end
    pg = polyshape(V(ci,1),V(ci,2));
    pg = intersect(pg, dom); pg = rmholes(pg);
    P(i) = pg;
end
cells.dom=dom; cells.polys=P; cells.xmin=xmin; cells.xmax=xmax; cells.ymin=ymin; cells.ymax=ymax; cells.n=n;
end

%==================== 画多边形着色 ====================%
function draw_cells(cells, values, clim)
P=cells.polys; n=cells.n; vmin=clim(1); vmax=clim(2); cmap=parula(256);
hold on;
for i=1:n
    if isempty(P(i).Vertices), continue; end
    t = (values(i)-vmin)/max(vmax-vmin, eps); t = max(0,min(1,t));
    c = cmap(1+floor(t*(size(cmap,1)-1)),:);
    ph = plot(P(i)); set(ph,'FaceColor',c,'EdgeColor','none'); % 想看网格线，把 'none' 改 'k'
end
colormap(cmap); caxis([vmin vmax]);
xlim([cells.xmin cells.xmax]); ylim([cells.ymin cells.ymax]);
end
