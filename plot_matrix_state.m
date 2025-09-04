function plot_matrix_state(outdir_matrix, step, field, outdir_fracture)
% 可视化矩阵（基岩）CSV某个标量场为“云图”，用 节点三角化 + 单元中心→节点插值
% outdir_matrix : e.g. 'out/matrix'
% step          : 整数，如 0
% field         : 列名，如 'T','p_w','S_w','C_eff' 等
% outdir_fracture : 可选，缺省则尝试 'out/fracture'

if nargin < 3, error('usage: plot_matrix_state(outdir_matrix, step, field, [outdir_fracture])'); end
if nargin < 4, outdir_fracture = 'out/fracture'; end

% === 背景网格（节点/边） ===
[nodes, faces] = load_mesh_();

% === 读矩阵 CSV（单元中心值） ===
csvM = fullfile(outdir_matrix, sprintf('state_matrix_step%05d.csv', step));
TM = readtable(csvM, 'Delimiter', ',', 'TextType','string');
if ~ismember(field, TM.Properties.VariableNames)
    error('字段 "%s" 不在 %s 中。\n可用列：\n%s', field, csvM, strjoin(TM.Properties.VariableNames, ', '));
end
time = read_time_(outdir_matrix, step);

Cxy  = [TM.cx, TM.cy];      % 单元中心坐标
Cval = TM.(field);          % 单元中心标量

% === 1) 节点三角化（覆盖整个域边界） ===
Vxy  = nodes(:,2:3);                % 节点坐标 (x,y)
triN = delaunay(Vxy(:,1), Vxy(:,2));

% === 2) 把“单元中心值”插值到“节点” ===
% 'natural' 自然邻插值更适合非规则散点；外推用 'nearest'
Fint = scatteredInterpolant(Cxy(:,1), Cxy(:,2), Cval, 'natural', 'nearest');
Vval = Fint(Vxy(:,1), Vxy(:,2));    % 节点上的值

% === 3) 画云图（插值着色，覆盖整个域），再叠加网格线、裂缝 ===
figure('Color','w'); hold on; axis equal; box on;
title(sprintf('Matrix "%s" @ step %d, t=%g', field, step, time), 'Interpreter','none');

% 连续着色云图（节点色 + 插值）
patch('Faces',triN, 'Vertices',Vxy, ...
      'FaceColor','interp', 'EdgeColor','none', ...
      'FaceVertexCData', Vval);

% 叠加网格边线（薄黑线）
if ~isempty(faces)
    id2row = id2row_map_(nodes(:,1));
    for i = 1:size(faces,1)
        r1 = id2row(faces(i,2)); r2 = id2row(faces(i,3));
        plot([nodes(r1,2),nodes(r2,2)], [nodes(r1,3),nodes(r2,3)], 'k-', 'LineWidth', 0.6);
    end
end

% 叠加同 step 的裂缝线段（若存在）
csvF = fullfile(outdir_fracture, sprintf('state_fracture_step%05d.csv', step));
if exist(csvF,'file')
    TF = readtable(csvF, 'Delimiter', ',', 'TextType','string');
    for i = 1:height(TF)
        plot([TF.x1(i), TF.x2(i)], [TF.y1(i), TF.y2(i)], 'r-', 'LineWidth', 1.8);
    end
end

colormap(parula); colorbar; xlabel('x'); ylabel('y');
set(gca,'DataAspectRatio',[1 1 1]); grid on; hold off;

end

% ---------- helpers ----------
function [nodes, faces] = load_mesh_()
    nodes = readmatrix('mesh_nodes.txt','FileType','text','NumHeaderLines',1);
    faces = readmatrix('mesh_faces.txt','FileType','text','NumHeaderLines',1);
    if isempty(nodes) || size(nodes,2) < 3
        error('mesh_nodes.txt 不存在或列数不足（需要: id x y ...）');
    end
    if isempty(faces) || size(faces,2) < 3
        warning('mesh_faces.txt 不存在或列数不足（需要: id n1 n2 ...）；将仅绘制云图。');
        faces = [];
    end
end

function id2row = id2row_map_(ids)
    mn = min(ids); mx = max(ids);
    id2row = NaN(mx - mn + 1, 1);
    for r = 1:numel(ids), id2row(ids(r) - mn + 1) = r; end
end

function t = read_time_(outdir, step)
    t = NaN;
    f = fullfile(outdir,'times.csv');
    if exist(f,'file')
        TT = readtable(f,'Delimiter',',','ReadVariableNames',false);
        ix = find(TT.Var1 == step, 1, 'first');
        if ~isempty(ix), t = TT.Var2(ix); end
    end
    if isnan(t), t = 0; end
end