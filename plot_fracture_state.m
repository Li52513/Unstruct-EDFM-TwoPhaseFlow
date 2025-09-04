function plot_fracture_state(outdir, step, field)
% 可视化裂缝 CSV 中某个标量场，按段着色（每段常色）
% outdir : e.g. 'out/fracture'
% step   : 整数，如 0
% field  : 列名，如 'pf_w','Sf_w','Tf'，或你以后新加的任意标量

if nargin < 3, error('usage: plot_fracture_state(outdir, step, field)'); end

% 背景网格（可选）
[nodes, faces] = load_mesh_();

% 裂缝 CSV
csv = fullfile(outdir, sprintf('state_fracture_step%05d.csv', step));
T = readtable(csv, 'Delimiter', ',', 'TextType','string');
if ~ismember(field, T.Properties.VariableNames)
    error('字段 "%s" 不在 %s 中。可用列：\n%s', field, csv, strjoin(T.Properties.VariableNames, ', '));
end
time = read_time_(outdir, step);

% 顶点/面/颜色
N = height(T);
V = [T.x1, T.y1; T.x2, T.y2];           % (2N)×2
F = [(1:2:2*N).', (2:2:2*N).'];         % N×2
Cseg = T.(field);                        % N×1
CperV = repelem(Cseg, 2);                % (2N)×1 —— 每个端点复制该段的颜色

% 绘图
figure('Color','w'); hold on; axis equal; box on;
title(sprintf('Fractures "%s" @ step %d, t=%g', field, step, time), 'Interpreter','none');

% 背景 Faces
if ~isempty(faces)
    id2row = id2row_map_(nodes(:,1));
    for i = 1:size(faces,1)
        r1 = id2row(faces(i,2)); r2 = id2row(faces(i,3));
        plot([nodes(r1,2),nodes(r2,2)], [nodes(r1,3),nodes(r2,3)], 'k-', 'LineWidth', 0.6);
    end
end

% 彩色线段（每段常色）
patch('Faces',F,'Vertices',V, ...
      'FaceColor','none', ...
      'EdgeColor','flat', ...               % 关键：按顶点色取边颜色
      'FaceVertexCData',CperV, ...          % 关键：长度=顶点数
      'LineWidth',2);

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
        warning('mesh_faces.txt 不存在或列数不足（需要: id n1 n2 ...）；将仅绘制裂缝。');
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
