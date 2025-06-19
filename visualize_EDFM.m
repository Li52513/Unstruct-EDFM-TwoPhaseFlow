
%% visualize_EDFM.m ──────────────────────────────────────────────────────────
% 可视化 EDFM 网格、裂缝段与裂缝-裂缝交点（含 Ghost Cell）
% 依赖文件（均跳过表头第一行）：
%   mesh_nodes.txt
%   mesh_faces.txt
%   mesh_cells.txt
%   mesh_info.txt       ← 仅含一行数字：ghostStartIndex（0-based 下标）
%   fractures_fractureSegments.txt
%   fractures_ffPoints.txt

close all; clear; clc;

%% ----------- 读取 txt 文件 -----------
nodes = readmatrix('mesh_nodes.txt'                 , 'FileType','text','NumHeaderLines',1);
faces = readmatrix('mesh_faces.txt'                 , 'FileType','text','NumHeaderLines',1);
cells = readmatrix('mesh_cells.txt'                 , 'FileType','text','NumHeaderLines',1);
segs  = readmatrix('fractures_fractureSegments.txt' , 'FileType','text','NumHeaderLines',1);
ffpt  = readmatrix('fractures_ffPoints.txt'         , 'FileType','text','NumHeaderLines',1);

% 读取 ghostStartIndex（C++ 导出，0-based 向量下标）
G0 = readmatrix('mesh_info.txt'                     , 'FileType','text','NumHeaderLines',0);

%% ----------- 区分真实单元 vs. 幽灵单元 -----------
Ncells  = size(cells,1);
rowIdx  = (1:Ncells)';            % 行号 1 对应 cells(1,:), 依次类推
isReal  = rowIdx <= G0;           % 行号 ≤ G0 是原始单元
realCells  = cells(isReal,   :);
ghostCells = cells(~isReal,  :);

%% ----------- 可视化参数设置 -----------
psz_node    = 18;       % 节点散点大小
psz_cell    = 36;       % 真单元中心散点大小
psz_ghost   = 36;       % 幽灵单元中心散点大小
psz_ffpt    = 50;       % 裂缝-裂缝交点散点大小
f_line_w    = 0.7;      % 网格边线宽度
seg_line_w  = 1.5;      % 裂缝段线宽
col_seg     = [0.80,0.0,0.0];
col_ffpt    = [0.0,0.6,0.0];
dx_node     = .01; dy_node = .01;
dx_cell     = .01; dy_cell = .01;

%% ----------- 绘图 -----------
figure('Color','w','Units','normalized','Position',[0.1,0.1,0.8,0.8]);
hold on; axis equal off;
title('EDFM mesh & fractures (with Ghost Cells)','FontWeight','bold');

% 1) 绘制 Nodes
scatter(nodes(:,2), nodes(:,3), psz_node, 'k', 'filled');
text(nodes(:,2)+dx_node, nodes(:,3)+dy_node, string(nodes(:,1)), ...
     'FontSize',10, 'Color',[.25 .25 .25], 'FontWeight','bold');

% 2) 绘制 Faces
for i = 1:size(faces,1)
    n1 = faces(i,2); 
    n2 = faces(i,3);
    plot([nodes(n1,2), nodes(n2,2)], [nodes(n1,3), nodes(n2,3)], ...
         'k-', 'LineWidth', f_line_w);
    text(faces(i,4), faces(i,5)+0.02, "F"+string(faces(i,1)), ...
         'FontSize',10, 'Color','r');
end

% 3) 绘制 Cell centres：真单元 → 蓝点，幽灵单元 → 红叉
scatter(realCells(:,2),  realCells(:,3),  psz_cell,  'b', 'filled');
scatter(ghostCells(:,2), ghostCells(:,3), psz_ghost, 'r', 'x', 'LineWidth',1.2);
% 只给真单元添加标签
text(realCells(:,2)+dx_cell, realCells(:,3)+dy_cell, string(realCells(:,1)), ...
     'FontSize',10, 'Color','b', 'FontWeight','bold');

% 4) 绘制 Fracture segments
for i = 1:size(segs,1)
    x1 = segs(i,3); y1 = segs(i,4);
    x2 = segs(i,5); y2 = segs(i,6);
    plot([x1,x2],[y1,y2], '-', 'LineWidth', seg_line_w, 'Color', col_seg);
    txt = segs(i,1) + "." + segs(i,2);       % 格式：fid.segID
    text(segs(i,7), segs(i,8), txt, 'FontSize',10, 'Color',[.1 .5 .1]);
end

% 5) 绘制 Fracture–Fracture intersection points
scatter(ffpt(:,2), ffpt(:,3), psz_ffpt, 'o', ...
        'MarkerEdgeColor','k', 'MarkerFaceColor',col_ffpt, 'LineWidth',1.2);
text(ffpt(:,2)+0.01, ffpt(:,3)+0.01, "X"+string(ffpt(:,1)), ...
     'FontSize',10, 'Color',col_ffpt, 'FontWeight','bold');

% 6) 添加图例
legend({ ...
    'Nodes', ...
    'Faces', ...
    'Real cell centres', ...
    'Ghost cell centres', ...
    'Frac segments', ...
    'Frac–Frac pts' ...
    }, 'Location','bestoutside');

hold off;
