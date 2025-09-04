%% visualize_EDFM.m ──────────────────────────────────────────────────────────
% 可视化 EDFM 网格、裂缝段、裂缝-裂缝交点（含 Ghost Cell）

close all; clear; clc;

%% ----------------- 读取文件 -----------------
nodes = readmatrix('mesh_nodes.txt'                 , 'FileType','text','NumHeaderLines',1);
faces = readmatrix('mesh_faces.txt'                 , 'FileType','text','NumHeaderLines',1);
cells = readmatrix('mesh_cells.txt'                 , 'FileType','text','NumHeaderLines',1);
segs  = readmatrix('fractures_fractureSegments.txt' , 'FileType','text','NumHeaderLines',1);
ffpt  = readmatrix('fractures_ffPoints.txt'         , 'FileType','text','NumHeaderLines',1);

%% ----------------- 基本维度 -----------------
Nnode = size(nodes,1);
Nface = size(faces,1);
Ncell = size(cells,1);

%% ----------------- 解析 G0（稳健版） -----------------
% 支持 "ghostStartIndex (0-based) = 42"、"42" 等格式
G0 = parse_ghost_index('mesh_info.txt', Ncell);

%% ----------------- 节点 ID → 行号 映射 -----------------
nid  = nodes(:,1);
nid_min = min(nid);  nid_max = max(nid);
id2row = nan(nid_max - nid_min + 1, 1);
for r = 1:Nnode
    id2row( nid(r) - nid_min + 1 ) = r;
end
id2row_fun = @(id) assert_return(id2row(id - nid_min + 1), ...
    sprintf('节点ID %d 在 nodes 表中不存在！', id));

%% ----------------- 区分 Real / Ghost Cells（优先用显式标志） -----------------
hasGhostFlag = size(cells,2) >= 6 && all(ismember(unique(cells(:,end)),[0 1]));
if hasGhostFlag
    isGhost = logical(cells(:,end));
    isReal  = ~isGhost;
else
    cid = cells(:,1);                 % 假定第1列为 cell_id
    zeroBased = (min(cid) == 0);      % 判断 0-based 还是 1-based
    if zeroBased
        isReal  = (cid <  G0);        % 0-based: real ids ∈ [0, G0-1]
    else
        isReal  = (cid <= G0);        % 1-based: real ids ∈ [1, G0]
    end
    isGhost = ~isReal;
end

% Sanity Check
if ~any(isReal) || ~any(isGhost)
    warning('Real/Ghost 判定异常：real=%d, ghost=%d, G0=%d（请核对 mesh_info.txt 与 mesh_cells.txt 约定）',...
        nnz(isReal), nnz(isGhost), G0);
end
realCells  = cells(isReal,  :);
ghostCells = cells(isGhost, :);

%% ----------------- 画图参数 -----------------
psz_node    = 18;  psz_cell = 36; psz_ghost = 36; psz_ffpt = 50;
f_line_w    = 0.7; seg_line_w = 1.5;
col_seg     = [0.80,0.0,0.0];
col_ffpt    = [0.0,0.6,0.0];
dx_node     = .01; dy_node = .01;
dx_cell     = .01; dy_cell = .01;

%% ----------------- 绘图 -----------------
figure('Color','w','Units','normalized','Position',[0.1,0.1,0.8,0.8]);
hold on; axis equal off;
title('EDFM mesh & fractures (with Ghost Cells)','FontWeight','bold');

% 1) Nodes
scatter(nodes(:,2), nodes(:,3), psz_node, 'k', 'filled');
text(nodes(:,2)+dx_node, nodes(:,3)+dy_node, string(nodes(:,1)), ...
     'FontSize',10, 'Color',[.25 .25 .25], 'FontWeight','bold');

% 2) Faces（用 ID → 行号 映射取端点）
for i = 1:Nface
    n1_id = faces(i,2);  n2_id = faces(i,3);
    r1 = id2row_fun(n1_id);  r2 = id2row_fun(n2_id);
    plot([nodes(r1,2), nodes(r2,2)], [nodes(r1,3), nodes(r2,3)], 'k-', 'LineWidth', f_line_w);
    if size(faces,2) >= 5
        text(faces(i,4), faces(i,5)+0.02, "F"+string(faces(i,1)), 'FontSize',10, 'Color','r');
    end
end

% 3) Cell centres：Real → 蓝点，Ghost → 红叉
scatter(realCells(:,2),  realCells(:,3),  psz_cell,  'b', 'filled');
scatter(ghostCells(:,2), ghostCells(:,3), psz_ghost, 'r', 'x', 'LineWidth',1.2);
text(realCells(:,2)+dx_cell, realCells(:,3)+dy_cell, string(realCells(:,1)), ...
     'FontSize',10, 'Color','b', 'FontWeight','bold');

% 4) Fracture segments
if ~isempty(segs)
    for i = 1:size(segs,1)
        x1 = segs(i,3); y1 = segs(i,4);
        x2 = segs(i,5); y2 = segs(i,6);
        plot([x1,x2],[y1,y2], '-', 'LineWidth', seg_line_w, 'Color', col_seg);
        if size(segs,2) >= 8
            txt = segs(i,1) + "." + segs(i,2);
            text(segs(i,7), segs(i,8), txt, 'FontSize',10, 'Color',[.1 .5 .1]);
        end
    end
end

% 5) Fracture–Fracture intersection points
if ~isempty(ffpt)
    scatter(ffpt(:,2), ffpt(:,3), psz_ffpt, 'o', 'MarkerEdgeColor','k', ...
            'MarkerFaceColor',col_ffpt, 'LineWidth',1.2);
    text(ffpt(:,2)+0.01, ffpt(:,3)+0.01, "X"+string(ffpt(:,1)), ...
         'FontSize',10, 'Color',col_ffpt, 'FontWeight','bold');
end

% 6) 图例与统计
legend({'Nodes','Faces','Real cell centres','Ghost cell centres','Frac segments','Frac–Frac pts'}, ...
       'Location','bestoutside');
fprintf('[viz] nodes=%d, faces=%d, cells=%d, real=%d, ghost=%d, G0=%d\n', ...
    Nnode, Nface, Ncell, nnz(isReal), nnz(isGhost), G0);
hold off;

%% ======================== 辅助函数 ========================
function x = assert_return(val, errmsg)
    if isempty(val) || (isnumeric(val) && any(isnan(val)))
        error(errmsg);
    end
    x = val;
end

function G0 = parse_ghost_index(fname, Ncells)
% 从文件中抓取所有整数；优先取落在 [1..Ncells] 区间内的“最后一个”，否则取最大值
    raw = fileread(fname);
    tok = regexp(raw, '[-+]?\d+', 'match');
    assert(~isempty(tok), 'mesh_info.txt 中未找到整数！');
    vals = str2double(tok);
    cand = vals(vals >= 1 & vals <= Ncells);
    if ~isempty(cand)
        G0 = cand(end);
    else
        G0 = max(vals);
    end
end
