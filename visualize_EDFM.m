% %% visualize_EDFM.m  ─────────────────────────────────────────────────────────
% % 可视化 EDFM 网格、裂缝段、裂缝-裂缝交点（无 ghost；严格按列；锁定视窗）
% 
% close all; clear; clc;
% 
% %% ----------------- 读取文件（严格按列） -----------------
% % mesh_* 由 C++ 导出：
% %   mesh_nodes.txt:  id x y z
% %   mesh_faces.txt:  id n1 n2 mx my mz
% %   mesh_cells.txt:  id cx cy cz
% nodes = readmatrix('mesh_nodes.txt' , 'FileType','text','NumHeaderLines',1);
% faces = readmatrix('mesh_faces.txt' , 'FileType','text','NumHeaderLines',1);
% cells = readmatrix('mesh_cells.txt' , 'FileType','text','NumHeaderLines',1);
% 
% % fractures_*（若无则跳过对应图层）：
% %   fractures_fractureSegments.txt: fid segid x1 y1 x2 y2 mx my cellID
% %   fractures_ffPoints.txt        : id  x  y   （若无 id 则自动生成）
% segs  = []; ffpt = [];
% if exist('fractures_fractureSegments.txt','file')
%     segs = readmatrix('fractures_fractureSegments.txt','FileType','text','NumHeaderLines',1);
% end
% if exist('fractures_ffPoints.txt','file')
%     tmp  = readmatrix('fractures_ffPoints.txt','FileType','text','NumHeaderLines',1);
%     if ~isempty(tmp)
%         if size(tmp,2) >= 3, ffpt = tmp(:,1:3);         % id x y
%         else,                ffpt = [(1:size(tmp,1))' tmp(:,1:2)]; end
%     end
% end
% 
% %% ----------------- 基本信息与域范围 -----------------
% Nnode = size(nodes,1); Nface = size(faces,1); Ncell = size(cells,1);
% xmin = min(nodes(:,2)); xmax = max(nodes(:,2));
% ymin = min(nodes(:,3)); ymax = max(nodes(:,3));
% % 视窗留一点边距
% mx = 0.02*(xmax-xmin + eps); my = 0.02*(ymax-ymin + eps);
% xlim_box = [xmin-mx, xmax+mx]; ylim_box = [ymin-my, ymax+my];
% 
% %% ----------------- 节点 ID → 行号 映射 -----------------
% nid     = nodes(:,1);
% nid_min = min(nid);  nid_max = max(nid);
% id2row  = nan(nid_max - nid_min + 1, 1);
% for r = 1:Nnode
%     id2row( nid(r) - nid_min + 1 ) = r;
% end
% id2row_fun = @(id) assert_return(id2row(id - nid_min + 1), ...
%     sprintf('节点ID %d 在 nodes 表中不存在！', id));
% 
% %% ----------------- 画图参数/开关 -----------------
% showNodeLabel   = true;
% showFaceLabel   = true;
% showCellLabel   = true;
% showSegLabel    = true;
% showFFPts       = true;
% 
% psz_node  = 18; psz_cell = 36; psz_ffpt = 50;
% f_line_w  = 0.8; seg_line_w = 1.5;
% col_seg   = [0.80,0.0,0.0];
% col_ffpt  = [0.0,0.6,0.0];
% dx_node   = .01; dy_node = .01;
% dx_cell   = .01; dy_cell = .01;
% 
% %% ----------------- 绘图 -----------------
% figure('Color','w','Units','normalized','Position',[0.08,0.08,0.82,0.82]);
% hold on; axis equal; box on;
% title('EDFM mesh & fractures (no ghost)','FontWeight','bold');
% 
% % 1) Faces（用 ID → 行号 映射取端点）
% for i = 1:Nface
%     n1_id = faces(i,2);  n2_id = faces(i,3);
%     r1 = id2row_fun(n1_id);  r2 = id2row_fun(n2_id);
%     plot([nodes(r1,2), nodes(r2,2)], [nodes(r1,3), nodes(r2,3)], ...
%          'k-', 'LineWidth', f_line_w);
%     if showFaceLabel && size(faces,2) >= 5
%         text(faces(i,4), faces(i,5)+0.02, "F"+string(faces(i,1)), ...
%              'FontSize',10, 'Color','r');
%     end
% end
% 
% % 2) Nodes
% scatter(nodes(:,2), nodes(:,3), psz_node, 'k', 'filled');
% if showNodeLabel
%     text(nodes(:,2)+dx_node, nodes(:,3)+dy_node, string(nodes(:,1)), ...
%          'FontSize',10, 'Color',[.25 .25 .25], 'FontWeight','bold');
% end
% 
% % 3) Cell centres（全部真实单元）
% scatter(cells(:,2), cells(:,3), psz_cell, 'b', 'filled');
% if showCellLabel
%     text(cells(:,2)+dx_cell, cells(:,3)+dy_cell, string(cells(:,1)), ...
%          'FontSize',10, 'Color','b', 'FontWeight','bold');
% end
% 
% % 4) Fracture segments：严格用 3:6 列，标签放 7:8 列
% if ~isempty(segs)
%     x1 = segs(:,3); y1 = segs(:,4);
%     x2 = segs(:,5); y2 = segs(:,6);
%     mxs = segs(:,7); mys = segs(:,8);        % 标注点（可不在中点）
%     fid = segs(:,1); sid = segs(:,2);
% 
%     % 可选：快速体检
%     fprintf('[seg check] x∈[%.3g, %.3g], y∈[%.3g, %.3g]\n', ...
%         min([x1;x2]), max([x1;x2]), min([y1;y2]), max([y1;y2]));
% 
%     % 逐段绘制
%     for k = 1:numel(x1)
%         plot([x1(k),x2(k)], [y1(k),y2(k)], '-', 'LineWidth', seg_line_w, 'Color', col_seg);
%         if showSegLabel
%             text(mxs(k), mys(k), sprintf('%d.%d', fid(k), sid(k)), ...
%                  'FontSize',10, 'Color',[.1 .5 .1]);
%         end
%     end
% end
% 
% % 5) Fracture–Fracture intersection points
% if showFFPts && ~isempty(ffpt)
%     scatter(ffpt(:,2), ffpt(:,3), psz_ffpt, 'o', ...
%             'MarkerEdgeColor','k', 'MarkerFaceColor',col_ffpt, 'LineWidth',1.2);
%     if size(ffpt,2) >= 3
%         text(ffpt(:,2)+0.01, ffpt(:,3)+0.01, "X"+string(ffpt(:,1)), ...
%              'FontSize',10, 'Color',col_ffpt, 'FontWeight','bold');
%     end
% end
% 
% % 6) 锁定视窗到域框
% xlim(xlim_box); ylim(ylim_box);
% 
% % 7) 图例与统计
% legendItems = {'Faces','Nodes','Cell centres'};
% if ~isempty(segs), legendItems{end+1} = 'Frac segments'; end
% if showFFPts && ~isempty(ffpt), legendItems{end+1} = 'Frac–Frac pts'; end
% legend(legendItems, 'Location','bestoutside');
% 
% fprintf('[viz] nodes=%d, faces=%d, cells=%d (no ghost). ', Nnode, Nface, Ncell);
% if ~isempty(segs), fprintf('segments=%d. ', size(segs,1)); end
% if ~isempty(ffpt), fprintf('ff-pts=%d. ', size(ffpt,1));  end
% fprintf('\n');
% hold off;
% 
% %% ======================== 辅助函数 ========================
% function x = assert_return(val, errmsg)
%     if isempty(val) || (isnumeric(val) && any(isnan(val)))
%         error(errmsg);
%     end
%     x = val;
% end



%% visualize_EDFM.m  ─────────────────────────────────────────────────────────
% 可视化 EDFM 网格、裂缝段、裂缝-裂缝交点
% 特点：无 ghost；严格按列；锁定视窗；裂缝段按宿主 cellID 上色

close all; clear; clc;

%% ----------------- 读取文件（严格按列） -----------------
% mesh_* 由 C++ 导出：
%   mesh_nodes.txt:  id x y z
%   mesh_faces.txt:  id n1 n2 mx my mz
%   mesh_cells.txt:  id cx cy cz
nodes = readmatrix('mesh_nodes.txt' , 'FileType','text','NumHeaderLines',1);
faces = readmatrix('mesh_faces.txt' , 'FileType','text','NumHeaderLines',1);
cells = readmatrix('mesh_cells.txt' , 'FileType','text','NumHeaderLines',1);

% fractures_*（若无则自动跳过）：
%   fractures_fractureSegments.txt: fid segid x1 y1 x2 y2 mx my cellID
%   fractures_ffPoints.txt        : id x y fracA segA fracB segB
segs  = [];
ffpt  = [];
if exist('fractures_fractureSegments.txt','file')
    segs = readmatrix('fractures_fractureSegments.txt','FileType','text','NumHeaderLines',1);
end
if exist('fractures_ffPoints.txt','file')
    tmp  = readmatrix('fractures_ffPoints.txt','FileType','text','NumHeaderLines',1);
    if ~isempty(tmp)
        if size(tmp,2) >= 3
            ffpt = tmp(:,1:3);                 % 只用 id/x/y
        end
    end
end

%% ----------------- 基本信息与域范围 -----------------
Nnode = size(nodes,1); Nface = size(faces,1); Ncell = size(cells,1);
xmin = min(nodes(:,2)); xmax = max(nodes(:,2));
ymin = min(nodes(:,3)); ymax = max(nodes(:,3));
% 视窗留 2% 边距
mx = 0.02*(xmax-xmin + eps); my = 0.02*(ymax-ymin + eps);
xlim_box = [xmin-mx, xmax+mx]; ylim_box = [ymin-my, ymax+my];

%% ----------------- 节点 ID → 行号 映射 -----------------
nid     = nodes(:,1);
nid_min = min(nid);  nid_max = max(nid);
id2row  = nan(nid_max - nid_min + 1, 1);
for r = 1:Nnode
    id2row( nid(r) - nid_min + 1 ) = r;
end
id2row_fun = @(id) assert_return(id2row(id - nid_min + 1), ...
    sprintf('节点ID %d 在 nodes 表中不存在！', id));

%% ----------------- 开关与样式 -----------------
showNodeLabel   = true;
showFaceLabel   = true;
showCellLabel   = true;
showSegLabel    = true;
showFFPts       = true;

psz_node  = 18; psz_cell = 36; psz_ffpt = 50;
f_line_w  = 0.8; seg_line_w = 1.8;
col_ffpt  = [0.0,0.6,0.0];
dx_node   = .01; dy_node = .01;
dx_cell   = .01; dy_cell = .01;

%% ----------------- 绘图 -----------------
figure('Color','w','Units','normalized','Position',[0.08,0.08,0.82,0.82]);
hold on; axis equal; box on;
title('EDFM mesh & fractures (no ghost)','FontWeight','bold');

% 1) Faces（用 ID → 行号 映射取端点）
for i = 1:Nface
    n1_id = faces(i,2);  n2_id = faces(i,3);
    r1 = id2row_fun(n1_id);  r2 = id2row_fun(n2_id);
    plot([nodes(r1,2), nodes(r2,2)], [nodes(r1,3), nodes(r2,3)], ...
         'k-', 'LineWidth', f_line_w);
    if showFaceLabel && size(faces,2) >= 5
        text(faces(i,4), faces(i,5)+0.02, "F"+string(faces(i,1)), ...
             'FontSize',10, 'Color','r');
    end
end

% 2) Nodes
scatter(nodes(:,2), nodes(:,3), psz_node, 'k', 'filled');
if showNodeLabel
    text(nodes(:,2)+dx_node, nodes(:,3)+dy_node, string(nodes(:,1)), ...
         'FontSize',10, 'Color',[.25 .25 .25], 'FontWeight','bold');
end

% 3) Cell centres（全部真实单元）
scatter(cells(:,2), cells(:,3), psz_cell, 'b', 'filled');
if showCellLabel
    text(cells(:,2)+dx_cell, cells(:,3)+dy_cell, string(cells(:,1)), ...
         'FontSize',10, 'Color','b', 'FontWeight','bold');
end

% —— 为图例预留基准句柄（后面统一 legend 用） ——
hFaceLegend = plot(nan,nan,'k-','LineWidth',f_line_w);    % Faces
hNodes      = scatter(nan,nan,psz_node,'k','filled');     % Nodes
hCells      = scatter(nan,nan,psz_cell,'b','filled');     % Cell centres

% 4) Fracture segments：严格用 3:6 列，并按 cellID 分组上色
segLegendHandles = gobjects(0);
segLegendLabels  = strings(0);

if ~isempty(segs)
    x1  = segs(:,3);  y1  = segs(:,4);
    x2  = segs(:,5);  y2  = segs(:,6);
    mxs = segs(:,7);  mys = segs(:,8);       % 标签位置
    fid = segs(:,1);  sid = segs(:,2);
    cid = segs(:,9);                          % 宿主 cellID

    % 便捷检查
    fprintf('[seg check] x∈[%.3g, %.3g], y∈[%.3g, %.3g]\n', ...
        min([x1;x2]), max([x1;x2]), min([y1;y2]), max([y1;y2]));

    % 按 cellID 分组上色
    [uCID,~,idxCID] = unique(cid,'stable');
    K = numel(uCID);
    cmap = lines(max(K,1));

    for g = 1:K
        mask = (idxCID == g);
        X = [x1(mask) x2(mask)]';
        Y = [y1(mask) y2(mask)]';
        line(X, Y, 'Color', cmap(g,:), 'LineWidth', seg_line_w);  % 该组所有段

        % 图例句柄（虚拟）
        h = plot(nan, nan, '-', 'Color', cmap(g,:), 'LineWidth', seg_line_w);
        segLegendHandles(end+1) = h; %#ok<SAGROW>
        segLegendLabels(end+1)  = "cell " + string(uCID(g)); %#ok<SAGROW>
    end

    % 段标签（fid.seg）
    if showSegLabel
        for k = 1:numel(x1)
            text(mxs(k), mys(k), sprintf('%d.%d', fid(k), sid(k)), ...
                 'FontSize',10, 'Color',[.1 .5 .1]);
        end
    end
end

% 5) Fracture–Fracture intersection points
if showFFPts && ~isempty(ffpt)
    hFF = scatter(ffpt(:,2), ffpt(:,3), psz_ffpt, 'o', ...
                  'MarkerEdgeColor','k', 'MarkerFaceColor',col_ffpt, 'LineWidth',1.2);
    hFFLegend = plot(nan,nan,'o','MarkerEdgeColor','k','MarkerFaceColor',col_ffpt);
else
    hFF = []; hFFLegend = [];
end

% 6) 锁定视窗到域框
xlim(xlim_box); ylim(ylim_box); axis equal; box on;

% 7) 统一图例（Faces / Nodes / Cell centres / 按 cellID 的段 / Frac–Frac pts）
legendHandles = [hFaceLegend, hNodes, hCells, segLegendHandles(:)'];
legendLabels  = ["Faces","Nodes","Cell centres", segLegendLabels];
if ~isempty(hFFLegend)
    legendHandles = [legendHandles, hFFLegend];
    legendLabels  = [legendLabels, "Frac–Frac pts"];
end
legend(legendHandles, legendLabels, 'Location','bestoutside');

% 统计输出
if ~isempty(segs)
    fprintf('[viz] nodes=%d, faces=%d, cells=%d; segments=%d; unique host cells=%d\n', ...
        Nnode, Nface, Ncell, size(segs,1), numel(unique(segs(:,9))));
else
    fprintf('[viz] nodes=%d, faces=%d, cells=%d; segments=0\n', Nnode, Nface, Ncell);
end

hold off;

%% ======================== 辅助函数 ========================
function x = assert_return(val, errmsg)
    if isempty(val) || (isnumeric(val) && any(isnan(val)))
        error(errmsg);
    end
    x = val;
end
