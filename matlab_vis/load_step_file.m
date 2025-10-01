function [fname, data] = load_step_file(folder, step)
% 读一个时间步的 txt/csv 文件（你的 out_txt_CO2 里这种行格式）
% 返回：
%  data.cell_id, data.cx, data.cy, data.cz, data.vol, data.p, data.T, data.step, data.time

% 匹配 step_******.txt / pT_step_******.csv
cand = dir(fullfile(folder, sprintf('*%06d*.txt', step)));
if isempty(cand)
    cand = dir(fullfile(folder, sprintf('*%06d*.csv', step)));
end
if isempty(cand)
    error('Cannot find step file for step=%d in %s', step, folder);
end
fname = fullfile(folder, cand(1).name);

% 读取（忽略以 # 开头的注释）
fid = fopen(fname, 'r');
if fid < 0, error('Cannot open %s', fname); end
C = textscan(fid, '%f %f %f %f %f %f %f', 'CommentStyle', '#', 'MultipleDelimsAsOne', true);
fclose(fid);

data.cell_id = C{1};
data.cx      = C{2};
data.cy      = C{3};
data.cz      = C{4};
data.vol     = C{5};
data.p       = C{6};
data.T       = C{7};

data.step = step;
data.time = parse_time_from_header(fname);
end
