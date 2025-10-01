function t = parse_time_from_header(fname)
% 从首行形如 "# step 1000  time 9.99999e+00" 的注释中解析 time；失败返回 NaN
fid = fopen(fname, 'r'); t = NaN;
if fid < 0, return; end
line = fgetl(fid); fclose(fid);
if ~ischar(line), return; end
tok = regexp(line, 'time\s+([0-9eE\+\-\.]+)', 'tokens', 'once');
if ~isempty(tok), t = str2double(tok{1}); end
end
