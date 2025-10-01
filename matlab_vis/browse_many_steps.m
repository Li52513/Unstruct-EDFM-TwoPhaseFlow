function browse_many_steps(folder, stepStart, stepEnd, stepStride)
% 逐步浏览（每隔 stepStride 步画一次）
if nargin < 4, stepStride = 10; end
for s = stepStart:stepStride:stepEnd
    try
        plot_one_step(folder, s);
        drawnow;
    catch ME
        warning('Failed step %d: %s', s, ME.message);
    end
end
end