%% File lists for sharing data with Ryan and Ben

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/Q3/*r*.mat');
Q3_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:26
    PreprocessSmellocatorData_basic_v2(Q3_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/Q4/*r*.mat');
Q4_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:42
    PreprocessSmellocatorData_basic_v2(Q4_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/Q5/*r*.mat');
Q5_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:32
    PreprocessSmellocatorData_basic_v2(Q5_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/Q8/*r*.mat');
Q8_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:39
    PreprocessSmellocatorData_basic_v2(Q8_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/Q9/*r*.mat');
Q9_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:25
    PreprocessSmellocatorData_basic_v2(Q9_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/O3/*r*.mat');
O3_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:23
    PreprocessSmellocatorData_basic_v2(O3_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/O8/*r*.mat');
O8_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:30
    PreprocessSmellocatorData_basic_v2(O8_Sessions{i},1);
end

AllFiles = dir('/mnt/grid-hs/pgupta/Behavior/O9/*r*.mat');
O9_Sessions = arrayfun(@(x) x.name, AllFiles, 'UniformOutput', false);
for i = 1:28
    if i~=4 % noisy lever signal
        PreprocessSmellocatorData_basic_v2(O9_Sessions{i},1);
    end
end