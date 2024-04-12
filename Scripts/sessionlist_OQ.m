%%
% open loop
PreprocessSmellocatorData('Q4_20221112_r0.mat',1);
PreprocessSmellocatorData('Q5_20221122_r0.mat',1);
PreprocessSmellocatorData('Q8_20221209_r0.mat',1);
PreprocessSmellocatorData('Q9_20221119_r0.mat',1);

% halts
PreprocessSmellocatorData('Q4_20221109_r0.mat',1);
PreprocessSmellocatorData('Q8_20221204_r0.mat',1);
PreprocessSmellocatorData('Q8_20221207_r0.mat',1);

% open loop
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O3/O3_20211005_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O7/O7_20220630_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O8/O8_20220702_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O9/O9_20220630_r0.mat',1);

% halts
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O3/O3_20210927_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O3/O3_20210929_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O7/O7_20220702_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O8/O8_20220704_r0.mat',1);
PreprocessSmellocatorData('/mnt/grid-hs/pgupta/Behavior/O9/O9_20220702_r1.mat',1);




%% bunch of sessions
Session{1,1} = '/mnt/grid-hs/pgupta/Behavior/Q3/Q3_20221019_r0.mat'; % halt
Session{2,1} = '/mnt/grid-hs/pgupta/Behavior/Q4/Q4_20221109_r0.mat'; % halt
Session{3,1} = '/mnt/grid-hs/pgupta/Behavior/Q5/Q5_20221115_r0.mat'; % halt
Session{4,1} = '/mnt/grid-hs/pgupta/Behavior/Q8/Q8_20221204_r0.mat'; % halt
Session{5,1} = '/mnt/grid-hs/pgupta/Behavior/Q9/Q9_20221112_r0.mat'; % halt

Session{6,1} = '/mnt/grid-hs/pgupta/Behavior/S1/S1_20230327_r0.mat'; % halt
Session{7,1} = '/mnt/grid-hs/pgupta/Behavior/S3/S3_20230321_r0.mat'; % OL replay
Session{8,1} = '/mnt/grid-hs/pgupta/Behavior/S6/S6_20230710_r0.mat'; % halt
Session{9,1} = '/mnt/grid-hs/pgupta/Behavior/S11/S11_20230801_r0.mat'; % halt
Session{10,1} = '/mnt/grid-hs/pgupta/Behavior/S12/S12_20230731_r0.mat'; % halt

for i = 1:10
    BasicSniffAnalysis(Session{i,1},1);
end