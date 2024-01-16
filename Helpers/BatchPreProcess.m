function [] = BatchPreProcess()
    % O3
    PreprocessSmellocatorData('O3_20210927_r0.mat',1);
    PreprocessSmellocatorData('O3_20210929_r0.mat',1);
    PreprocessSmellocatorData('O3_20211005_r0.mat',1);
    
    % O8
    PreprocessSmellocatorData('O8_20220702_r0.mat',1);
    PreprocessSmellocatorData('O8_20220704_r0.mat',1);
    PreprocessSmellocatorData('O8_20220707_r0.mat',1);
    
    % O9
    PreprocessSmellocatorData('O9_20220630_r0.mat',1);
    PreprocessSmellocatorData('O9_20220702_r1.mat',1);
    
    % S1
    PreprocessSmellocatorData('S1_20230312_r0.mat',1);
    PreprocessSmellocatorData('S1_20230314_r0.mat',1);
    PreprocessSmellocatorData('S1_20230327_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S1_20230329_r0.mat';'S1_20230329_r2.mat'],1);
    PreprocessSmellocatorData('S1_20230403_r0.mat',1);
    PreprocessSmellocatorData('S1_20230405_r0.mat',1);
    PreprocessSmellocatorData('S1_20230407_r0.mat',1);
    PreprocessSmellocatorData('S1_20230413_r0.mat',1);
    
    % S3
    PreprocessSmellocatorData('S3_20230316_r0.mat',1);
    PreprocessSmellocatorData('S3_20230321_r0.mat',1);
    PreprocessSmellocatorData('S3_20230327_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S3_20230329_r0.mat';'S3_20230329_r1.mat'],1);
    PreprocessSmellocatorData('S3_20230403_r0.mat',1);
    PreprocessSmellocatorData('S3_20230405_r0.mat',1);
    PreprocessSmellocatorData('S3_20230407_r0.mat',1);
    
    % S6
    PreprocessSmellocatorData('S6_20230622_r0.mat',1);
    PreprocessSmellocatorData('S6_20230710_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S6_20230713_r0.mat';'S6_20230713_r1.mat'],1);
    PreprocessSmellocatorData('S6_20230718_r0.mat',1);
    PreprocessSmellocatorData('S6_20230721_r0.mat',1);
    PreprocessSmellocatorData('S6_20230724_r0.mat',1);
    PreprocessSmellocatorData('S6_20230727_r0.mat',1);
    
    % S7
    PreprocessSmellocatorData('S7_20230530_r0.mat',1);
    PreprocessSmellocatorData('S7_20230608_r0.mat',1);
    PreprocessSmellocatorData('S7_20230616_r0.mat',1);
    PreprocessSmellocatorData('S7_20230622_r0.mat',1);
    PreprocessSmellocatorData('S7_20230627_r0.mat',1);
    PreprocessSmellocatorData('S7_20230629_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S7_20230705_r0.mat';'S7_20230705_r1.mat'],1);
    PreprocessSmellocatorData('S7_20230707_r0.mat',1);
    
    % S11
    PreprocessSmellocatorData('S11_20230801_r0.mat',1);
    PreprocessSmellocatorData('S11_20230805_r0.mat',1);
    PreprocessSmellocatorData('S11_20230807_r0.mat',1);
    PreprocessSmellocatorData('S11_20230812_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S11_20230813_r0.mat';'S11_20230813_r1.mat'],1);
    
    % S12
    PreprocessSmellocatorData('S12_20230727_r0.mat',1);
    PreprocessSmellocatorData('S12_20230731_r0.mat',1);
    PreprocessSmellocatorData('S12_20230804_r0.mat',1);
    PreprocessSmellocatorData_Multi(['S12_20230809_r0.mat';'S12_20230809_r1.mat'],1);
    
end