load('/mnt/data/Sorted/T2/_2025-05-21_09-18-56_2025-05-21_11-11-12_2025-05-21_11-23-51/OdorMaps/SniffOdorSummaryFirstPulseParams.mat')
figure;

AirPSTH = []; 
OdorPSTH = [];
for x = 1:size(AllunitsPSTH,2)
    AirPSTH = vertcat(AirPSTH,AllunitsPSTH{x}(1:4,:));
    meanAir = mean(AllunitsPSTH{x}(1:4,:),1);
    OdorPSTH = vertcat(OdorPSTH,(AllunitsPSTH{x}(5:19,:)-meanAir));
end

subplot(1,4,1);
imagesc(AirPSTH./max(AirPSTH,[],2),[-1 1]);

subplot(1,4,3);
imagesc(OdorPSTH./max(abs(OdorPSTH),[],2),[-1 1]);

load('/mnt/data/Sorted/T3/_2025-05-16_13-48-44_2025-05-16_15-40-38_2025-05-16_15-49-31/OdorMaps/SniffOdorSummaryFirstPulseParams.mat')

AirPSTH = []; 
OdorPSTH = [];
for x = 1:size(AllunitsPSTH,2)
    AirPSTH = vertcat(AirPSTH,AllunitsPSTH{x}(1:4,:));
    meanAir = mean(AllunitsPSTH{x}(1:4,:),1);
    OdorPSTH = vertcat(OdorPSTH,AllunitsPSTH{x}(5:19,:)-meanAir);
end

subplot(1,4,2);
imagesc(AirPSTH./max(AirPSTH,[],2),[-1 1]);

subplot(1,4,4);
imagesc(OdorPSTH./max(abs(OdorPSTH),[],2),[-1 1]);

colormap(brewermap([100],'RdBu'))
