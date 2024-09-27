
SessionName = 'S12_20230731_r0';
SessionName = 'O3_20211005_r0';

%% load a particular fit config
Paths = WhichComputer();
for whichset = 1:6
    FitPath = fullfile(Paths.Wolf.processed,'sniffPSTHPredictions',SessionName(1:regexp(SessionName,'_','once')-1),SessionName);

    switch whichset
        case 1
            FitType = '_lsqcurvefit_1_20ms_norectify.mat';
        case 2
            FitType = '_lsqcurvefit_2_20ms_norectify.mat';
        case 3
            FitType = '_lsqcurvefit_3_20ms_norectify.mat';
        case 4
            FitType = '_lsqcurvefit_1_20ms_rectify.mat';
        case 5
            FitType = '_lsqcurvefit_2_20ms_rectify.mat';
        case 6
            FitType = '_lsqcurvefit_3_20ms_rectify.mat';
    end

    FitPath = fullfile(FitPath,[SessionName,FitType]);
    load(FitPath, 'InputVector','PSTHs','fittedkernel','PSTHbinsize');

    nUnits = size(PSTHs,1);

    %% get predictions
    for i = 1:nUnits
        % get kernels
        [baseline,kernels,locationcoef] = ParseSniffKernels(fittedkernel{i},'independentcoeffs', true);
        [zdata] = SniffKernels2continuousPSTH(baseline,kernels,locationcoef,InputVector(1:7,:));
        zdata(zdata<0) = 0;
        PredictedPSTH{whichset}(i,:) = zdata;
    end

end

%% get the timestamps from the original processed file
WhereSession = fullfile(Paths.Wolf.processed,'forWDW',[SessionName, '_processed.mat']);
load(WhereSession, 'TracesOut');
downsample  = PSTHbinsize/(1000/500);
TS_temp     = TracesOut.Timestamps{1};
TS     = TS_temp(1):PSTHbinsize/1000:TS_temp(end);

%% get fits from wulf
wulfPath = fullfile(Paths.Wolf.processed,'sniffPSTHPredictions',SessionName(1:regexp(SessionName,'_','once')-1),SessionName,[SessionName,'_wdw.mat']);
WulfPSTHs = LoadWulfGLMOutput(wulfPath,nUnits,TS);

%% Looking at a particular unit
window = [-0.1 1]; % in seconds
window = window/(PSTHbinsize/1000); % in downsampled bins
CLsniffs = [];
for whichUnit = 1:nUnits
    %whichUnit = 58;
    if isempty(CLsniffs) % only do this once
        [CLsniffs] = PlotSniffRaster_wdw(SessionName, whichUnit);
        CLsniffs(1:5,:) = [];
        CLsniffs(end-4:end,:) = [];
    end

    %% make the psths
    for sniffType = 1:5
        whichSniffs = find(CLsniffs(:,7)==sniffType);
        summedTrace = 0; predTrace = zeros(7,diff(window)+1); predictedTraces = [];
        for n = 1:numel(whichSniffs)
            sniffStart = CLsniffs(whichSniffs(n),8);
            [~,idx] = min(abs(TS-sniffStart));
            idx = idx + window;
            summedTrace = summedTrace + PSTHs(whichUnit,idx(1):idx(2));
            for whichset = 1:6
                predictedTraces(whichset,:) = PredictedPSTH{whichset}(whichUnit,idx(1):idx(2));
            end
            predictedTraces(7,:) = WulfPSTHs(whichUnit,idx(1):idx(2));
            predTrace = predTrace + predictedTraces;
        end
        summedTrace = summedTrace/n;
        predTrace   = predTrace/n;
        SniffPSTH{whichUnit}{sniffType} = summedTrace;
        ExpectedPSTH{whichUnit}{sniffType} = predTrace;

        for whichset = 1:7
            mycorr = corrcoef(summedTrace,predTrace(whichset,:));
            PSTHcorrs{whichUnit}{sniffType}(whichset) = mycorr(1,2);
        end
    end
    
    
    %% 

    %%
%     figure;
%     for sniffType = 1:5
%         subplot(1,5,sniffType);
%         plot(sgolayfilt(SniffPSTH{sniffType},1,3),'k');
%         hold on
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(1,:),1,3),'r');
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(2,:),1,3),'g');
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(3,:),1,3),'b');
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(4,:),1,3),':r');
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(5,:),1,3),':g');
%         plot(sgolayfilt(ExpectedPSTH{sniffType}(6,:),1,3),':b');
%         %plot(sgolayfilt(ExpectedPSTH{sniffType}',1,3));
%     end
end
%% compare sniffPSTHs for a particular unit
whichUnit = 1;
figure; 
for sniffType = 1:5
    subplot(1,5,sniffType); 
    hold on;
    plot(sgolayfilt(SniffPSTH{whichUnit}{sniffType},1,3),'k');
    plot(ExpectedPSTH{whichUnit}{sniffType}(7,:),'r'); 
    plot(ExpectedPSTH{whichUnit}{sniffType}(1,:),'b'); 
    plot(ExpectedPSTH{whichUnit}{sniffType}(4,:),'g'); 
end


%% comparing R2 for psth
figure;
for sniffType = 1:5
    M(:,:,sniffType) = cell2mat(cellfun(@(x) x{sniffType}, PSTHcorrs, 'UniformOutput', false)'); 
    subplot(1,5,sniffType);
    errorbar(1:7,mean(M(:,:,sniffType)),std(M(:,:,sniffType)));
end


