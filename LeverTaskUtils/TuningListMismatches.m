function [OdorMismatches] = TuningListMismatches(TrialSequence,EphysTuningTrials)
lastrow = min(size(TrialSequence,1),size(EphysTuningTrials,1));
OdorMismatches = intersect(find(TrialSequence(2:lastrow,1)-EphysTuningTrials(2:lastrow,1)),...
                find(~isnan(TrialSequence(2:lastrow,1)-EphysTuningTrials(2:lastrow,1))));