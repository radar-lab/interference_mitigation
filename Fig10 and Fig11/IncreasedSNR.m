function [increasedSNR] = IncreasedSNR(priCH,e,TargetIdx,RefCells,GuardCells)
    % Calculate SNR before ANC
    NoisePwr = sum(abs(priCH((TargetIdx-RefCells/2-GuardCells/2) : (TargetIdx-GuardCells/2-1))).^2);
    NoisePwr = NoisePwr + sum(abs(priCH((TargetIdx+GuardCells/2+1) : (TargetIdx+RefCells/2+GuardCells/2))).^2);
    NoisePwr = NoisePwr/RefCells;
    TargetPwr           = abs(priCH(TargetIdx))^2;
    SIRBeforeANC        = 10*log10( TargetPwr / NoisePwr);
    % Calculate SNR after ANC
    NoisePwr = sum(abs(e((TargetIdx - RefCells/2 - GuardCells/2) : (TargetIdx - GuardCells/2 - 1))).^2);
    NoisePwr = NoisePwr + sum(abs(e((TargetIdx+GuardCells/2+1) : (TargetIdx+RefCells/2+GuardCells/2))).^2);
    NoisePwr = NoisePwr/RefCells;
    TargetPwr           = abs(e(TargetIdx))^2;
    SIRAfterANC        = 10*log10( TargetPwr / NoisePwr);
    % Calculate increased SNR
    increasedSNR = SIRAfterANC - SIRBeforeANC;
end

