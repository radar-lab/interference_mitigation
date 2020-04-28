function [SIRValue] = SIRCacl(data,TargetIdx,RefCells,GuardCells)
    % Calculate SIR
    NoisePwr    = sum(abs(data((TargetIdx-RefCells/2-GuardCells/2) : (TargetIdx-GuardCells/2-1))).^2);
    NoisePwr    = NoisePwr + sum(abs(data((TargetIdx+GuardCells/2+1) : (TargetIdx+RefCells/2+GuardCells/2))).^2);
    NoisePwr    = NoisePwr/RefCells;
    TargetPwr   = abs(data(TargetIdx))^2;
    SIRValue    = 10*log10( TargetPwr / NoisePwr);
end

