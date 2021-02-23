MeasuredGeneratorSpeedVec = [0:0.1:1.5]*97;
GeneratorTorqueDemandVec = nan(size(MeasuredGeneratorSpeedVec));


for idx = 1: length(MeasuredGeneratorSpeedVec)
    
    MeasuredGeneratorSpeed = MeasuredGeneratorSpeedVec(idx);
    GeneratorTorqueDemand = max(min(Control_Torque_OptGain*MeasuredGeneratorSpeed^2, Torque_Demanded), T_MinTorque);
    
    % Calculate torques at beginning and end of CP_max operation
    TorqueOmegaB = Control_Torque_OptGain*(omega_ABCD(2)*pi/30)^2;
    TorqueOmegaC = Control_Torque_OptGain*(omega_ABCD(3)*pi/30)^2;
    
    % Apply maximum demanded torque when pitch is 1 degree above fine-pitch and
    % the wind turbine is not in start-up mode.
    % if (Event == 5 && Clock >= EventTime) && (MeasuredPitch > (FineMaxPitch(2) - 1.0)*pi/180) %Emergency shutdown
    %     GeneratorTorqueDemand = 0.0;
    % else   
    if MeasuredGeneratorSpeed < omega_ABCD(1)*pi/30
        GeneratorTorqueDemand = T_MinTorque;
    elseif MeasuredGeneratorSpeed < omega_ABCD(2)*pi/30
        GeneratorTorqueDemand = (TorqueOmegaB - T_MinTorque)/((omega_ABCD(2) - omega_ABCD(1))*pi/30)*(MeasuredGeneratorSpeed - omega_ABCD(1)*pi/30) + T_MinTorque;
    elseif MeasuredGeneratorSpeed >= omega_ABCD(3)*pi/30
        GeneratorTorqueDemand = (Torque_Demanded - TorqueOmegaC)/((omega_ABCD(4) - omega_ABCD(3))*pi/30)*(MeasuredGeneratorSpeed - omega_ABCD(3)*pi/30) + TorqueOmegaC;
    end
    GeneratorTorqueDemandVec(idx) = GeneratorTorqueDemand;
end

plot(MeasuredGeneratorSpeedVec, GeneratorTorqueDemandVec)