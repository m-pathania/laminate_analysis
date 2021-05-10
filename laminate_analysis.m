%% laminate_analysis(laminate, loads, delT, delC)
% laminate = object of Laminate Class
% loads = column vector of mechanical loads = [Nx Ny Nxy Mx My Mxy].'
% delT = Temprature Change
% delC = Change in moisture
% degradation = "complete" or "partial" depending on analysis
% degradation = "complete" by default
% Calculates First Ply Failure Load and Last Ply failure load recursiverly
% Lamina is degraded by considering complete degradation or partial
% degradation.


% f_load = First Ply Failure Load
% l_load = Last Ply Failure Load
function [f_load, l_load] = laminate_analysis(laminate, loads, delT, delC, degradation)
    % default method of degradation is complete degradation
    if nargin == 4
        degradation = "complete";
    end
    
    f_load_set = false;
    while laminate.is_failed() == false
        [sr, failed_ply, failed_mode] = ply_failure_analysis(laminate, loads, delT, delC);
        for ii = failed_ply
            laminate.status(ii) = 0;
        end
        
        if f_load_set == false
            f_load = loads/sr;
            f_load_set = true;
            l_load = f_load;   % case when all plies fail together
        else
            l_load = loads/sr;
        end
        loads = l_load;
        
        % Complete Degradation
        if strcmpi(degradation, "complete")
            for ii = failed_ply
                laminate.laminas(ii).Qbar = zeros(3);
            end
        % Partial Degradation
        elseif strcmpi(degradation, "partial")
            for ii = 1:size(failed_ply,2)
                if(failed_mode(ii) == 1)
                    laminate.laminas(failed_ply(ii)).elastic = laminate.laminas(failed_ply(ii)).elastic.*[0 1 1 0];
                else
                    laminate.laminas(failed_ply(ii)).elastic = laminate.laminas(failed_ply(ii)).elastic.*[1 0 1 0];
                end
            end
        else
            error("degradation method can be either 'complete' or 'partial'")
        end
    end
end


%% ply_failure_analysis(laminate, loads, delT, delC)
% laminate = object of Laminate Class
% loads = column vector of mechanical loads = [Nx Ny Nxy Mx My Mxy]
% delT = Temprature Change
% delC = Change in moisture
% returns the strength_ratio, plyies that failed and the mode of failure
% failure_mode_ ==> 1 = longitudinal
%                   2 = Transverse
%                   3 = Shear

function [strength_ratio, failed_ply_, failure_mode_] = ply_failure_analysis(laminate, loads, delT, delC)

    % Mechanical Loads
    strain_M_xy = global_strain(laminate, loads);
    stress_M_12 = local_stress(laminate, strain_M_xy);
    
    % Strength Ratio
    strength_ratio_table = zeros(size(stress_M_12));
    for ii = 1:length(laminate.laminas)
        % Longitudinal Stress
        if(stress_M_12(1,ii) > 0)
            strength_ratio_table(1,ii) = stress_M_12(1,ii)/laminate.laminas(ii).strength(1);
        else
            strength_ratio_table(1,ii) = -stress_M_12(1,ii)/laminate.laminas(ii).strength(2);
        end
        
        % Transverse Stress
        if(stress_M_12(2,ii) > 0)
            strength_ratio_table(2,ii) = stress_M_12(2,ii)/laminate.laminas(ii).strength(3);
        else
            strength_ratio_table(2,ii) = -stress_M_12(2,ii)/laminate.laminas(ii).strength(4);
        end
        
        % Shear Stress
        strength_ratio_table(3,ii) = abs(stress_M_12(3,ii)/laminate.laminas(ii).strength(5));
    end
    
    strength_ratio = -inf;
    for ii = 1:length(laminate.laminas)
        if (laminate.status(ii) ~= 0)       % checking if laminate has not failed already
            for jj = 1:size(strength_ratio_table, 1)
                if abs(strength_ratio - strength_ratio_table(jj, ii)) < 1e-6
                    failed_ply = [failed_ply, ii];
                    failure_mode = [failure_mode, jj];
                elseif (strength_ratio < strength_ratio_table(jj, ii))
                    strength_ratio = strength_ratio_table(jj, ii);
                    failed_ply = ii;
                    failure_mode = jj;
                end
            end
        end
    end
    
    % Residual Thermal Stress
    strain_T_xy = global_strain(laminate, thermal_load(laminate, delT));
    for ii = 1:size(strain_T_xy, 2)
        strain_T_xy(:, ii) = strain_T_xy(:, ii) - delT*laminate.laminas(ii).thermal;
    end
    stress_T_12 = local_stress(laminate, strain_T_xy);
    
    % Residual hygro stress
    strain_C_xy = global_strain(laminate, hygro_load(laminate, delC));
    for ii = 1:size(strain_T_xy, 2)
        strain_C_xy(:, ii) = strain_C_xy(:, ii) - delT*laminate.laminas(ii).hygro;
    end
    stress_C_12 = local_stress(laminate, strain_C_xy);
    
    % Residual Stresses
    residual_stress = stress_T_12 + stress_C_12;
    
    % Again Calculating strength Ratio for the lamina that failed from
    % mechanical loading after include effect of residual stressed by
    % hygrothermal effects.
    strength_ratio_table = zeros(size(failed_ply));
    for ii = 1:size(failed_ply,2)
        if(failure_mode(ii) == 1)               % failed due to longitudinal stress
            if(stress_M_12(1,failed_ply(ii)) > 0)
                strength_ratio_table(ii) = stress_M_12(1,failed_ply(ii))/(laminate.laminas(ii).strength(1) - residual_stress(1, failed_ply(ii)));
            else
                strength_ratio_table(ii) = abs(stress_M_12(1,failed_ply(ii))/(laminate.laminas(ii).strength(2) - residual_stress(1, failed_ply(ii))));
            end
        elseif(failure_mode(ii) == 2)           % failed due to transverse stress
            if(stress_M_12(2,failed_ply(ii)) > 0)
                strength_ratio_table(ii) = stress_M_12(2,failed_ply(ii))/(laminate.laminas(ii).strength(3) - residual_stress(2, failed_ply(ii)));
            else
                strength_ratio_table(ii) = abs(stress_M_12(2,failed_ply(ii))/(laminate.laminas(ii).strength(4) - residual_stress(2, failed_ply(ii))));
            end
        else                                    % failed due to shear stress
            strength_ratio_table(ii) = abs(stress_M_12(3,failed_ply(ii))/(laminate.laminas(ii).strength(5) - residual_stress(3, failed_ply(ii))));
        end
    end
    
    strength_ratio = -inf;
    for ii = 1:length(strength_ratio_table)
        if abs(strength_ratio - strength_ratio_table(ii)) < 1e-6
            failed_ply_ = [failed_ply_, failed_ply(ii)];
            failure_mode_ = [failure_mode_, failure_mode(ii)];
        elseif (strength_ratio < strength_ratio_table(ii))
            strength_ratio = strength_ratio_table(ii);
            failed_ply_ = failed_ply(ii);
            failure_mode_ = failure_mode(ii);
        end
    end
end


%% global_strain(laminate, loads)
% laminate = object of Laminate Class
% loads = column vector of (mechanical, thermal or hygro) 
% loads = [Nx Ny Nxy Mx My Mxy].'
% return a metrix containing global strains as rows
% columns represent the lamina sequence

function ex = global_strain(laminate, loads)
    strains_mid = laminate.ABD()\loads;
    ex = zeros(3, size(laminate.laminas,2));
    for ii = 1:size(laminate.laminas,2)
        ex(:,ii) = strains_mid(1:3) + 0.5*(laminate.z(ii) + laminate.z(ii + 1))*strains_mid(4:6);
    end
end


%% global_stress(laminate, ex)
% laminate = object of Laminate Class
% ex = global strains for each lamina
% return a metrix containing global stresses as rows
% columns represent the lamina sequence

function sx = global_stress(laminate, ex)
    sx = zeros(3, size(laminate.laminas,2));
    for ii = 1:size(laminate.laminas,2)
        sx(:,ii) = laminate.laminas(ii).Qbar*ex(:,ii);
    end
end


%% local_strain(laminate, ex)
% laminate = object of Laminate Class
% ex = global strains for each lamina
% return a metrix containing local strains as rows
% columns represent the lamina sequence

function e12 = local_strain(laminate, ex)
    e12 = zeros(size(ex));
    R = diag([1; 1; 2]);
    for ii = 1:size(laminate.laminas,2)
        e12(:,ii) = R*laminate.laminas(ii).T()*laminate.laminas(ii).Qbar*(R\ex(:,ii));
    end
end

%% local_stress(laminate, ex)
% laminate = object of Laminate Class
% ex = global strains for each lamina
% return a metrix containing local stresses as rows
% columns represent the lamina sequence

function s1 = local_stress(laminate, ex)
    sxy = global_stress(laminate, ex);
    s1 = zeros(size(sxy));
    for ii = 1:length(laminate.laminas)
        s1(:,ii) = laminate.laminas(ii).T()*sxy(:, ii);
    end
end


%% thermal_load(laminate, delT)
% laminate = object of Laminate Class
% delT = Temprature Change
% Returns Thermal equivaled load as a column vector
% out = [Nx Ny Nxy Mx My Mxy].'

function out = thermal_load(laminate, delT)
    Nt = zeros(3,1);
    Mt = zeros(3,1);
    for ii = 1:length(laminate.laminas)
        Nt = Nt + laminate.laminas(ii).Qbar*laminate.laminas(ii).thermal*(laminate.z(ii + 1) - laminate.z(ii));
        Mt = Mt + laminate.laminas(ii).Qbar*laminate.laminas(ii).thermal*(laminate.z(ii + 1)^2 - laminate.z(ii)^2);
    end
    out = [delT*Nt; 0.5*delT*Mt];
end


%% hygro_load(laminate, delC)
% laminate = object of Laminate Class
% delT = Moisture Change
% Returns equivaled load due to moisture change as a column vector
% out = [Nx Ny Nxy Mx My Mxy].'

function out = hygro_load(laminate, delC)
    Nc = zeros(3,1);
    Mc = zeros(3,1);
    for ii = 1:length(laminate.laminas)
        Nc = Nc + laminate.laminas(ii).Qbar*laminate.laminas(ii).hygro*(laminate.z(ii + 1) - laminate.z(ii));
        Mc = Mc + laminate.laminas(ii).Qbar*laminate.laminas(ii).hygro*(laminate.z(ii + 1)^2 - laminate.z(ii)^2);
    end
    out = [delC*Nc; 0.5*delC*Mc];
end


%% First and Last Ply Failure Loads

% Coding Assignment
% ME607: Introduction to Composite Materials

% Name = Mayank Pathania
% Roll No. = 204103314
% Specialization = Machine Design
% Indian Institute of Technology, Guwahati
