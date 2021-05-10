%% Lamina Class
% geometric = [thickness, angle]
% elastic = [E1, E2, v12, G12]
% strength = [s1tu, s1cu, s2tu, s2cu, t12u]
% hygro = [betax, betay, betaxy].'
% thermal = [alphax, alphay, alphaxy].'
% Qbar = Reduced Transformed Stiffness Matrix


classdef Lamina
    properties
        geometric
        elastic
        strength
        hygro
        thermal
        Qbar
    end
    methods
        function obj = Lamina(geo, elas, stren, hygro)
            obj.geometric = [geo(1), geo(2)];
            obj.elastic = elas;
            obj.strength = stren;

            R = diag([1 1 2]);  % Reuters Matrix
            T = obj.T();        % Transformation Matrix
            
            obj.Qbar = (((T\obj.Q())*R)*T)/R;
            obj.thermal = (R*(T\[hygro(1); hygro(2); 0]));
            obj.hygro = (R*(T\[hygro(3); hygro(4); 0]));
        end
        function out = T(obj)
            % Returns Transformation Matrix
            a = obj.geometric(2);
            out = [cos(a)^2 sin(a)^2 2*sin(a)*cos(a);
                   sin(a)^2 cos(a)^2 -2*sin(a)*cos(a);
                   -sin(a)*cos(a) sin(a)*cos(a) cos(a)^2 - sin(a)^2];
        end
        function out = S(obj)
            % returns compliance matrix
            E1 = obj.elastic(1);
            E2 = obj.elastic(2);
            v12 = obj.elastic(3);
            G12 = obj.elastic(4);
            out = [1/E1 -v12/E1 0;
                -v12/E1 1/E2 0;
                0 0 1/G12];
        end
        function out = Q(obj)
            % returns stiffness matrix
            out = inv(obj.S());
        end
    end
end


%% First and Last Ply Failure Loads

% Coding Assignment
% ME607: Introduction to Composite Materials

% Name = Mayank Pathania
% Roll No. = 204103314
% Specialization = Machine Design
% Indian Institute of Technology, Guwahati
