%% First and Last Ply Failure Loads

% Coding Assignment
% ME607: Introduction to Composite Materials

% Name = Mayank Pathania
% Roll No. = 204103314
% Specialization = Machine Design
% Indian Institute of Technology, Guwahati

clear; clc;
%% Reading Input Data
fid = fopen("output.dat","w");

disp("Input Data (respective units as in csv files) ===>")
properties_geometric = readData("properties_geometry.csv")
properties_elastic = readData("properties_elastic.csv")
properties_strength = readData("properties_strength.csv")
properties_hygrothermal = readData("properties_hygrothermal.csv")
external_load = readData("external_loads.csv").'
%% Converting Data to SI Units

properties_elastic = properties_elastic.*[1e9 1e9 1 1e9];
properties_geometric = properties_geometric.*[1e-3 (pi/180)];
properties_strength = properties_strength.*[1 1 1 1 1]*1e6;
properties_hygrothermal = properties_hygrothermal.*[1e-6 1e-6 1 1];

delT = external_load(7);
delC = external_load(8);

external_load = external_load(1:6);
external_load = external_load.*([1 1 1 1 1 1].'*1e3);
%% Creating Laminate Object

fprintf("\n\n-------------------------------------------------------------------------\n")
laminate = Laminate(properties_geometric, properties_elastic, properties_strength, properties_hygrothermal);

fprintf("\nNumber of lamina =\t%d\n\n", size(properties_geometric, 1))
disp("Surfaces of laminate (m)")
disp(laminate.z)
fprintf("\n\nABD Matrix of laminate ===>\n")
disp(laminate.ABD())

fprintf(fid, "\nNumber of lamina =\t%d\n\n", size(properties_geometric, 1));
fprintf(fid, "\nSurfaces of laminate (m)\n");
fprintf(fid,"\t%f\n",laminate.z);

fprintf("-------------------------------------------------------------------------\n")
%% Mid Surface Analysis

fprintf("\n\n-------------------------------------------------------------------------\n")
strains_mid = laminate.ABD()\external_load;
disp(" ")
disp("Mid Surface Strains ===>")
disp(strains_mid(1:3))
disp("Mid Surface Curvature ===>")
disp(strains_mid(4:6))

fprintf(fid, "\nMid Surface Strains ===>\n");
fprintf(fid, "\t%f\n",strains_mid(1:3));
fprintf(fid, "\nMid Surface Curvature ===>\n");
fprintf(fid, "\t%f\n", strains_mid(4:6));
fprintf("-------------------------------------------------------------------------\n")

%% Analysis of each lamina

% Global Strains
fprintf("\n\n-------------------------------------------------------------------------\n")
fprintf(fid, "\n\n-------------------------------------------------------------------------\n");
global_strains = global_strain(laminate, external_load);
for ii = 1:size(global_strains,2)
    fprintf("\nglobal strains in lamina %d ==>\n",ii)
    fprintf("\t\t%f\n",global_strains(:,ii))
    fprintf(fid, "\nglobal strains in lamina %d ==>\n",ii);
    fprintf(fid, "\t\t%f\n",global_strains(:,ii));
end
fprintf("-------------------------------------------------------------------------\n")
fprintf(fid, "-------------------------------------------------------------------------\n");

% Global Stresses
fprintf("\n\n-------------------------------------------------------------------------\n")
fprintf(fid, "\n\n-------------------------------------------------------------------------\n");
global_stresses = global_stress(laminate, global_strains);
for ii = 1:size(global_strains,2)
    fprintf("\nglobal stresses in lamina %d ==>\n",ii)
    fprintf("\t\t%f MPa\n",global_stresses(:,ii))
    fprintf(fid, "\nglobal stresses in lamina %d ==>\n",ii);
    fprintf(fid, "\t\t%f MPa\n",global_stresses(:,ii));
end
fprintf("-------------------------------------------------------------------------\n")
fprintf(fid, "-------------------------------------------------------------------------\n");

% Local Strains
fprintf("\n\n-------------------------------------------------------------------------\n")
fprintf(fid, "\n\n-------------------------------------------------------------------------\n");
local_strains = local_strain(laminate, global_strains);
for ii = 1:size(global_strains,2)
    fprintf("\nlocal strains in lamina %d ==>\n",ii)
    fprintf("\t\t%f \n",local_strains(:,ii))
    fprintf(fid, "\nlocal strains in lamina %d ==>\n",ii);
    fprintf(fid, "\t\t%f \n",local_strains(:,ii));
end
fprintf("-------------------------------------------------------------------------\n")
fprintf(fid, "-------------------------------------------------------------------------\n");

% Local Stresses
fprintf("\n\n-------------------------------------------------------------------------\n")
fprintf(fid, "\n\n-------------------------------------------------------------------------\n");
local_stresses = local_stress(laminate, global_strains);
for ii = 1:size(global_strains,2)
    fprintf("\nlocal stresses in lamina %d ==>\n",ii)
    fprintf("\t\t%f MPa\n",local_stresses(:,ii))
    fprintf(fid, "\nlocal stresses in lamina %d ==>\n",ii);
    fprintf(fid, "\t\t%f MPa\n",local_stresses(:,ii));
end
fprintf("-------------------------------------------------------------------------\n")
fprintf(fid, "-------------------------------------------------------------------------\n");


%% Laminate Analysis

[first_ply_failure_load, last_ply_failure_load_complete] = laminate_analysis(laminate, external_load, delT, delC, "complete");
[~, last_ply_failure_load_partial] = laminate_analysis(laminate, external_load, delT, delC, "partial");

fprintf("\n\nFirst Ply Failure Load ===>\n")
fprintf("Nx =\t%f N/m\n",first_ply_failure_load(1))
fprintf("Ny =\t%f N/m\n",first_ply_failure_load(2))
fprintf("Nxy =\t%f N/m\n",first_ply_failure_load(3))
fprintf("Mx =\t%f N\n",first_ply_failure_load(4))
fprintf("My =\t%f N\n",first_ply_failure_load(5))
fprintf("Mxy =\t%f N\n",first_ply_failure_load(6))

fprintf("\n\nLast Ply Failure Load (Complete Degradation) ===>\n")
fprintf("Nx =\t%f N/m\n",last_ply_failure_load_complete(1))
fprintf("Ny =\t%f N/m\n",last_ply_failure_load_complete(2))
fprintf("Nxy =\t%f N/m\n",last_ply_failure_load_complete(3))
fprintf("Mx =\t%f N\n",last_ply_failure_load_complete(4))
fprintf("My =\t%f N\n",last_ply_failure_load_complete(5))
fprintf("Mxy =\t%f N\n",last_ply_failure_load_complete(6))

fprintf("\n\nLast Ply Failure Load (Partial Degradation) ===>\n")
fprintf("Nx =\t%f N/m\n",last_ply_failure_load_partial(1))
fprintf("Ny =\t%f N/m\n",last_ply_failure_load_partial(2))
fprintf("Nxy =\t%f N/m\n",last_ply_failure_load_partial(3))
fprintf("Mx =\t%f N\n",last_ply_failure_load_partial(4))
fprintf("My =\t%f N\n",last_ply_failure_load_partial(5))
fprintf("Mxy =\t%f N\n",last_ply_failure_load_partial(6))



fprintf(fid, "\n\nFirst Ply Failure Load ===>\n");
fprintf(fid, "Nx =\t%f N/m\n",first_ply_failure_load(1));
fprintf(fid, "Ny =\t%f N/m\n",first_ply_failure_load(2));
fprintf(fid, "Nxy =\t%f N/m\n",first_ply_failure_load(3));
fprintf(fid, "Mx =\t%f N\n",first_ply_failure_load(4));
fprintf(fid, "My =\t%f N\n",first_ply_failure_load(5));
fprintf(fid, "Mxy =\t%f N\n",first_ply_failure_load(6));

fprintf(fid, "\n\nLast Ply Failure Load (Complete Degradation) ===>\n");
fprintf(fid, "Nx =\t%f N/m\n",last_ply_failure_load_complete(1));
fprintf(fid, "Ny =\t%f N/m\n",last_ply_failure_load_complete(2));
fprintf(fid, "Nxy =\t%f N/m\n",last_ply_failure_load_complete(3));
fprintf(fid, "Mx =\t%f N\n",last_ply_failure_load_complete(4));
fprintf(fid, "My =\t%f N\n",last_ply_failure_load_complete(5));
fprintf(fid, "Mxy =\t%f N\n",last_ply_failure_load_complete(6));

fprintf(fid, "\n\nLast Ply Failure Load (Partial Degradation) ===>\n");
fprintf(fid, "Nx =\t%f N/m\n",last_ply_failure_load_partial(1));
fprintf(fid, "Ny =\t%f N/m\n",last_ply_failure_load_partial(2));
fprintf(fid, "Nxy =\t%f N/m\n",last_ply_failure_load_partial(3));
fprintf(fid, "Mx =\t%f N\n",last_ply_failure_load_partial(4));
fprintf(fid, "My =\t%f N\n",last_ply_failure_load_partial(5));
fprintf(fid, "Mxy =\t%f N\n",last_ply_failure_load_partial(6));
%% Function to read data from files
% readData(fname)
% fname = file name of  csv file
% out = matrix containing data from csv file

function out = readData(fname)
    out = table2array(readtable(fname, "VariableNamingRule","preserve"));
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

