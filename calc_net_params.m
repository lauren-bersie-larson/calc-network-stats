function [stats] = calc_net_params()
clear all
clc
close all
% Written by Lauren Bersie-Larson, 11/20/2018
% Put in your networks folder to know everything about them!
% Calculate volume fraction (and other network parameters) of networks to
% see if they change with network orientation
files = dir('*.mat');
stats = [];
% Sort files so that they are in numerical order
filescell = struct2cell(files);
filescell = filescell(1,:);
files2 = natsortfiles(filescell);

for i = 1:length(files2)
    filename = files2{i};
    load(filename)
    nodes = round(nodes, 10);
    tot_num_nodes = length(nodes);
    boundaries = [-0.5 0.5 -0.5 0.5 -0.5 0.5];

    num_fibs = length(fibers);
    omega = calc_orient(nodes, fibers);
    omega_xx = omega(1,1);
    
    nodes = conv_2D_2_lin(nodes);
    fibers = conv_2D_2_lin(fibers);
    avg_len = mean(calc_lens(nodes, fibers));
    total_len = sum(calc_lens(nodes, fibers));
    vol_fract = (total_len*(15e-6)*pi*((100e-9)^2))/((15e-6)^3);
    
    bnd_nodes = find_boundary_nodes(nodes,boundaries);
    int_nodes = find_int_nodes(nodes, boundaries); 
    num_bnd_nodes = length(bnd_nodes);
    num_int_nodes = length(int_nodes);
    
    fibers = conv_lin_2_2D(fibers, 2);
    connectivity = zeros(length(nodes),1);
    for j = 1:length(int_nodes)
        i = int_nodes(j);
        connectivity(i) = nnz(fibers(:) == i);
    end
    connectivity(~any(connectivity,2),:) = [];
    conn = mean(connectivity);
    
    % Calculate entropy of network
    nodes = conv_lin_2_2D(nodes, 3);
    ent = calc_entropy(nodes(int_nodes, :));
    
    stats = [stats; omega_xx, vol_fract, tot_num_nodes, num_bnd_nodes, num_int_nodes, num_fibs, conn, avg_len, ent];
    clear nodes fibers omega connectivity conn num_int_nodes avg_len ent
end

end
