% This script reads element topology and nodal coordinate files exported from
% ANSYS and solves the "hole in an infinite plate" problem for comparison
% to ANSYS and analytical results.
%
clear; clc;

E = 1e4;                             % Young's Modulus and Poisson's ratio
nu =0.3;                             % for aluminum

% The two matrices defined below are used repeatedly and don't change from
% one calculation to the other and so are defined outside the element loop.

E_matrix = E/(1 - nu^2) * [ 1  nu  0 ; nu  1  0 ; 0  0  (1 - nu)/2 ];
B_1 = [ 1 0 0 0 ; 0 0 0 1 ; 0 1 1 0 ];

% Load element topology and nodal coordinate files from ANSYS

load element_topology_hinfplt.txt
%elem_map_nodes = element_topology_hinfplt(:,7:14);
elem_map_nodes = element_topology_hinfplt(:,7:10);
load nodal_coordinates_hinfplt.txt
number_of_nodes = max(max(elem_map_nodes));
[N_elem, number_nodes_per_elem] = size(elem_map_nodes);

xc = zeros(number_of_nodes,1);      xc_def = zeros(number_of_nodes,1);
yc = zeros(number_of_nodes,1);      yc_def = zeros(number_of_nodes,1);

for i = 1:number_of_nodes
    xc(i) = nodal_coordinates_hinfplt(i,2);
    yc(i) = nodal_coordinates_hinfplt(i,3);
end

nodes_restrained_x = zeros(number_of_nodes,1); 
nodes_restrained_y = zeros(number_of_nodes,1);

restraint_counter_x = 0; 
restraint_counter_y = 0;

for i = 1:number_of_nodes
    if xc(i) == 0.
        restraint_counter_x = restraint_counter_x + 1;
        nodes_restrained_x(restraint_counter_x) = i;

        if yc(i) == 1.00
            restraint_counter_y = restraint_counter_y + 1;
            nodes_restrained_y(restraint_counter_y) = i;
        end
    end
end
dofs_restrained = zeros(restraint_counter_x + restraint_counter_y,1);

for i = 1:restraint_counter_x
    dofs_restrained(i) = 2*nodes_restrained_x(i)-1;
end
for i = 1:restraint_counter_y
    dofs_restrained(i+restraint_counter_x,1) = 2*nodes_restrained_y(i);
end

dofs_loaded = [16 ];%  46  ;  48 ;  50 ;  52 ; 54  ; 56  ; 58  ; 60  ; 62  ; 44 ];
loads =   [ -10];% ; 2/3 ; 1/3 ; 2/3 ; 1/3 ; 2/3 ; 1/3 ; 2/3 ; 1/3 ; 2/3 ; 1/6 ];

% Allocate space for global stiffness matrix.  Note that we are making no
% attempt in this illustration to optimize storage/speed parameters by
% implementing partial solutions during the assembly.

K = zeros(2*number_of_nodes,2*number_of_nodes);

%elem_map_dofs = zeros(N_elem,16);

elem_map_dofs = zeros(N_elem,8);

for i = 1:N_elem
%   for j = 1:8
    for j = 1:4
        elem_map_dofs(i,2*j-1) = 2*elem_map_nodes(i,j) - 1;
        elem_map_dofs(i,2*j)   = 2*elem_map_nodes(i,j);
    end
end

% For comparison with ANSYS, a 2 x 2 Gaussian quadrature is used.  If you
% want to see the difference in results between the 2 x 2 and 3 x 3, comment out the
% next two lines and uncomment the two following those.     
gs(1) = -1/sqrt(3); gs(2) = 1/sqrt(3);
wt(1) = 1; wt(2) = 1;
% gs(1) =  -sqrt(0.6);  gs(2) = 0;       gs(3) = sqrt(0.6);
% wt(1) =        5/9 ;  wt(2) = 8/9;     wt(3) = 5/9;

for i = 1:N_elem
%   x_el(1:8) = xc(elem_map_nodes(i,1:8));
%   y_el(1:8) = yc(elem_map_nodes(i,1:8));

    x_el(1:4) = xc(elem_map_nodes(i,1:4));
    y_el(1:4) = yc(elem_map_nodes(i,1:4));
    
%   k_el = zeros(16,16);        % k_el is an element stiffness matrix
    k_el = zeros(8,8);        % k_el is an element stiffness matrix
    for j = 1:length(gs)
        xi = gs(j);
        for k = 1:length(gs)
            eta = gs(k);
            
            % We need to evaluate the partial derivatives of the
            % interpolation functions for use in both the Jacobian as well
            % as the B_3 component of the strain/displacement matrix.
%             p_f5_xi = -xi*(1 - eta);          p_f5_eta = -0.5*(1 - xi^2);
%             p_f6_xi = 0.5*(1 - eta^2);        p_f6_eta =  -eta*(1 + xi);
%             p_f7_xi = -xi*(1 + eta);          p_f7_eta =  0.5*(1 - xi^2);
%             p_f8_xi =-0.5*(1 - eta^2);        p_f8_eta =  -eta*(1 - xi);
            
            p_f1_xi = -0.25*(1 - eta); %- 0.5*p_f8_xi - 0.5*p_f5_xi;
            p_f2_xi =  0.25*(1 - eta); %- 0.5*p_f6_xi - 0.5*p_f5_xi;
            p_f3_xi =  0.25*(1 + eta); %- 0.5*p_f6_xi - 0.5*p_f7_xi;
            p_f4_xi = -0.25*(1 + eta); %- 0.5*p_f7_xi - 0.5*p_f8_xi;
            
            p_f1_eta = -0.25*(1 - xi); % - 0.5*p_f8_eta - 0.5*p_f5_eta;
            p_f2_eta = -0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f5_eta;
            p_f3_eta =  0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f7_eta;
            p_f4_eta =  0.25*(1 - xi); % - 0.5*p_f7_eta - 0.5*p_f8_eta;

                        p_f1_xi = -0.25*(1 - eta); 
            p_f2_xi =  0.25*(1 - eta); 
            p_f3_xi =  0.25*(1 + eta); 
            p_f4_xi = -0.25*(1 + eta); 
            
            p_f1_eta = -0.25*(1 - xi); % - 0.5*p_f8_eta - 0.5*p_f5_eta;
            p_f2_eta = -0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f5_eta;
            p_f3_eta =  0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f7_eta;
            p_f4_eta =  0.25*(1 - xi); % - 0.5*p_f7_eta - 0.5*p_f8_eta;

            
            % calculate Jacobian elements:
            p_x_xi = x_el(1)*p_f1_xi + x_el(2)*p_f2_xi + x_el(3)*p_f3_xi + x_el(4)*p_f4_xi;
 %                  + x_el(5)*p_f5_xi + x_el(6)*p_f6_xi + x_el(7)*p_f7_xi + x_el(8)*p_f8_xi;
            
            p_x_eta= x_el(1)*p_f1_eta + x_el(2)*p_f2_eta + x_el(3)*p_f3_eta + x_el(4)*p_f4_eta ;
  %                 + x_el(5)*p_f5_eta + x_el(6)*p_f6_eta + x_el(7)*p_f7_eta + x_el(8)*p_f8_eta;
               
            p_y_xi = y_el(1)*p_f1_xi + y_el(2)*p_f2_xi + y_el(3)*p_f3_xi + y_el(4)*p_f4_xi ;
  %                 + y_el(5)*p_f5_xi + y_el(6)*p_f6_xi + y_el(7)*p_f7_xi + y_el(8)*p_f8_xi;
            
            p_y_eta= y_el(1)*p_f1_eta + y_el(2)*p_f2_eta + y_el(3)*p_f3_eta + y_el(4)*p_f4_eta ;
   %                + y_el(5)*p_f5_eta + y_el(6)*p_f6_eta + y_el(7)*p_f7_eta + y_el(8)*p_f8_eta;
           
            % calculate Jacobian determinant as well as Gamma elements:
            detJ = p_x_xi * p_y_eta - p_x_eta * p_y_xi;
            Gamma_11 =  p_y_eta/detJ;  Gamma_12 = -p_y_xi/detJ;
            Gamma_21 = -p_x_eta/detJ;  Gamma_22 =  p_x_xi/detJ;
            
            % calculate B_2 and B_3 matrices, then B = B_1 * B_2 * B_3
            B_2 = [ Gamma_11  Gamma_12      0        0    ;
                    Gamma_21  Gamma_22      0        0    ;
                        0        0      Gamma_11  Gamma_12;
                        0        0      Gamma_21  Gamma_22 ];
                    
            B_3 = [ p_f1_xi 0  p_f2_xi 0 p_f3_xi  0 p_f4_xi  0 ; % p_f5_xi  0 p_f6_xi  0 p_f7_xi  0 p_f8_xi 0 ;
                   p_f1_eta 0 p_f2_eta 0 p_f3_eta 0 p_f4_eta 0 ; %p_f5_eta 0 p_f6_eta 0 p_f7_eta 0 p_f8_eta 0 ;
                      0 p_f1_xi  0 p_f2_xi  0  p_f3_xi 0  p_f4_xi ; %0  p_f5_xi 0  p_f6_xi 0  p_f7_xi 0  p_f8_xi ;
                      0 p_f1_eta 0 p_f2_eta 0 p_f3_eta 0 p_f4_eta ]; %0 %p_f5_eta 0 p_f6_eta 0 p_f7_eta 0 p_f8_eta];
            
            B = B_1*B_2*B_3;
            
            % integrate
            k_el = k_el + wt(j)*wt(k)*B'*E_matrix*B*detJ;
        end
    end
    
    % Map element stiffness matrix into global stiffness matrix:
%     for m = 1:16
%         for n = 1:16
     for m = 1:8
        for n = 1:8
            global_m = elem_map_dofs(i,m);
            global_n = elem_map_dofs(i,n);
            K(global_m,global_n) = K(global_m,global_n) + k_el(m,n);
        end
    end
    
    % repeat for next element
end

% Global stiffness matrix is now assembled.. impose boundary conditions:
for i = 1:length(dofs_restrained)
    K(dofs_restrained(i),:) = 0;
    K(:,dofs_restrained(i)) = 0;
    K(dofs_restrained(i),dofs_restrained(i)) = 1;
end

% Apply loads:
R = zeros(2*number_of_nodes,1);
for i = 1:length(dofs_loaded)
    R(dofs_loaded(i),1) = loads(i);
end

S = eye(2*number_of_nodes);
for i = 1:2*number_of_nodes
    S(i,i) = 1/sqrt(K(i,i));
end
K_S = S*K*S;
rslt = eig(K_S);
cond_nmbr = max(rslt)/min(rslt);

% Solve:
D = K\R;
XD =zeros(number_of_nodes/2,1);
YD =zeros(number_of_nodes/2,1);

for i = 1:number_of_nodes
XD(i)= D(2*i-1)
end

for i = 1:number_of_nodes
YD(i)= D(2*i)
end    


    
% For the purpose of amplifying the displacement field to inspect results,
% find a magnification factor based on the largest displacement, DMX:

dmx = zeros(number_of_nodes,1);
for i = 1:number_of_nodes
    dmx(i) = sqrt( D(2*i-1)^2 + D(2*i)^2 );
end
dmx_max = max(dmx);
dmx_mag = 0.05/dmx_max;

dmx_mag = dmx_mag*5;

for i = 1:number_of_nodes
    xc_def(i) = xc(i) + dmx_mag*D(2*i - 1);
    yc_def(i) = yc(i) + dmx_mag*D(2*i);
end
    
plot(xc,yc,'ko',xc_def,yc_def,'ro')
axis([-.5 3.1 -.5 2.1])
title('Amplified displacements for Beam problem')

% Now generate stresses at gauss points within elements.  We will use the
% extrapolation procedure described in Section 6.10 to estimate stresses at
% nodal locations.

gs(1) = -1/sqrt(3); gs(2) = 1/sqrt(3);
wt(1) = 1; wt(2) = 1;
sigma_x = zeros(length(gs),length(gs),N_elem);
sigma_y = zeros(length(gs),length(gs),N_elem);
tau_xy  = zeros(length(gs),length(gs),N_elem);
for i = 1:N_elem
    x_el(1:4) = xc(elem_map_nodes(i,1:4));
    y_el(1:4) = yc(elem_map_nodes(i,1:4));
    u_el(1:4) = D(2*elem_map_nodes(i,1:4) - ones(1,4));
    v_el(1:4) = D(2*elem_map_nodes(i,1:4));
    dofs_el = zeros(8,1);
    dofs_el(1:2:7) = u_el(1:4);
    dofs_el(2:2:8) = v_el(1:4);
    for j = 1:length(gs)
        xi = gs(j);
        for k = 1:length(gs)
            eta = gs(k);
            % We need to evaluate the partial derivatives of the
            % interpolation functions for use in both the Jacobian as well
            % as the B_3 component of the strain/displacement matrix.
%             p_f5_xi = -xi*(1 - eta);          p_f5_eta = -0.5*(1 - xi^2);
%             p_f6_xi = 0.5*(1 - eta^2);        p_f6_eta =  -eta*(1 + xi);
%             p_f7_xi = -xi*(1 + eta);          p_f7_eta =  0.5*(1 - xi^2);
%             p_f8_xi =-0.5*(1 - eta^2);        p_f8_eta =  -eta*(1 - xi);
            
            p_f1_xi = -0.25*(1 - eta); %- 0.5*p_f8_xi - 0.5*p_f5_xi;
            p_f2_xi =  0.25*(1 - eta); %- 0.5*p_f6_xi - 0.5*p_f5_xi;
            p_f3_xi =  0.25*(1 + eta); %- 0.5*p_f6_xi - 0.5*p_f7_xi;
            p_f4_xi = -0.25*(1 + eta); %- 0.5*p_f7_xi - 0.5*p_f8_xi;
            
            p_f1_eta = -0.25*(1 - xi); % - 0.5*p_f8_eta - 0.5*p_f5_eta;
            p_f2_eta = -0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f5_eta;
            p_f3_eta =  0.25*(1 + xi); % - 0.5*p_f6_eta - 0.5*p_f7_eta;
            p_f4_eta =  0.25*(1 - xi); % - 0.5*p_f7_eta - 0.5*p_f8_eta;
           
            % calculate Jacobian elements:
            p_x_xi = x_el(1)*p_f1_xi + x_el(2)*p_f2_xi + x_el(3)*p_f3_xi + x_el(4)*p_f4_xi;
 %                  + x_el(5)*p_f5_xi + x_el(6)*p_f6_xi + x_el(7)*p_f7_xi + x_el(8)*p_f8_xi;
            
            p_x_eta= x_el(1)*p_f1_eta + x_el(2)*p_f2_eta + x_el(3)*p_f3_eta + x_el(4)*p_f4_eta ;
  %                 + x_el(5)*p_f5_eta + x_el(6)*p_f6_eta + x_el(7)*p_f7_eta + x_el(8)*p_f8_eta;
               
            p_y_xi = y_el(1)*p_f1_xi + y_el(2)*p_f2_xi + y_el(3)*p_f3_xi + y_el(4)*p_f4_xi ;
  %                 + y_el(5)*p_f5_xi + y_el(6)*p_f6_xi + y_el(7)*p_f7_xi + y_el(8)*p_f8_xi;
            
            p_y_eta= y_el(1)*p_f1_eta + y_el(2)*p_f2_eta + y_el(3)*p_f3_eta + y_el(4)*p_f4_eta ;
   %                + y_el(5)*p_f5_eta + y_el(6)*p_f6_eta + y_el(7)*p_f7_eta + y_el(8)*p_f8_eta;
           
            % calculate Jacobian determinant as well as Gamma elements:
            detJ = p_x_xi * p_y_eta - p_x_eta * p_y_xi;
            Gamma_11 =  p_y_eta/detJ;  Gamma_12 = -p_y_xi/detJ;
            Gamma_21 = -p_x_eta/detJ;  Gamma_22 =  p_x_xi/detJ;
            
            % calculate B_2 and B_3 matrices, then B = B_1 * B_2 * B_3
            B_2 = [ Gamma_11  Gamma_12      0        0    ;
                    Gamma_21  Gamma_22      0        0    ;
                        0        0      Gamma_11  Gamma_12;
                        0        0      Gamma_21  Gamma_22 ];
                    
            B_3 = [ p_f1_xi 0  p_f2_xi 0 p_f3_xi  0 p_f4_xi  0 ; % p_f5_xi  0 p_f6_xi  0 p_f7_xi  0 p_f8_xi 0 ;
                   p_f1_eta 0 p_f2_eta 0 p_f3_eta 0 p_f4_eta 0 ; %p_f5_eta 0 p_f6_eta 0 p_f7_eta 0 p_f8_eta 0 ;
                      0 p_f1_xi  0 p_f2_xi  0  p_f3_xi 0  p_f4_xi ; %0  p_f5_xi 0  p_f6_xi 0  p_f7_xi 0  p_f8_xi ;
                      0 p_f1_eta 0 p_f2_eta 0 p_f3_eta 0 p_f4_eta ]; %0 %p_f5_eta 0 p_f6_eta 0 p_f7_eta 0 p_f8_eta];
            
            B = B_1*B_2*B_3;
            sigma = E_matrix*B*dofs_el;
            sigma_x(j,k,i) = sigma(1);
            sigma_y(j,k,i) = sigma(2);
            tau_xy(j,k,i)  = sigma(3);
        end
    end
end

% plot stress results on an element-by-element basis to see discontinuities
% in stresses.  The contour values for sigma_x below are used to make a
% comparison with an ANSYS plot.
cs(1) = -1 ; cs(2) = 0; cs(3) = 1;
N_cont = length(cs)^2;
x_cont = zeros(N_cont,1); y_cont = zeros(N_cont,1); 
sigx = zeros(N_cont,1); sigy = zeros(N_cont,1); t_xy = zeros(N_cont,1);

figure
%cont_vals_sigy = [ -0.1 0.25 0.6 0.95 1.3 1.65 2.0 2.35 2.7 3.05 ];
cont_vals_sigx = [ -37.032 -28.8 -20.57 -12.3 -4.11 4.11 12.34 20.57 28.8 37.032];

StressX = zeros(24,1);

for i = 1:N_elem
    x_el(1:4) = xc(elem_map_nodes(i,1:4));
    y_el(1:4) = yc(elem_map_nodes(i,1:4));
    i_counter = 1;
    for j = 1:length(cs)
        xi = cs(j);
        r  = sqrt(3)*xi;
        for k = 1:length(cs)
            eta = cs(k);
            s   = sqrt(3)*eta;
            % In addition to calculating strains, we will need to know the
            % coordinates of each point where we calculate stresses so that
            % we can make contour plots later.
%             f5 = 0.5*(1 - xi^2)*(1 - eta);
%             f6 = 0.5*(1 - eta^2)*(1 + xi);
%             f7 = 0.5*(1 - xi^2)*(1 + eta);
%             f8 = 0.5*(1 - eta^2)*(1 - xi);
            f1 = 0.25*(1 - xi)*(1 - eta);% - 0.5*f8 - 0.5*f5;
            f2 = 0.25*(1 + xi)*(1 - eta);% - 0.5*f6 - 0.5*f5;
            f3 = 0.25*(1 + xi)*(1 + eta);% - 0.5*f6 - 0.5*f7;
            f4 = 0.25*(1 - xi)*(1 + eta);% - 0.5*f7 - 0.5*f8;
  
            x_v = x_el(1)*f1 + x_el(2)*f2 + x_el(3)*f3 + x_el(4)*f4;% + x_el(5)*f5 + x_el(6)*f6 + x_el(7)*f7 + x_el(8)*f8;
            y_v = y_el(1)*f1 + y_el(2)*f2 + y_el(3)*f3 + y_el(4)*f4;% + y_el(5)*f5 + y_el(6)*f6 + y_el(7)*f7 + y_el(8)*f8;

            x_cont(i_counter) = x_v;
            y_cont(i_counter) = y_v;

            fr_1 = 0.25*(1 - r)*(1 - s);
            fr_2 = 0.25*(1 + r)*(1 - s);
            fr_3 = 0.25*(1 + r)*(1 + s);
            fr_4 = 0.25*(1 - r)*(1 + s);

            sigx(i_counter) = sigma_x(1,1,i)*fr_1 + sigma_x(1,2,i)*fr_4 + sigma_x(2,1,i)*fr_2 + sigma_x(2,2,i)*fr_3;
            sigy(i_counter) = sigma_y(1,1,i)*fr_1 + sigma_y(1,2,i)*fr_4 + sigma_y(2,1,i)*fr_2 + sigma_y(2,2,i)*fr_3;
            t_xy(i_counter) = tau_xy(1,1,i)*fr_1  + tau_xy(1,2,i)*fr_4  + tau_xy(2,1,i)*fr_2  + tau_xy(2,2,i)*fr_3;  
            i_counter = i_counter + 1;
        end
    end
    
    x1 = [ x_cont(1) x_cont(4) x_cont(7) ; x_cont(2) x_cont(5) x_cont(8) ; x_cont(3) x_cont(6) x_cont(9) ];
    x2 = [ y_cont(1) y_cont(4) y_cont(7) ; y_cont(2) y_cont(5) y_cont(8) ; y_cont(3) y_cont(6) y_cont(9) ];

    sigma_x_plot= [sigx(1) sigx(4) sigx(7); sigx(2) sigx(5) sigx(8) ; sigx(3) sigx(6) sigx(9)];

    StressX(4*i)    =   sigx(3);    StressX(4*i-1)  =   sigx(9);    StressX(4*i-2)  =   sigx(7);    StressX(4*i-3)  =   sigx(1);

    hold on
    contourf(x1,x2,sigma_x_plot,cont_vals_sigx)
    title('Sig_X_X Stress for Beam problem')

%    contourf(x1,x2,sigma_x_plot)
    caxis([ -37 37])
%    axis([0 3 0 2])
axis([-0.1 3.1 -0.1 2.1])
    %colorbar('vert')
end
%--------------------------------------------------------------------
figure
%cont_vals_sigy = [ -0.1 0.25 0.6 0.95 1.3 1.65 2.0 2.35 2.7 3.05 ];
cont_vals_sigy = [ -13.6334 -10.6 -7.57 -4.5 -1.5 1.5 4.54 7.57 10.60 13.6334];

StressY = zeros(24,1);

for i = 1:N_elem
    x_el(1:4) = xc(elem_map_nodes(i,1:4));
    y_el(1:4) = yc(elem_map_nodes(i,1:4));
    i_counter = 1;
    for j = 1:length(cs)
        xi = cs(j);
        r  = sqrt(3)*xi;
        for k = 1:length(cs)
            eta = cs(k);
            s   = sqrt(3)*eta;
            % In addition to calculating strains, we will need to know the
            % coordinates of each point where we calculate stresses so that
            % we can make contour plots later.
%             f5 = 0.5*(1 - xi^2)*(1 - eta);
%             f6 = 0.5*(1 - eta^2)*(1 + xi);
%             f7 = 0.5*(1 - xi^2)*(1 + eta);
%             f8 = 0.5*(1 - eta^2)*(1 - xi);
            f1 = 0.25*(1 - xi)*(1 - eta);% - 0.5*f8 - 0.5*f5;
            f2 = 0.25*(1 + xi)*(1 - eta);% - 0.5*f6 - 0.5*f5;
            f3 = 0.25*(1 + xi)*(1 + eta);% - 0.5*f6 - 0.5*f7;
            f4 = 0.25*(1 - xi)*(1 + eta);% - 0.5*f7 - 0.5*f8;
  
            x_v = x_el(1)*f1 + x_el(2)*f2 + x_el(3)*f3 + x_el(4)*f4;% + x_el(5)*f5 + x_el(6)*f6 + x_el(7)*f7 + x_el(8)*f8;
            y_v = y_el(1)*f1 + y_el(2)*f2 + y_el(3)*f3 + y_el(4)*f4;% + y_el(5)*f5 + y_el(6)*f6 + y_el(7)*f7 + y_el(8)*f8;

            x_cont(i_counter) = x_v;
            y_cont(i_counter) = y_v;

            fr_1 = 0.25*(1 - r)*(1 - s);
            fr_2 = 0.25*(1 + r)*(1 - s);
            fr_3 = 0.25*(1 + r)*(1 + s);
            fr_4 = 0.25*(1 - r)*(1 + s);

            sigx(i_counter) = sigma_x(1,1,i)*fr_1 + sigma_x(1,2,i)*fr_4 + sigma_x(2,1,i)*fr_2 + sigma_x(2,2,i)*fr_3;
            sigy(i_counter) = sigma_y(1,1,i)*fr_1 + sigma_y(1,2,i)*fr_4 + sigma_y(2,1,i)*fr_2 + sigma_y(2,2,i)*fr_3;
            t_xy(i_counter) = tau_xy(1,1,i)*fr_1  + tau_xy(1,2,i)*fr_4  + tau_xy(2,1,i)*fr_2  + tau_xy(2,2,i)*fr_3;  
            i_counter = i_counter + 1;
        end
    end
    
    x1 = [ x_cont(1) x_cont(4) x_cont(7) ; x_cont(2) x_cont(5) x_cont(8) ; x_cont(3) x_cont(6) x_cont(9) ];
    x2 = [ y_cont(1) y_cont(4) y_cont(7) ; y_cont(2) y_cont(5) y_cont(8) ; y_cont(3) y_cont(6) y_cont(9) ];

    sigma_y_plot= [sigy(1) sigy(4) sigy(7); sigy(2) sigy(5) sigy(8) ; sigy(3) sigy(6) sigy(9)];

    StressY(4*i)    =   sigy(3);    StressY(4*i-1)  =   sigy(9);    StressY(4*i-2)  =   sigy(7);    StressY(4*i-3)  =   sigy(1);

    hold on

    contourf(x1,x2,sigma_y_plot,cont_vals_sigy)
    title('Sig_Y_Y Stress for Beam problem')

    caxis([ -13.664 13.664])
    axis([-0.1 3.1 -0.1 2.1])
    %colorbar('vert')
end
%----------------------Figure 4 tauxy----------------------------------

figure
cont_vals_tauxy = [ -13.5093 -11.6183 -9.72 -7.83 -5.94 -4.05 -2.16  -.272 1.61 3.5];

StressXY = zeros(24,1);

for i = 1:N_elem
    x_el(1:4) = xc(elem_map_nodes(i,1:4));
    y_el(1:4) = yc(elem_map_nodes(i,1:4));
    i_counter = 1;
    for j = 1:length(cs)
        xi = cs(j);
        r  = sqrt(3)*xi;
        for k = 1:length(cs)
            eta = cs(k);
            s   = sqrt(3)*eta;
            % In addition to calculating strains, we will need to know the
            % coordinates of each point where we calculate stresses so that
            % we can make contour plots later.
%             f5 = 0.5*(1 - xi^2)*(1 - eta);
%             f6 = 0.5*(1 - eta^2)*(1 + xi);
%             f7 = 0.5*(1 - xi^2)*(1 + eta);
%             f8 = 0.5*(1 - eta^2)*(1 - xi);
            f1 = 0.25*(1 - xi)*(1 - eta);% - 0.5*f8 - 0.5*f5;
            f2 = 0.25*(1 + xi)*(1 - eta);% - 0.5*f6 - 0.5*f5;
            f3 = 0.25*(1 + xi)*(1 + eta);% - 0.5*f6 - 0.5*f7;
            f4 = 0.25*(1 - xi)*(1 + eta);% - 0.5*f7 - 0.5*f8;
  
            x_v = x_el(1)*f1 + x_el(2)*f2 + x_el(3)*f3 + x_el(4)*f4;% + x_el(5)*f5 + x_el(6)*f6 + x_el(7)*f7 + x_el(8)*f8;
            y_v = y_el(1)*f1 + y_el(2)*f2 + y_el(3)*f3 + y_el(4)*f4;% + y_el(5)*f5 + y_el(6)*f6 + y_el(7)*f7 + y_el(8)*f8;

            x_cont(i_counter) = x_v;
            y_cont(i_counter) = y_v;

            fr_1 = 0.25*(1 - r)*(1 - s);
            fr_2 = 0.25*(1 + r)*(1 - s);
            fr_3 = 0.25*(1 + r)*(1 + s);
            fr_4 = 0.25*(1 - r)*(1 + s);

            sigx(i_counter) = sigma_x(1,1,i)*fr_1 + sigma_x(1,2,i)*fr_4 + sigma_x(2,1,i)*fr_2 + sigma_x(2,2,i)*fr_3;
            sigy(i_counter) = sigma_y(1,1,i)*fr_1 + sigma_y(1,2,i)*fr_4 + sigma_y(2,1,i)*fr_2 + sigma_y(2,2,i)*fr_3;
            t_xy(i_counter) = tau_xy(1,1,i)*fr_1  + tau_xy(1,2,i)*fr_4  + tau_xy(2,1,i)*fr_2  + tau_xy(2,2,i)*fr_3;  
            i_counter = i_counter + 1;
        end
    end
    
    x1 = [ x_cont(1) x_cont(4) x_cont(7) ; x_cont(2) x_cont(5) x_cont(8) ; x_cont(3) x_cont(6) x_cont(9) ];
    x2 = [ y_cont(1) y_cont(4) y_cont(7) ; y_cont(2) y_cont(5) y_cont(8) ; y_cont(3) y_cont(6) y_cont(9) ];

    tau_xy_plot= [t_xy(1) t_xy(4) t_xy(7); t_xy(2) t_xy(5) t_xy(8) ; t_xy(3) t_xy(6) t_xy(9)];
    StressXY(4*i)    =   t_xy(3);    StressXY(4*i-1)  =   t_xy(9);    StressXY(4*i-2)  =   t_xy(7);    StressXY(4*i-3)  =   t_xy(1);

    hold on

    contourf(x1,x2,tau_xy_plot,cont_vals_tauxy)
    title('Tau_X_Y Stress for Beam problem')
%    contourf(x1,x2,sigma_x_plot)
    caxis([ -13.5 3.5])
%    axis([0 3 0 2])
axis([-0.1 3.1 -0.1 2.1])
    %colorbar('vert')
end

%----------------------------------------------------------------------


XD
YD
StressX 
StressY
StressXY

% end of program