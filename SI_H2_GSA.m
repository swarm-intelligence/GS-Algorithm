clc;
clear;
close all;

%Number of Population & Dimensions
pop_size = 100;
dim = 2;

%Iteration Condition
max_iter = 200;

%Domain of Benchmarks
from = -5.12;
to = -1*from;

%Results for n Times Execution
num_of_result = 5;
%Columns of total Result : dim,(gbest_fitness),(time)
total_result = zeros(num_of_result,dim+2);

for n=1:num_of_result
    tic;
    
    %Nfe Condition
    nfe = 1;
    max_nfe = 20000;
    
    %Initializing Random Population
    X = unifrnd(from,to,[pop_size dim]);
    
    %Calculating Fitness
    F_result = zeros(1, pop_size);
    for i = 1:pop_size
        F_result(1,i) = F1(X(i,:));
        nfe = nfe + 1;
    end
    
    %Initializing Best & Worst Fitness
    fit_best = min(F_result(1,:));
    fit_worst = max(F_result(1,:));
    
    G = zeros(1, max_iter);
	%Initializing Value of G0 & Alpha
    G0 = 0.9;
    alpha = 2;
    
    empty = zeros(1,pop_size);
    emptyX = zeros(pop_size,dim);
    
    m = empty;
    M = empty;
    
    F_ij = emptyX;
    F = emptyX;
    
    acceleration = emptyX;
    velocity = emptyX;
    
    %Starting the Main Loop
    for j = 1:max_iter
		%Calculating Gravitational Coefficient
        G(1,j) = G0 * exp((alpha*j)/max_iter);
        
		%Calculating Gravitational Mass for each Material
        m(1,:) = (F_result(1,:)-fit_worst)/(fit_best-fit_worst);
        M(1,:) = (m(1,:)/sum(m));
        
		%Computing the Force of each Material on the other Materials
        for l = 1:2:pop_size
            den = sqrt(sum((X(l+1,:)-X(l,:)).^2)) + eps;
            if den < 0
                den = den * -1;
            end
            F_ij(l,:) = G(1,j) ...
                * ((M(1,l) * M(1,l+1)) / den) ...
                * (X(l+1,:)-X(l,:));
        end
        
		%The Total Force on all Materials
        F = rand() * F_ij;
        
		%Calculating the Acceleration
        for z = 1:pop_size
            acceleration(z,:) = F(z,:) / (M(1,z) + eps);
        end
        
		%Calculating the Velocity
        velocity = rand() * velocity + acceleration;
		
		%Calculating the new Position
        X = rand() * X + velocity;
        
		%Computing the Fitness
        for o = 1:pop_size
            F_result(1,o) = F1(X(o,:));
            nfe = nfe + 1;
            if (nfe == max_nfe)%Termination Condition
                break;
            end
        end
        
		%Finding the Best and Worst Fintnesses
        [fit_best,index] = min(F_result(1,:));
        fit_worst = max(F_result(1,:));
        
		%Termination Condition
        if (nfe == max_nfe)
            break;
        end
        
    end
    
    total_result(n,1) = toc;
    total_result(n,2) = fit_best;
    total_result(n,3:end) = X(index,:);
end

min_fitness = min(total_result(:,2));
max_fitness = max(total_result(:,2));
mean_fitness = mean(total_result(:,2));
std_fitness = std(total_result(:,2));
mean_time = mean(total_result(:,1));

disp(strcat('Popsize:', num2str(pop_size), ', Dim:', num2str(dim)));
disp(strcat('mean fitness: ', num2str(mean_fitness)));
disp(strcat('max fitness: ', num2str(max_fitness)));
disp(strcat('min fitness: ', num2str(min_fitness)));
disp(strcat('std fitness: ', num2str(std_fitness)));
disp(strcat('mean time: ', num2str(mean_time)));