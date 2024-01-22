%{

Given the x (centers) and y (counts, or normalized values, or percent
change, etc.) of some histogram, calculate th

Code assumes the bins at the edges do not have zero values

Calculate the cdf at the edges of the bins, then use linear interpolation
to insverse sample the cdf at random points from 0 to 1

%}

function random_samples = inverse_transform_sample_histogram(...
    bin_centers,bin_counts,rng_seed,N)

normalized_proportions = bin_counts/sum(bin_counts);

bin_size = mode(diff(bin_centers));
% scale up the last bin's size slight to ensure intended number of edges
buffer = 1.1;
edges = (bin_centers(1)-bin_size/2):bin_size:(bin_centers(end)+buffer*bin_size/2);

cdf_y = [0,cumsum(normalized_proportions)];

ARBITRARY_PRECISION = 6;

% Use the indices where cdf is changing since interp1 needs unique x (domain)
% points. (In other words, if you keep point pairs where y does not change over
% x, then that would interp1 with invalid input). 
valid_indices = find(round(diff(cdf_y),ARBITRARY_PRECISION));
% +1 to include the last point, since diff would otherwise cut it out
valid_indices = [valid_indices,max(valid_indices)+1];

cdf_x = edges(valid_indices);
cdf_y = cdf_y(valid_indices);

% Initialize random number generator
rng(rng_seed);

% Calculate the random samples using inverse transform sampling
random_samples = interp1(cdf_y,cdf_x,rand(N,1), 'linear', NaN);

% Assert no NaN values occurred
assert(all(~isnan(random_samples)));

%{
% Calculate new proportions of the bars such that linearly connecting the
% bars and the min and max bounds of the bars produces a distribution with
% equivalent area within any given bar but such that the distribution does
% not have aburpt transitions that a histogram
%%% Generate the 'A' matrix
A_diag_vector = (0.75)*ones(length(normalized_proportions),1);
A_diag_vector([1 end]) = 5/8;
A_off_diag_vector = 0.125*ones(length(normalized_proportions)-1,1);
A = diag(A_diag_vector) + diag(A_off_diag_vector,1) + diag(A_off_diag_vector,-1);
b = normalized_proportions';

base_vector = [0.25 0.5 0.25 zeros(1,2*length(normalized_proportions)+1-3)];
A = [];
for i = 1:length(normalized_proportions)
    A(i,:) = circshift(base_vector,(i-1)*2);
end

A_diag_vector([1 end]) = 0.75;
% A = [
%     0.75 0.25
%%% Solve the equation Ax = b, where b is the proportion vector, and x will be
%%% the interpolated proportions
% Use lsqnonneg to ensure proportions are all postive
calculated_proportions = lsqnonneg(A,b); 

calculated_proportions = (A\(normalized_proportions'))';
%%% Add zeros at the upper and lower bounds
calculated_proportions = [0, calculated_proportions, 0];
%%% Define the x values for interpolating proportions
bin_size = mode(diff(bin_centers));
edges = (bin_centers(1)-bin_size/2):bin_size:(bin_centers(end)+bin_size/2);
calculated_proportions_x_values = [edges(1),bin_centers,edges(end)];

% Sample the distribution very finely
N_fine = 1e3;
calculated_proportions_x_values_fine = linspace(min(calculated_proportions_x_values),max(calculated_proportions_x_values),N_fine);
calculated_proportions_fine = interp1(calculated_proportions_x_values,calculated_proportions,calculated_proportions_x_values_fine);

% Normalize the newly sampled proportions
normalized_calculated_proportions_fine = calculated_proportions_fine/sum(calculated_proportions_fine);


% Get random sample from the distribution using inverse transform sampling
%%% Calculate CDF
ARBITRARY_PRECISION = 6;
assert(round(sum(normalized_calculated_proportions_fine),ARBITRARY_PRECISION)==1)
cdf = round(cumsum(normalized_calculated_proportions_fine),ARBITRARY_PRECISION);
%%% Keep only values from the last point where cdf = 0 to the first point
%%% where cdf = 1; define those as the bounds. This is necessary because
%%% the sampling method below uses the inverse of the CDF, so the inverse
%%% of the CDF must have unique domain points.
post_last_zero_cdf_index = find(cdf<=0,1,'last')+1;
pre_first_unity_cdf_index = find(cdf>=1,1,'first')-1;
cdf_final = cdf(post_last_zero_cdf_index:pre_first_unity_cdf_index);
%%% Initialize random number generator
rng(rng_seed)
%%% Calculate the random samples using inverse transform sampling
random_samples  = interp1(cdf(post_last_zero_cdf_index:pre_first_unity_cdf_index),...
    calculated_proportions_x_values_fine(post_last_zero_cdf_index:pre_first_unity_cdf_index),...
    rand(N,1), 'linear', NaN);
%%% Remove any NaN values; need to fix this later on to prevent NaN values
%%% from occurring at all
if (any(isnan(random_samples)))
    n_nan_removed = sum(isnan(random_samples));
    random_samples(isnan(random_samples)) = [];
    
    warning('removing %d fiber since they had NaN diameter; remaining fibers: %d\n',...
        n_nan_removed,length(random_samples)); pause(3)
    
end
%}