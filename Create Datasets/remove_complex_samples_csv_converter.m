clear all
close all;
clc;

%%

load('New Databases/rand_data_iid_noise_v3_worse_k_i.mat');
data = new_data;
data(:,1:2) = new_data(:,1:2).*10^2;
original_size = length(data);
new_data = zeros(size(data));
complex = false;
count = 1;

for i=1:length(data)
    complex = false;
    for j=3:5   %DoA Values
        if imag(data(i,j)) ~=0
            complex = true;
            break;
        end
    end
    if complex == true
        continue;
    end
    new_data(count,:) = data(i,:);
    count = count+1;
end

last_element = find(new_data(:,1),1,'last');
new_data = new_data(1:last_element,:);
new_size = length(new_data);
(original_size-new_size)/original_size

% titles = {'x_tag' 'y_tag' 'doa_12' 'doa_23' 'D_phi_12' 'D_phi_23'};
% titles = {'x_tag' 'y_tag' 'doa_12' 'doa_34' 'D_phi_12' 'D_phi_34' 'D_phi_13' 'D_phi_14' 'D_phi_23' 'D_phi_24'};
save ("New Databases/v4/rand_data_worse_k_i_v4.mat", "new_data")

csvwrite('New Databases/v4/rand_data_worse_k_i_v4.csv', new_data)
