function [y, ii, jj, ll, index, len] = gety(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4, orient_x, gyro_data_x, orient_y, gyro_data_y, orient_z, gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, RFtime, time, i, j, l, m, N, factor)
    

    y = []; length = 0; len = [0]; index = [];
    
    for o = 1:1:N
        
        [y1, i, j, l, index1, len1] = getyNPVA(r_sim1, rdot_sim1, r_sim2, rdot_sim2, r_sim3, rdot_sim3, r_sim4, rdot_sim4,orient_x,gyro_data_x, orient_y,gyro_data_y, orient_z,gyro_data_z, acc_data_x, acc_data_y, acc_data_z, acc_time, gyro_time, RFtime, time, i, j, l, m+o-1, factor);
        
        if o == 1
            ii = i; jj = j; ll = l; 
        end
        y = [y;y1];
        
        length = length + len1;
        len = [len; length];

        index = [index; index1];
        
    end
  
end