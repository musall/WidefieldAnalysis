% This code is designed to align each animal's images to a reference image.
% The reference image has already been aligned to the Allen Atlas beforehand.
% Thus the code needs to find the pathway of reference image
% This code output a series of opts3 which contains the alignment information of respective image.
% and output a opts3_all file which contains all the alignment information.

function [] = Output_AlignToAllen(datapath, usfac, animal_num)  % animal_num only for the mice which were named "mSM + number"


for i = animal_num : animal_num      % The serial number of animal be used
    full_datapath = strcat(datapath, "\mSM", num2str(i), "\SpatialDisc");
    opts2 = [];
    
    if exist(full_datapath, 'dir') == 7
        
        
        folder_list = dir(fullfile(full_datapath));     % the name of subfolder begins from dir_list(3).name
        subfolder_num = size(folder_list, 1);
        
        for a = 3 : subfolder_num
            opts2_path = strcat(full_datapath, '\', folder_list(a).name, '\opts2.mat');
            if exist(opts2_path, 'file') == 2
                
                opts2 = load(opts2_path);
                opts2 = opts2(1).opts;
                opts2_location = opts2_path;
                ref = load(strcat(full_datapath, '\', folder_list(a).name, '\Snapshot_3.mat'));
                ref = ref(1).snap;     
                
            end
            
        end
        
        
        for n = 3 : subfolder_num     % Initial value should be 3
            
            image_folder = strcat(full_datapath, '\', folder_list(n).name);
            if exist(image_folder, 'dir') == 7
                
                output_final = [100,100,100,100];
                Greg_final = [];
                image_path_final = [];
                angle_final = 0;
                test_final = [];
                for c = 1 : 4
                    image_path = strcat(image_folder, '\Snapshot_', num2str(c), '.mat');
                    if exist(image_path, 'file') == 2
                        
                        b = load(image_path);
                        test = b(1).snap;
                        
                        [output, Greg, angle, error_curve] = Widefield_dftregistration_rotationway(ref,test,usfac);
                        
                        if output(1) < output_final(1)
                            output_final = output;
                            Greg_final = Greg;
                            image_path_final = image_path;
                            angle_final = angle;
                            
                            test_final = test;
                        end
                        
                    end
                end
                
                if isempty(Greg_final) == 1
                    disp(image_folder);
                    disp('contains no images!');
                else
                    
                    %% making a triangle marker to find the total rotation
                    % and translation
                    blank_panel = zeros(size(Greg_final));
                    blank_panel(270, 210) = 10;  blank_panel(270, 430) = 100000; 
                    blank_panel = imrotate(blank_panel, angle_final);
                    blank_panel = imtranslate(blank_panel, [output_final(4), output_final(3)]);
                    
                    if angle_final ~= 0
                        blank_panel = blank_panel((size(blank_panel, 1)/2 - size(Greg_final, 1)/2) : (size(blank_panel, 1)/2 + ...
                            size(Greg_final, 1)/2)-1, (size(blank_panel, 2)/2 - size(Greg_final, 2)/2) : (size(blank_panel, 2)/2 + ...
                            size(Greg_final, 2)/2)-1);
                    end
                    %%
                    
                    Greg_final_toAllen = alignAllenTransIm(Greg_final, opts2(1).transParams);
                    blank_panel = alignAllenTransIm(blank_panel, opts2(1).transParams);
                    blank_panel(isnan(blank_panel) == 1) = 0;                  

                    [y_10, x_10] = find(blank_panel <= 10 & blank_panel > 0.001);  % Sometimes, this command will induce many really tiny numbers,
                                                                                   % like 10^(-8). Thus using blank_panel>0 may induce bugs.
                    [y_100000, x_100000] = find(blank_panel > 10);
                    [tC_y, tC_x, angleD] = MoveJudge(x_10, y_10, x_100000, y_100000);   % !Pay attention here! The rotation & translation of image
                                                                                        % will divide the marker point into a few points. It's 
                                                                                        % pretty annoying. If there are bugs, this may be the root
                                                                                             
                    opts2_temporal = opts2;
                    opts2_temporal(1).fPath = image_path_final;
                    opts2_temporal(1).animal = strcat("mSM", num2str(i));
                    [~, d] = fileparts(fileparts(image_path_final));
                    opts2_temporal(1).rec = d;
                    opts2_temporal(1).transParams(1).angleD = angleD;
                    opts2_temporal(1).transParams(1).tC = [tC_x, tC_y]';

                    
                    opts3 = opts2_temporal;
                    folder_list(n).name(folder_list(n).name == '-') = [];
                    save((strcat(fileparts(image_path_final),'\opts3')), 'opts3');
                    opts3_all.(strcat('aligned_', folder_list(n).name)) = opts2_temporal;     % Use opts3_all first.
                    
                    try
                        for_save = alignAllenTransIm(test_final, opts3.transParams);
                        
                        savepath_local = strcat('X:\smusall\WF_Alignment\', 'mSM', num2str(i), '\');
                        imwrite(for_save, strcat(savepath_local, 'aligned', num2str(n-2), '.png'));
                        savepath_lab = strcat(fileparts(image_path_final),'\alignedtoAllen');
                        save((savepath_lab), 'Greg_final_toAllen');
                        
                    catch
                        disp(strcat('aligned_', folder_list(n).name));
                        disp(' has problem with alignAllenTransIm with the saved parameter!');
                        opts3_all.(strcat('aligned_', folder_list(n).name)) = 'Varibles do not work';
                        
                    end
                    
                end
            end
        end
        
        save((strcat(fileparts(fileparts(image_folder)), '\opts3_all')), 'opts3_all');
    end
    
end


end

%%
function [tC_y, tC_x, angleD] = MoveJudge(x_10, y_10, x_100000, y_100000)

x_10 = mode(x_10);  y_10 = mode(y_10);
x_100000 = mode(x_100000);  y_100000 = mode(y_100000);

vector1 = [y_100000 - y_10, x_100000 - x_10, 0];       % Must be 3D vector
vector2 = [0, 220, 0];    % [270 - 270, 430 - 210, 0]

angleD = atan2d(norm(cross(vector1,vector2)),dot(vector1,vector2));
tC_y = (y_10 + y_100000)/2 - 270;
tC_x = (x_10 + x_100000)/2 - 320;


end