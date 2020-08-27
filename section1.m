clear all

[xy_origin ima_pattern] = get_real_points_checkerboard_vmmc (9, 14 * 8, 1);

image_name = 'Images/section1/720p/im';
xy_target = zeros ([10, 2, 9]);
H = zeros ([10, 3, 3]);

figure(1)
hold on
for i = 1:10
    image = imread (sprintf('%s%d.jpg', image_name, i));
    xy_target(i, :, :) = get_user_points_vmmc (image);
    H(i, :, :) = homography_solve_vmmc (xy_origin', squeeze(xy_target(i, :, :)));
    H(i, :, :) = homography_refine_vmmc (xy_origin', squeeze(xy_target(i, :, :)), H(i, :, :));
    subplot (10, 2, 2 * i - 1)
    imshow (image)
    
    tr_image = imtransform (ima_pattern, maketform ('projective', squeeze(H(i, :, :))'), 'XData', [1 size(image, 2)], 'YData', [1 size(image, 1)]);  
    subplot (10, 2, 2 * i)
    imshow (tr_image)
end

H2 = {};
for i = 1:10
    H2{i} = squeeze(H(i, :, :));
end
A = internal_parameters_solve_vmmc (H2');