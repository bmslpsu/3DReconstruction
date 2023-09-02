function [vec_align, R_aligned, angle_align]=AlignVectors(vec_original,vec_rotation,vec_goal)

rot_angles=-pi/2:0.001:pi/2;
for ii=1:length(rot_angles)
    R_fix=makehgtform('axisrotate',vec_rotation,rot_angles(ii));
    vec_fixed=R_fix*[vec_original ;1];
    vec_fixed=vec_fixed(1:end-1);
    vec_dot_product(ii)=dot(vec_fixed,vec_goal);
end
[~,index_min]=max(vec_dot_product); %finds the index of the largest dot product (smallest angle)

R_aligned=makehgtform('axisrotate',vec_rotation,rot_angles(index_min));
vec_fixed=R_aligned*[vec_original ;1];
vec_align=vec_fixed(1:end-1);
angle_align=rot_angles(index_min);