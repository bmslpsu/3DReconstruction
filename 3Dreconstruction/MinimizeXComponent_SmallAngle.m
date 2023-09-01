function [vec_min,R_min, angle_min]=MinimizeXComponent_SmallAngle(vec_original,vec_rotation)

rot_angles=-pi/2:0.001:pi/2;
for ii=1:length(rot_angles)
    R_fix=makehgtform('axisrotate',vec_rotation,rot_angles(ii));
    vec_fixed=R_fix*[vec_original ;1];
    vec_fixed=vec_fixed(1:end-1);
    vec_x_comp(ii)=abs(vec_fixed(1));
end
[~,index_min]=min(vec_x_comp); %finds the index of the smallest x-component

R_min=makehgtform('axisrotate',vec_rotation,rot_angles(index_min));
vec_fixed=R_min*[vec_original ;1];
vec_min=vec_fixed(1:end-1);
angle_min=rot_angles(index_min);