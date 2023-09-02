function rotation_angle=EstimateWingRotation(span2Hat,wingNormVecWael2,z_wing2)
% finds the rotation that minimizes the angle between the wing normal and a
% projection axis along the span axis of the wing
theta_fix=-pi+0.1:0.001:pi-0.1;
for ww=1:length(theta_fix)
     R_fix=makehgtform('axisrotate',span2Hat,theta_fix(ww));
     vec_fixed=R_fix*[wingNormVecWael2 ;1];
     vec_fixed=vec_fixed(1:end-1);
     vec_dot_product(ww)=dot(vec_fixed,z_wing2);
      
end
[~, index_rotation_angle]=max(vec_dot_product);
rotation_angle=theta_fix(index_rotation_angle)*180/pi;