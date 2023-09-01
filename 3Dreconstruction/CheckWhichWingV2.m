function [wing1_checked, wing2_checked]=CheckWhichWingV2(wing1, wing2, wing1Old, wing2Old)
% quick function that verifies that newly generated reconstructions have
% the same name as the old ones. For example using kmeans might cluster two
% different reconstructions in a various ways based on the size of each
% cluster.

del_centers=abs(norm(mean(wing1Old)-mean(wing1)));

del_centers_oppositeWings=abs(norm(mean(wing1Old)-mean(wing2)));

if del_centers>del_centers_oppositeWings %if true then wing 1 and wing 2 have been flipped
    wing1_checked=wing2;
    wing2_checked=wing1;
else
    wing1_checked=wing1;
    wing2_checked=wing2;
end