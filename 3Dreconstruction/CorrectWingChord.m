function [chord_wing, chord_wing_alt]=CorrectWingChord(chord_wing,chord_wing_alt)

while 1
    right_wing_chord_input=input(['Select the following options:' newline '1) change chord sign.'  newline...
        '2) Flip chord and alt chord.' newline ...
        '3) Stop' newline]);
    
    switch right_wing_chord_input
        case 1
            chord_wing=-chord_wing;
        case 2
            chord1_resserve=chord_wing;
            chord_wing=chord_wing_alt;
            chord_wing_alt=chord1_resserve;
            clear chord1_resserve
        case 3
            break
    end
end