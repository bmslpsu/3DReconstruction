# Improvements
## Bug 1 improvement
It is able to regenerate body angles for specific frames within the original frames. There is an option to overwrite the whole body angles file or correct some specific frames in the file if the body angles file exists.  
If only need to calculate body angle, it will never overwrite frames now.
Now there is a limitation that the frames must be within the original frame range. More functions can be added if it is necessary.

## Bug 2
A mixed color problem is fixed.  
This is not a bug. It is just because the 3D model is show in a mirrored way. Red shows on the left but it is the right in fact, same for the blue. 

# Bugs
## ~~Bug 1:~~
~~When generating/regenerating new frames, all data is refreshed, which makes the hull reconstruction error
something need to be fixed after ```initialize directory``` and ```initialize(root, fly)```.~~
  
~~Because Bug 1 can cause a huge problem when running the code accidentally, the code now is temporally separated into two '.m' files.
The ```run_mask.m``` is to create masks for a video only, ```run_main.m``` is for hull reconstruction.~~

## Bug 2: 
Colors flip at ```analyze the hull``` in ```run_main.m```. Left wing should be blue but it shows red, the same issue with right wing.