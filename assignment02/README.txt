First go to the scene directory.

For the Uniform Bezier Spline Subdivision part
run ../bin/mk/02_model 01_spline.json

For the Sphere Subdivision
run ../bin/mk/02_model 02_subdivspheres.json

For the Normal Smoothing
run ../bin/mk/02_model 02_subdivspheres.json

For the Catmull-Clark Subdivision
run ../bin/mk/02_model 05_subdivmonkey.json

For the De Casteljau
run ../bin/mk/02_model 01_spline.json


For the displacement
I set a variable(displacement_mapping)in the surface struct, which is false in default. You can change it to true to see the result of displacement.
After that, run ../bin/mk/02_model 02_subdivspheres.json and ../bin/mk/02_model 03_subdivquads.json
