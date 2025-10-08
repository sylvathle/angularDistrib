# Angular Distribution

Minimalistic code to measure angular distribution from macro.
A sphere of 20cm is placed in the center of the world, only material is G4_Galactic to avoid any interaction.
The sphere is only used to measure the angles of the particles crossing it.

Only 4 classes are used, of which SBG4PSSphereSurfaceFlux (modified from G4PSSphereSurfaceFlux) capture the angular distribution of the incoming particles (30 bins).
The angle corresponds to the incidence over the sphere, i.e. the angle between  momentum of the particle and the point of contact with the sphere. 
e.g. if the angle is 0 it means the particle momentum is exactly perpendicular to the surface of the sphere at the point of interaction, and goes towards its center.

Also relevant is the macro (test.mac) at the root of the project where the angular distribution is defined.

The code only print the normalized angular distribution captured on the sphere as a form of a list that can be directly copied in a python script.

## Install & run

```
 git clone https://github.com/sylvathle/angularDistrib.git
 cd angularDistrib
 mkdir build
 cd build
 cmake ..
 make
 ./sim ../test.mac
```

To plot the angular distribution there is a python script *plot.py* where the list printed by the ./sim run can be copied to the variable *distrib_angles*, then run with
```
python3 plot.py
```

and an image is created.



