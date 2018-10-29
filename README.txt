The codes used for this project are found in the MATLAB directory. This is a running list of the important codes (hopefully updated)

---Original codes & Descriptions---

RadialRingMagnetFieldGen.m

This code requires the inner radius, outer radius, and depth of the ring magnet used. 
It is fed the r and z data points of a given mesh and spits out the magnetic field on that slice.
The code itself has comments that explain the way it works



pipeMeshForRingMagCalc.m

This code loops over different parameters for the problem (it does not know about the magnet--yet.)
It then creates cylindrical positional meshes, and the corresponding magnetic field at each point (each cartesian direction--x, y, z--has its own file).
These are saved in a directory that has a name corersponding to the important inputs to make those meshes.



flowProfile_RingMagnetOmegaCalc102318.m

This code takes in the meshes and magnetic fields created by pipeMeshForRingMagCalc.m and calculates the corresponding value of omega (angular velocity) of a ring magnet positioned next to the pipe.
This is a one-shot code, future versions (later dates) will have extra capabilities, and hopefully are included in this document.



flowProfile_RingMagnetOmegaCalc102518.m

This code searches for directories containing magnetic field mesh and spatial data, and then calculates the speed at which the flowmeter would be spinning for flow profiles with different values of a in 1-(r/R)^a shaping (volumetric flow rate conserved).
Future versions will have to be able to be more selective with which files to look through.



