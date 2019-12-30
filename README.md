
# Simulating a Jello Cube

**Features:**

- Implemented Acceleration Functionality.
- Implemented Hooke's Law, Damping, & Collision Detection.
	- Jello Cube bounces off the wall as expected.
- Implemented the Force Field.
- Renders the Inclined Plane from file, to fit within the bounding box.
- Implemented collision detection with the inclined plane.
	- Jello cube collides and reacts appropriately when in contact with the plane.


The Jello Cube is elastic, and can be bent, stretched, and squeezed.


If **L** is a vector from a mass point A to a mass point B,\
and **R** is the rest length of the spring between them,


## Hooke's Law (3D)

The elastic force exerted on A :

**F = -k<sub>Hook</sub> (|L| - R) (L / |L|)**


## Damping (3D)

Absorption of energy, thus gradually decreasing velocity.

**F = -k<sub>d</sub> ( ( (v<sub>A</sub> - v<sub>B</sub>) • L ) / |L| ) (L / |L|)** 

( d - damping coefficient )

<br/>

**Command Argument Format:**

	world/[name]

	e.g. world/jello.w
	
	
	
	
## Result

![jello simulation demo](result/jellosimulationdemo.gif)

***Note***

Different display modes to visualize the different types of springs :

* Structural Springs
* Shear Springs
* Bend Springs

Implementation of the 3  lead to a more stable, realistic simulation.



*base world code by Jernej Barbic*
