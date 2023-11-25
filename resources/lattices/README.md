# Lattice File format .ltc

The lattice file describes the computational domain of the simulation. A .ltc file is composed of a header that contains

1. The name of the domain without any spaces

2. The number of dimensions

3. The extension of each dimension

	- 2D case:

		- width
	
		- height

and of a data part in which the domain is described as a set nodes.

Each node has a value between 0 and 2 with

- 0 referring to a fluid node

- 1 referring to a boundary node

- 2 referring to a solid node