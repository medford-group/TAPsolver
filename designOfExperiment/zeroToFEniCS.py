#### This is a condensed and adjusted tutorial on using FEniCS and Dolfin-adjoint
####
#### Further examples can be found at:
####
#### FEniCS Basics: http://www.dolfin-adjoint.org/en/latest/documentation/maths/3-gradients.html
#### Dolfin-adjoint Basics: http://www.dolfin-adjoint.org/en/latest/documentation/examples.html
####
#### The primary difference is that all discussions are based around chemical 
#### engineering applications
#### 
#### Written by Adam C. Yonge
#### January 31, 2020
####

## Both FEniCS and fenics_adjoint (dolfin_adjoint) must be imported. fenics_adjoint overwrites
## some functions in fenics, so swapping the order would result in errors. Though you can delete 
## any of the imports below 'fenics' 'fenics_adjoint', I've generally found that importing all
## of these packages is helpful, so I usually just copy and paste them to the top of each fenics
## script I generate

from fenics import *
from fenics_adjoint import *
import mpmath
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import math as mp
import dijitso
import time
import imageio
import csv
import ast
import shutil
import sys
import os
import scipy
import pip
import pkg_resources
import ufl
from ufl import sqrt,exp,ln
import warnings

#### First the mesh of interest can be defined. Several options for mesh generation are available, 
#### and can be found on: 
#### https://fenicsproject.org/docs/dolfin/1.4.0/python/demo/documented/built-in_meshes/python/documentation.html

## Unit Interval mesh
#mesh = UnitIntervalMesh(10)

## cells in horizontal direction
#nx = 10
## cells in vertical direction
#ny = 15
## Generate the mesh
#mesh = UnitSquareMesh(nx, ny)

#x1 = Point(0,0)
#y1 = Point(3,4)
#nx = 7
#ny = 10
#diagDirection = "left/right"
## Other possible diagonal directions =  “left”, “right”, “left/right”, “crossed”
#mesh = RectangleMesh(x1, y1, nx, ny, diagonal=diagDirection)

## Unit Cube Mesh
nx = 3
ny = 3
nz = 3
mesh = UnitCubeMesh(nx, ny, nz)

## Box mesh
p1 = Point(0,0,0)
p2 = Point(2,2,2)
nx = 3
ny = 3
nz = 3
mesh = BoxMesh(p1, p2, nx, ny, nz)

## To visualize the plot you just created, uncomment out the following command.
print("Plotting a UnitIntervalMesh")
plot(mesh, title="Unit interval")
plt.show()

