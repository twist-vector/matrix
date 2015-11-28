# Simple Matrix module for Nim

This module implements a linear algebra (matrix) class.

The algorithms in this module are not necessarily the fastest
possible or memory efficient.  Some optimizations have been
made, hoever.

## Limitations
  -  supports only dense matrices
  -  not the most memory efficient algorithms (returns new
     matrices rather than performing operations in place)
  - no "inverse" function is provided at this point

## Examples

### Creation

A "constructor" for matrix types is provided via the generic `newMatrix[T]`
function.  Two flavors are provided, one which generates a zero-filled matrix of
the requested size, and one that creates the matrix and populates it with th
provided data.
```nim
import matrix
let x = newMatrix[float](2,3,[0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
let y = newMatrix[float](3,2,[0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
let z = newMatrix[float](2,3)
echo x
echo y
echo z
```
An overloaded "stringify" operator `$` is provided to pretty-print a matrix
allowing the easy use of `echo`.  The output of the above code would be similar
to:
```
|0.0, 1.0, 2.0|
|3.0, 4.0, 5.0|

|0.0, 1.0|
|2.0, 3.0|
|4.0, 5.0|

|0.0, 0.0, 0.0|
|0.0, 0.0, 0.0|
```

Note that the data provided in the creation of both `x` and `y` is identical but
it is interpreted or reshaped differently.  The number of elements in the
supplied array must equal the total number of elements in the to-be-created
matrix or an exception will be thrown.

Convenience functions are provided for common matrices include ones, zeros, and
the identity:
```nim
import matrix
let x = ones[float](2,3)
let y = zeros[float](3,2)
let z = eye[float](3)
echo x
echo y
echo z
```
producing
```
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|

|0.0, 0.0|
|0.0, 0.0|
|0.0, 0.0|

|1.0, 0.0, 0.0|
|0.0, 1.0, 0.0|
|0.0, 0.0, 1.0|
```
Of course a matrix of other types are supported:
```nim
import matrix, complex
let w = rand[float64](2,3,5)
let x = ones[float64](2,3)
let y = zeros[int](3,2)
let z = eye[Complex](3)
echo w
echo x
echo y
echo z
```
reulting in
```
|1.982323868801377, 4.202426847057126, 1.766680486226218|
|2.23291717398272, 1.593463861559403, 4.432142166115156|

|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|

|0, 0|
|0, 0|
|0, 0|

|(1.0, 0.0), (0.0, 0.0), (0.0, 0.0)|
|(0.0, 0.0), (1.0, 0.0), (0.0, 0.0)|
|(0.0, 0.0), (0.0, 0.0), (1.0, 0.0)|
```
of course, the random values of `w` may be different.


### Math functions

The arithmetic operators are overloaded to support matrices.  As such `+`, `-`,
and `*` are available between matrices of identical types.  Division of matrices
is not (yet) supported.  In addition to the linear algebra matrix multiply `*`
and element-by-element multiply is provided via the new operator `*.` as seen
below:
```nim
import matrix
let x = ones[float](3,3)
echo x
echo (x+x)
echo (x-x)
echo (x*x)
echo (x*.x)
```
resulting in
```
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|

|2.0, 2.0, 2.0|
|2.0, 2.0, 2.0|
|2.0, 2.0, 2.0|

|0.0, 0.0, 0.0|
|0.0, 0.0, 0.0|
|0.0, 0.0, 0.0|

|3.0, 3.0, 3.0|
|3.0, 3.0, 3.0|
|3.0, 3.0, 3.0|

|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|
```
Scalar-matrix operations are supported:
```nim
import matrix
let x = ones[float](3,3)
echo x
echo x+10.0
echo 10.0+x
echo x-10.0
echo 10.0-x
echo 10.0*x
echo x*10.0
```
resulting in:
```
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|
|1.0, 1.0, 1.0|

|11.0, 11.0, 11.0|
|11.0, 11.0, 11.0|
|11.0, 11.0, 11.0|

|11.0, 11.0, 11.0|
|11.0, 11.0, 11.0|
|11.0, 11.0, 11.0|

|-9.0, -9.0, -9.0|
|-9.0, -9.0, -9.0|
|-9.0, -9.0, -9.0|

|9.0, 9.0, 9.0|
|9.0, 9.0, 9.0|
|9.0, 9.0, 9.0|

|10.0, 10.0, 10.0|
|10.0, 10.0, 10.0|
|10.0, 10.0, 10.0|

|10.0, 10.0, 10.0|
|10.0, 10.0, 10.0|
|10.0, 10.0, 10.0|
```

## Caveats and TODO
- At this point there is no matrix inversion nor is division by matrices
supported.
- All types of floats, ints and complex have been tested.  User defined types,
obviously, have not been verified.
