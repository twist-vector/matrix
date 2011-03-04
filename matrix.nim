#
#
#            Nimrod's Runtime Library
#        (c) Copyright 2011 Tom Krauss
#
#    See the file "copying.txt", included in this
#    distribution, for details about the copyright.
#
#
# This module implements a linear algebra (matrix) class.
#
# The algorithms in this module are not necessarily the fastest
# possible or memory efficient.  Some optimizations have been
# made, hoever.
#
# Limitations
# -----------
#   -  supports only dense matrices
#
#   -  not the most memory efficient algorithms (returns new
#      matrices rather than performing operations in place)
#
#

import
  math,
  complex
  
  
type 
  TMatrix*[T] = object
    transposed: bool
    dataRows:   int
    dataCols:   int
    data:       seq[T]
  
    
const
  EPS = 5.0e-10
  

proc index[T](x: TMatrix[T], r,c: int): int {.inline.} = 
  if r<0  or  r>(x.rows()-1):
    raise newException(EInvalidIndex, "matrix index out of range")
  if c<0  or  c>(x.cols()-1):
    raise newException(EInvalidIndex, "matrix index out of range")
  result = if x.transposed: c*x.dataCols+r else: r*x.dataCols+c
  

proc rows*[T](x: TMatrix[T]): int {.inline.} = 
  ## Returns the number of rows in the matrix `x`.
  result = if x.transposed: x.dataCols else: x.dataRows
  
proc cols*[T](x: TMatrix[T]): int {.inline.}  = 
  ## Returns the number of columns in the matrix `x`.
  result = if x.transposed: x.dataRows else: x.dataCols



proc matrix*[T](rows, cols: int, d: openarray[T]): TMatrix[T] = 
  ## Constructor.  Initializes the matrix by allocating memory
  ## for the data and setting the number of rows and columns
  ## and sets the data to the values specified in `d`.
  result.dataRows = rows
  result.dataCols = cols
  newSeq(result.data, rows*cols)
  if len(d)>0:
    if len(d)<(rows*cols):
      raise newException(EInvalidIndex, "insufficient data supplied in matrix constructor")

    for i in countup(0,rows*cols-1):
      result.data[i] = d[i]


proc setSize*[T](x: var TMatrix[T], rows, cols: int) = 
  ## Initializes the matrix by allocating memory
  ## for the data and setting the number of rows and columns.
  x.dataRows = rows
  x.dataCols = cols
  newSeq(x.data, rows*cols)
 

proc eye*[t](N: int): TMatrix[T] = 
  ## Returns the NxN square identity matrix.
  result = matrix[T](N,N)
  for i in countup(0,N-1):
    result[i,i] = T(1.0)

proc ones*[t](r,c: int): TMatrix[T] = 
  ## Returns a matrix with `r` rows of `c` columns which
  ## has all elements set to 1.
  result = matrix[T](r,c)
  for i in countup(0,r*c-1):
    result.data[i] = T(1.0)


proc `[,]`*[T](x: TMatrix[T], r,c: int): T = 
  ## Element access.  Returns the element at row `r` column `c`.
  result = x.data[x.index(r,c)]


proc `[,]=`*[T](x: var TMatrix[T], r,c: int, a: T) = 
  ## Sets the value of the element at row `r` column `c` to
  ## the value supplied in `a`.
  x.data[x.index(r,c)] = a


proc `$`*[T](x: TMatrix[T]): string = 
  ## "Pretty prints" the matrix.  All elements of the matrix are 
  ## included so large matrices will result in large strings.
  result = ""
  for r in countup(0,x.rows()-1):
    result = result & "|"
    for c in countup(0,x.cols()-1):
      result = result & $x[r,c]
      if c != (x.cols()-1):
        result = result & ", "
    result = result & "|\n"
  

proc transpose*[T](x: var TMatrix[T]): TMatrix[T] = 
  ## Transpose matrix (does not actually transpose data)
  result = x
  result.transposed = not x.transposed


proc `==`*[T](a: TMatrix[T], b: TMatrix[T]): bool = 
  ## Compare two matrices `x` and `y` for equality.  The
  ## matrices must be of the same size and all elements
  ## must have the same value or the return result will 
  ## be `false`.
  if (a.rows()==b.rows()) and (a.cols()==b.cols()):
    result = true
    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        if a[r,c] != b[r,c]:
          result = false
          break
      if not result:
        break
  else:
    result=false



proc `=~`*[T](a, b: TMatrix[T]): bool =
  ## Compare two matrices `a` and `b` approximately.
  if (a.rows()==b.rows()) and (a.cols()==b.cols()):
    result = true
    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        if abs(a[r,c]-b[r,c])>EPS:
          result = false
          break
      if not result:
        break
  else:
    result=false



proc `*`*[T](a: TMatrix[T], b: TMatrix[T]): TMatrix[T] = 
  ## Matrix multiply
  assert( a.cols()==b.rows() )
  result.setSize(a.rows(), b.cols())
  for i in countup(0,a.rows()-1):
    for j in countup(0,b.cols()-1):
      result[i,j] = 0.0
      for k in countup(0,a.cols()-1):
        result[i,j] = result[i,j] + a[i,k]*b[k,j]


proc `*.`*[T](a: TMatrix[T], b: TMatrix[T]): TMatrix[T] = 
  ## Element-by-element multiply
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]*b[r,c]


proc `-`*[T](a: TMatrix[T], b: TMatrix[T]): TMatrix[T] = 
  ## Element-by-element subtraction (a-b)
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]-b[r,c]

proc `-`*[T](a: TMatrix[T], b: T): TMatrix[T] = 
  ## Subtraction of type T `b` from matrix `a` (a-b)
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]-b

proc `-`*[T](a: T, b: TMatrix[T]): TMatrix[T] = 
  ## Subtraction of matrix `a` from type T `b` (a-b)
  result.setSize(b.rows(), b.cols())
  for r in countup(0,b.rows()-1):
    for c in countup(0,b.cols()-1):
      result[r,c] = a-b[r,c]


proc `+`*[T](a: TMatrix[T], b: TMatrix[T]): TMatrix[T] = 
  ## Element-by-element addition (a+b)
  assert( a.rows()==b.rows() )
  assert( a.cols()==b.cols() )
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]+b[r,c]

proc `+`*[T](a: TMatrix[T], b: T): TMatrix[T] = 
  ## Addition of type T to matrix (a+b)
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]+b

proc `+`*[T](a: T, b: TMatrix[T]): TMatrix[T] = 
  ## Addition of type T to matrix (a+b)
  result.setSize(b.rows(), b.cols())
  for r in countup(0,b.rows()-1):
    for c in countup(0,b.cols()-1):
      result[r,c] = a+b[r,c]

proc `/`*[T](a: TMatrix[T], b: T): TMatrix[T] = 
  ## Division of matrix `a` by type T `b` (a/b)
  result.setSize(a.rows(), a.cols())
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      result[r,c] = a[r,c]/b


template apply (a, b: expr): stmt = 
  result.setSize(a.rows(), a.cols())
  for i in countup(0,len(a.data)-1):
    result.data[i] = b(a.data[i])

proc sin*[T](a: TMatrix[T]):    TMatrix[T] = apply(a, sin)
proc cos*[T](a: TMatrix[T]):    TMatrix[T] = apply(a, cos)
proc tan*[T](a: TMatrix[T]):    TMatrix[T] = apply(a, tan)
proc arcsin*[T](a: TMatrix[T]): TMatrix[T] = apply(a, arcsin)
proc arccos*[T](a: TMatrix[T]): TMatrix[T] = apply(a, arccos)
proc arctan*[T](a: TMatrix[T]): TMatrix[T] = apply(a, arctan)
proc sinh*[T](a: TMatrix[T]):   TMatrix[T] = apply(a, sinh)
proc cosh*[T](a: TMatrix[T]):   TMatrix[T] = apply(a, cosh)
proc tanh*[T](a: TMatrix[T]):   TMatrix[T] = apply(a, tanh)
proc sqrt*[T](a: TMatrix[T]):   TMatrix[T] = apply(a, sqrt)
proc ln*[T](a: TMatrix[T]):     TMatrix[T] = apply(a, ln)
proc log10*[T](a: TMatrix[T]):  TMatrix[T] = apply(a, log10)
proc log2*[T](a: TMatrix[T]):   TMatrix[T] = apply(a, log2)
proc exp*[T](a: TMatrix[T]):    TMatrix[T] = apply(a, exp)
proc abs*[T](a: TMatrix[T]):    TMatrix[T] = apply(a, abs)
  



proc qr*[T](AM: TMatrix[T], Q: var TMatrix[T], R: var TMatrix[T]) =
  var A = AM
  var n = A.rows()
  var m = A.cols()
  
  var d: TMatrix[T] = matrix[T](n,m)
  for j in countup(0,n-1):
    var s = T(0)
    for i in countup(j,m-1):
      s = s + A[i,j]*A[i,j]

    s = sqrt(s) 
    if A[j,j]>T(0):
      d[j,0] = -s
    else:
      d[j,0] = s

    var fak = sqrt(s * (s+abs(A[j,j])))
    A[j,j] = A[j,j] - d[j,0]
    for k in countup(j,m-1):
      A[k,j] = A[k,j]/fak

    for i in countup((j+1),n-1):
      s = T(0)
      for k in countup(j,m-1):
        s = s + A[k,j]*A[k,i]

      for k in countup(j,m-1):
        A[k,i] = A[k,i] - A[k,j]*s

  # Reconstruct Q and R.  The diagonals of R are in d, the upper
  # triangle of R is in the upper triangle of A.  The lower 
  # triangle of A holds the Householder vectors which we use to
  # rebuild Q (via Q transpose).
  var Qt: TMatrix[T] = eye[T](n)
  for i in countup(0,min(n,m)-1):
    R[i,i] = d[i,0]

    var w: TMatrix[T] = matrix[T](n,1)
    
    for k in countup(i,n-1):
      w[k,0] = A[k,i]
    
    Qt = (eye[T](n) - w*transpose(w))*Qt
    for j in countup((i+1),n-1):
      R[i,j] = A[i,j]

  Q = transpose(Qt)





when isMainModule:
  #
  # Set up the various matricies for checking matrix
  # arithmetic.
  #
  var x:    TMatrix[float] = matrix[float](2,3,[0.0, 1.0, 2.0, 10.0, 11.0, 12.0])
  var x2:   TMatrix[float] = matrix[float](2,3,[0.0, 1.0, 4.0, 100.0, 121.0, 144.0])
  var xpx:  TMatrix[float] = matrix[float](2,3,[0.0, 2.0, 4.0, 20.0, 22.0, 24.0])
  var xp10: TMatrix[float] = matrix[float](2,3,[10.0, 11.0, 12.0, 20.0, 21.0, 22.0])
  var xm10: TMatrix[float] = matrix[float](2,3,[-10.0, -9.0, -8.0, 0.0, 1.0, 2.0])
  var tmx:  TMatrix[float] = matrix[float](2,3,[10.0, 9.0, 8.0, 0.0, -1.0, -2.0])
  var id:   TMatrix[float] = matrix[float](3,3,[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
  var o:    TMatrix[float] = matrix[float](2,3,[1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
  var z:    TMatrix[float] = matrix[float](2,3) # Automatically zeroed
  assert( x==x )
  assert( x*.x==x2 )
  assert( x+x==xpx )
  assert( x-x==z )
  assert( x+10.0==xp10 )
  assert( 10.0+x==xp10 )
  assert( x-10.0==xm10 )
  assert( 10.0-x==tmx )
  assert( eye[float](3)==id )
  assert( ones[float](2,3)==o )

  #
  # Check matrix multiply
  #
  var y = transpose(x)
  var xy = matrix[float](2,2,[5.0, 35.0, 35.0, 365.0])
  var xty = x*y
  assert( xty==xy )
  
  #
  # Check int matrix.  If we can do simple arithmetic
  # we'll assume it works for everything (not checking
  # all possible functions for all possible types).
  #
  var intmat  = matrix[int](2,3,[0, 1, 2, 10, 11, 12])
  var intmat2 = matrix[int](2,3,[0, 2, 4, 20, 22, 24])
  assert( (intmat+intmat)==intmat2 ) 
  
  
  #
  # Check complex matrix.  If we can do simple arithmetic
  # we'll assume it works for everything (not checking
  # all possible functions for all possible types).
  #
  var cmat: TMatrix[TComplex]
  var cmat2: TMatrix[TComplex]
  cmat.setSize(3,2)
  cmat2.setSize(3,2)
  for r in countup(0,cmat.rows()-1):
    for c in countup(0,cmat.cols()-1):
      cmat[r,c] = (float(10*r), float(c))
      cmat2[r,c] = cmat[r,c] + cmat[r,c]
  assert( (cmat+cmat)==cmat2 )

  # 
  # Check some of the math functions
  #
  var a,b: TMatrix[float]
  a.setSize(3,2)
  b.setSize(3,2)
  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      a[r,c] = 0.1*float(r*c) + 0.1
      b[r,c] = sin(a[r,c])
  assert( sin(a)==b )

  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      b[r,c] = sqrt(a[r,c])
  assert( sqrt(a)==b )

  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      b[r,c] = ln(a[r,c])
  assert( ln(a)==b )

  for r in countup(0,a.rows()-1):
    for c in countup(0,a.cols()-1):
      b[r,c] = exp(a[r,c])
  assert( exp(a)==b )
  
  
  var AM: TMatrix[float] = matrix(3,3,[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
  var Q, R: TMatrix[float] = matrix[float](3,3)
  qr(AM,Q,R)
  assert( (Q*R)=~AM )
  
  
  
