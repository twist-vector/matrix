import
  unittest,
  complex,
  math,
  matrix


suite "matrix checks":

  test "creation and access checks":
    # Check direct access to elements
    # A 2x3 matrix |  0  1  2 |
    #              | 10 11 12 |
    var x = newMatrix[float](2,3,[0.0, 1.0, 2.0, 10.0, 11.0, 12.0])
    check( x.rows() == 2 )
    check( x.cols() == 3 )
    # access...
    check( x[0,0] == 0 )
    check( x[0,1] == 1 )
    check( x[0,2] == 2 )
    check( x[1,0] == 10 )
    check( x[1,1] == 11 )
    check( x[1,2] == 12 )
    # modification...
    x[0,0] = 100.0
    check( x[0,0] == 100.0 )
    x[0,0] = 0.0
    check( x[0,0] == 0.0 )

    # verify resizing
    var t = newMatrix[float](2,3)
    check( t.rows() == 2 )
    check( t.cols() == 3 )
    t.setSize(20,10)
    check( t.rows() == 20 )
    check( t.cols() == 10 )


    # Verify transpose
    # A 3x2 matrix | 0 10 |
    #              | 1 11 |
    #              | 2 12 |
    let xt = newMatrix[float](3,2,[0.0, 10.0, 1.0, 11.0, 2.0, 12.0])
    check( transpose(x) == xt )

    # Check creation of ones, zeros, and identity matrices
    let o  = newMatrix[float](2,3,[1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    let o2 = ones[float](2,3)
    check( o == o2 )

    let z  = newMatrix[float](2,3,[0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    let z2 = zeros[float](2,3)
    check( z == z2 )

    let e  = newMatrix[float](3,3,[1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])
    let e2 = eye[float](3)
    check( e == e2 )


  test "simple math checks":
    let x    = newMatrix[float](2,3,[0.0, 1.0, 2.0, 10.0, 11.0, 12.0])
    let x2   = newMatrix[float](2,3,[0.0, 1.0, 4.0, 100.0, 121.0, 144.0])
    let xpx  = newMatrix[float](2,3,[0.0, 2.0, 4.0, 20.0, 22.0, 24.0])
    let xp10 = newMatrix[float](2,3,[10.0, 11.0, 12.0, 20.0, 21.0, 22.0])
    let xm10 = newMatrix[float](2,3,[-10.0, -9.0, -8.0, 0.0, 1.0, 2.0])
    let xt10 = newMatrix[float](2,3,[0.0, 10.0, 20.0, 100.0, 110.0, 120.0])
    let xd10 = newMatrix[float](2,3,[0.0, 0.1, 0.2, 1.0, 1.1, 1.2])
    let tmx  = newMatrix[float](2,3,[10.0, 9.0, 8.0, 0.0, -1.0, -2.0])
    let z    = newMatrix[float](2,3)
    # element-by-element multiply
    check( x*.x==x2 )
    # matrix addition
    check( x+x==xpx )
    # matrix subtraction
    check( x-x==z )
    # matrix-scalar addition and subtraction
    check( x+10.0 == xp10 )
    check( 10.0+x == xp10 )
    check( x-10.0 == xm10 )
    check( 10.0-x == tmx )
    # matrix-scalar multiply and divide
    check( x*10.0 == xt10 )
    check( 10.0*x == xt10 )
    check( x/10.0 == xd10 )
    # matrix multiply
    var y = transpose(x)
    var xy = newMatrix[float](2,2,[5.0, 35.0, 35.0, 365.0])
    var xty = x*y
    check( xty == xy )

    # also make a quick check of integer matrix
    var intmat  = newMatrix[int](2,3,[0, 1, 2, 10, 11, 12])
    var intmat2 = newMatrix[int](2,3,[0, 2, 4, 20, 22, 24])
    check( (intmat+intmat) == intmat2 )


  test "complex matrix checks":
    var cmat: Matrix[Complex]
    var cmat2: Matrix[Complex]
    cmat.setSize(3,2)
    cmat2.setSize(3,2)
    for r in countup(0,cmat.rows()-1):
      for c in countup(0,cmat.cols()-1):
        cmat[r,c] = (float(10*r), float(c))
        cmat2[r,c] = cmat[r,c] + cmat[r,c]
    check( (cmat+cmat) == cmat2 )


  test "math function checks":
    var a,b: Matrix[float]
    a.setSize(3,2)
    b.setSize(3,2)
    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        a[r,c] = 0.1*float(r*c) + 0.1
        b[r,c] = sin(a[r,c])
    check( sin(a) == b )

    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        b[r,c] = sqrt(a[r,c])
    check( sqrt(a) == b )

    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        b[r,c] = ln(a[r,c])
    check( ln(a) == b )

    for r in countup(0,a.rows()-1):
      for c in countup(0,a.cols()-1):
        b[r,c] = exp(a[r,c])
    check( exp(a) == b )


  test "QR decomposition check":
    let m: Matrix[float] = newMatrix(3,3,[1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0])
    var q, r: Matrix[float] = newMatrix[float](3,3)
    qr(m,q,r)
    assert( (q*r) =~ m )
