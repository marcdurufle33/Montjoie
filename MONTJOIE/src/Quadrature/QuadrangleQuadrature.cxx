#ifndef MONTJOIE_FILE_QUADRANGLE_QUADRATURE_CXX

namespace Montjoie
{

  //! constructs quadrature formulas over unit square [0, 1]^2
  /*!
    \param[in] p polynomials of degree at most p must be integrated exactly
    \param[out] points 3-D integration points \xi_i
    \param[out] weights 3-D integration weights \omega_i
    \param[in] type_quad type of quadrature (see Globatto<Real_wp>
    \param[in] alpha parameter for Gauss-Jacobi rules
    \param[in] beta parameter for Gauss-Jacobi rules
    The quadrature formulas can be used to integrate :
    \int [(1-x)(1-y)]^alpha (x y)^beta f(x, y) = \sum \omega_i f(\xi_i)
   */
  void QuadrangleQuadrature
  ::ConstructQuadrature(int p, VectR2& points, VectReal_wp& weights, int type_quad,
			Real_wp alpha, Real_wp beta)
  {
    int r = p/2;
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(r, type_quad, alpha, beta);
    
    Matrix<int> NumQuad2D, coor;
    QuadrangleQuadrature::ConstructQuadrilateralNumbering(r, NumQuad2D, coor);
    
    weights.Reallocate((r+1)*(r+1));
    points.Reallocate((r+1)*(r+1));
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
	{
	  weights(NumQuad2D(i,j)) = gauss.Weights(i)*gauss.Weights(j);
	  points(NumQuad2D(i,j)) = R2(gauss.Points(i), gauss.Points(j));
	}
  }
  

  //! constructing numbering used in Montjoie for nodes of a square
  void QuadrangleQuadrature
  ::ConstructQuadrilateralNumbering(int order, Matrix<int>& NumNodes2D,
				    Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes2D.Clear();
	CoordinateNodes.Clear();
	return;
      }
    else if (order == 0)
      {
	NumNodes2D.Reallocate(1, 1);
	CoordinateNodes.Reallocate(1, 2);
	NumNodes2D(0, 0) = 0;
	CoordinateNodes(0, 0) = 0;
	CoordinateNodes(0, 1) = 0;
	return;
      }

    NumNodes2D.Reallocate(order+1, order+1);
    
    // Num Nodes makes a bijection between a tensorial numbering (i,j)
    // four vertices of the square
    NumNodes2D(0, 0) = 0;
    NumNodes2D(order,0) = 1;
    NumNodes2D(order, order) = 2;
    NumNodes2D(0, order) = 3;
    
    // four edges
    int nb = 4;
    for (int i = 1; i < order; i++)
      NumNodes2D(i, 0) = nb++;

    for (int i = 1; i < order; i++)
      NumNodes2D(order, i) = nb++;

    for (int i = 1; i < order; i++)
      NumNodes2D(order-i, order) = nb++;

    for (int i = 1; i < order; i++)
      NumNodes2D(0, order-i) = nb++;
    
    // inside the quadrilateral
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes2D(i, j) = nb++;
    
    CoordinateNodes.Reallocate((order+1)*(order+1), 2);
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	{
	  CoordinateNodes(NumNodes2D(i,j),0) = i;
	  CoordinateNodes(NumNodes2D(i,j),1) = j;
	}
  }

}

#define MONTJOIE_FILE_QUADRANGLE_QUADRATURE_CXX
#endif
