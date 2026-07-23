#ifndef MONTJOIE_FILE_HEXAHEDRON_QUADRATURE_CXX

namespace Montjoie
{
  //! constructs quadrature formulas over unit cube [0, 1]^3
  /*!
    \param[in] p polynomials of degree at most p must be integrated exactly
    \param[out] points 3-D integration points \xi_i
    \param[out] weights 3-D integration weights \omega_i
    \param[in] type_quad type of quadrature (see Globatto<Real_wp>
    \param[in] alpha parameter for Gauss-Jacobi rules
    \param[in] beta parameter for Gauss-Jacobi rules
    The quadrature formulas can be used to integrate :
    \int [(1-x)(1-y)(1-z)]^alpha (x y z)^beta f(x, y, z) = \sum \omega_i f(\xi_i)
   */
  void HexahedronQuadrature
  ::ConstructQuadrature(int p, VectR3& points, VectReal_wp& weights, int type_quad,
			Real_wp alpha, Real_wp beta)
  {
    int r = p/2;
    Globatto<Real_wp> gauss;
    gauss.ConstructQuadrature(r, type_quad, alpha, beta);
    
    Array3D<int> NumQuad3D;
    Matrix<int> coor;
    HexahedronQuadrature::
      ConstructHexahedralNumbering(r, NumQuad3D, coor);
    
    weights.Reallocate((r+1)*(r+1)*(r+1));
    points.Reallocate((r+1)*(r+1)*(r+1));
    
    for (int i = 0; i <= r; i++)
      for (int j = 0; j <= r; j++)
        for (int k = 0; k <= r; k++)
          {
            weights(NumQuad3D(i, j, k)) = gauss.Weights(i)*gauss.Weights(j)*gauss.Weights(k);
            points(NumQuad3D(i, j, k)).Init(gauss.Points(i), gauss.Points(j), gauss.Points(k));
          }
  }

  
  //! fills the 12 edges of the cube
  void HexahedronQuadrature::GetEdgeCube(Matrix<int>& hex_edge)
  {
    string data("0 1 \n 1 2 \n 3 2 \n 0 3 \n 0 4 \n 1 5 \n 2 6 \n 3 7 \n 4 5 \n 5 6 \n 7 6 \n 4 7");
    istringstream stream_data(data);
    
    hex_edge.ReadText(stream_data);
  }
  

  //! Numbering of nodes in a cube
  void HexahedronQuadrature
  ::ConstructHexahedralNumbering(int order, Array3D<int>& NumNodes3D,
				 Matrix<int>& CoordinateNodes)
  {
    if (order < 0)
      {
	NumNodes3D.Clear();
	CoordinateNodes.Clear();
	return;
      }
    else if (order == 0)
      {
	NumNodes3D.Reallocate(1, 1, 1);
	CoordinateNodes.Reallocate(1, 3);
	NumNodes3D(0, 0, 0) = 0;
	CoordinateNodes(0, 0) = 0;
	CoordinateNodes(0, 1) = 0;
	CoordinateNodes(0, 2) = 0;
	return;
      }
    
    NumNodes3D.Reallocate(order+1, order+1, order+1);
    NumNodes3D.Fill(-1);
    // CoordinateNodes(NumNodes3D(i,j,k),0) = i
    // pour 1 et 2 à la place de l'indice 0 CoordinateNodes renverra j et k
    CoordinateNodes.Reallocate(NumNodes3D.GetSize(), 3);
        
    int node = 0;
    // nodes on vertices
    NumNodes3D(0, 0, 0) = node++;
    NumNodes3D(order, 0, 0) = node++;
    NumNodes3D(order, order, 0) = node++;
    NumNodes3D(0, order, 0) = node++;
    NumNodes3D(0, 0, order) = node++;
    NumNodes3D(order, 0, order) = node++;
    NumNodes3D(order, order, order) = node++;
    NumNodes3D(0, order, order) = node++;
    
    for (int i = 0; i <= order; i += order)
      for (int j = 0; j <= order; j += order)
	for (int k = 0; k <= order; k += order)
	  {
	    node = NumNodes3D(i, j, k);
	    CoordinateNodes(node, 0) = i;
	    CoordinateNodes(node, 1) = j;
	    CoordinateNodes(node, 2) = k;
	  }
    
    int num_vertex_1, num_vertex_2;
    node = 8;
    Matrix<int> hex_edge; GetEdgeCube(hex_edge);
    // nodes on edges
    for (int ne = 0; ne < 12; ne++)
      {
	// we get the two extremities of the edge
	num_vertex_1 = hex_edge(ne, 0);
	num_vertex_2 = hex_edge(ne, 1);
	// we get the (i,j,k) coordinate of the first vertex
	int i1 = CoordinateNodes(num_vertex_1, 0);
	int j1 = CoordinateNodes(num_vertex_1, 1);
	int k1 = CoordinateNodes(num_vertex_1, 2);
	// the second vertex
	int i2 = CoordinateNodes(num_vertex_2, 0);
	int j2 = CoordinateNodes(num_vertex_2, 1);
	int k2 = CoordinateNodes(num_vertex_2, 2);
	
	// we are deducing the coordinates for any point on the edge
	for (int l = 1; l < order; l++)
	  {
	    int i = (l*i2 + (order-l)*i1)/order;
	    int j = (l*j2 + (order-l)*j1)/order;
	    int k = (l*k2 + (order-l)*k1)/order;
	    
	    NumNodes3D(i, j, k) = node;
	    node++;
	  }
      }
    
    // nodes on first face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(0, i, j) = node++;
    
    // nodes on second face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, 0, j) = node++;

    // nodes on third face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, j, 0) = node++;

    // nodes on fourth face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, j, order) = node++;

    // nodes on fifth face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(i, order, j) = node++;

    // nodes on sixth face
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	NumNodes3D(order, i, j) = node++;
    
    // nodes inside the hexahedron
    for (int i = 1; i < order; i++)
      for (int j = 1; j < order; j++)
	for (int k = 1; k < order; k++)
	  NumNodes3D(i, j, k) = node++;
    
    // generation of CoordinateNodes
    for (int i = 0; i <= order; i++)
      for (int j = 0; j <= order; j++)
	for (int k = 0; k <= order; k++)
	  {
	    node = NumNodes3D(i, j, k);
	    CoordinateNodes(node,0) = i;
	    CoordinateNodes(node,1) = j;
	    CoordinateNodes(node,2) = k;
	  }
  }
  
}

#define MONTJOIE_FILE_HEXAHEDRON_QUADRATURE_CXX
#endif
