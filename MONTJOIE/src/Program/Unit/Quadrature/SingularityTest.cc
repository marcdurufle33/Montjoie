#define MONTJOIE_WITH_TWO_DIM
#define MONTJOIE_WITH_THREE_DIM

#define MONTJOIE_WITH_NODAL_H1

#include "FiniteElement/MontjoieFiniteElement.hxx"

#include "Quadrature/SingularIntegration2D.cxx"
#include "Quadrature/SingularIntegration3D.cxx"

using namespace Montjoie;

int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);
  
  //Real_wp threshold = 1e4*epsilon_machine;
  
  VectR2 points; VectReal_wp weights;
  TriangleQuadrature::ConstructQuadrature(100, points, weights, TriangleQuadrature::QUADRATURE_TENSOR);
  
  R3 A1(0.2, 0.3, 0.24), B1(1.5, 0.8, -0.5), C1(0.3, 1.6, 0.4);
  //R3 A2(0.8, 1.2, 0.3), B2(A1), C2(-0.3, 0.76, 2.3);
  //  R3 A2(0.8, 1.2, 0.3), B2(A1), C2(B1);
  R3 A2(C1), B2(A1), C2(B1);
  Real_wp gamma = 0.5;
  R3 vec_u, vec_v, vec_w;
  vec_u = B1 - A1; vec_v = C1-A1;
  TimesProd(vec_u, vec_v, vec_w);
  Real_wp jacobT1 = Norm2(vec_w);
  
  vec_u = B2 - A2; vec_v = C2-A2;
  TimesProd(vec_u, vec_v, vec_w);
  Real_wp jacobT2 = Norm2(vec_w);
  
  Real_wp IntR(0);
  for (int i = 0; i < points.GetM(); i++)
    for (int j = 0; j < points.GetM(); j++)
      {
	R3 pt_x = A1 + points(i)(0)*(B1-A1) + points(i)(1)*(C1-A1);
	R3 pt_y = A2 + points(j)(0)*(B2-A2) + points(j)(1)*(C2-A2);
	Real_wp poids_x = weights(i)*jacobT1;
	Real_wp poids_y = weights(j)*jacobT2;
	
	IntR += pow(pt_x.Distance(pt_y), gamma)*poids_x*poids_y;
      }
  
  DISP(IntR);
  
  SingularDoubleQuadratureGalerkin_Base<Real_wp, Dimension2> var;
  
  VectR3 PtsAi(3), PtsAj(3);
  PtsAi(0) = A1; PtsAi(1) = B1; PtsAi(2) = C1;
  PtsAj(0) = A2; PtsAj(1) = B2; PtsAj(2) = C2;
  bool curved = false;
  VectReal_wp Dsi(1), Dsj(1);
  Dsi(0) = jacobT1; Dsj(0) = jacobT2;
  VectR3 PointsQuadI, PointsQuadJ, NormaleQuadI, NormaleQuadJ;
  EuclidianDistanceClass_Base<Dimension2> var_distance;
  Matrix<Real_wp> mat_elem, mat_II, mat_JJ;
  //var.ComputeElemMatrix(PtsAi, PtsAj, curved, Dsi, Dsj, PointsQuadI, PointsQuadJ,
  //			PointsQuadI, PointsQuadJ, NormaleQuadI, NormaleQuadJ,
  // NormaleQuadI, NormaleQuadJ, var_distance, mat_elem, mat_II, mat_JJ);
  
  return 0;
}
