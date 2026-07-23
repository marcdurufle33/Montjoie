#include "FiniteElement/MontjoieFiniteElement.hxx"

using namespace Montjoie;

class MyFiniteElement : public TetrahedronReference<1>
{
   public :
   
   // number of dofs on each edge, triangle, tetrahedron, etc ?
   void ConstructNumberMap(NumberMap& map, int dg) const
  {
  }
   
   // main method used to construct the finite element object
   // r : order of approximation for basis functions
   // rgeom : order of approximation for geometry (shape functions)
   // rquad : order of quadrature
   // if rgeom is equal to 0 (and similarly rquad equal to 0)
   // we are considering that rgeom = r, and rquad = r
   // if type_quad is equal to -1, we use the default quadrature
   // formulas (probably Gauss formulas)
   // rsurf_tri, rsurf_quad : orders for boundary integrals
   // type_surf_tri, type_surf_quad : type of integration for triangles and quadrangles
   void ConstructFiniteElement(int r, int rgeom = 0, int rquad = 0, int type_quad = -1,
                               int rsurf_tri = 0, int rsurf_quad = 0, int type_surf_tri = -1,                            int type_surf_quad = -1, int gauss_z = -1)
  {
  }

   // computation of phi_i(x) for a given point x
   // the values phi_i(x) have to be placed in the array phi
   void ComputeValuesPhiRef(const R3& x, VectReal_wp& phi) const
  {
  }

   // computation of  \nabla phi_i(x) for a given point x
   void ComputeGradientPhiRef(const R3& x, VectR3& grad_phi) const
  {
  }

};


int main(int argc, char** argv)
{
  InitMontjoie(argc, argv);

  int r = 3;
  VectR2 pts((r+1)*(r+1));
  pts(0).Init(0.0, 0.0);
  pts(1).Init(1.0, 0.0);
  pts(2).Init(1.0, 1.0);
  pts(3).Init(0.0, 1.0);

  pts(4).Init(0.3, 0.0);
  pts(5).Init(0.7, 0.0);
  pts(6).Init(1.0, 0.3);
  pts(7).Init(1.0, 0.7);
  pts(8).Init(0.7, 1.0);
  pts(9).Init(0.3, 1.0);
  pts(10).Init(0.0, 0.7);
  pts(11).Init(0.0, 0.3);

  pts(12).Init(0.3, 0.3);
  pts(13).Init(0.7, 0.3);
  pts(14).Init(0.3, 0.7);
  pts(15).Init(0.7, 0.7);

  
  QuadrangleGeomReference quad;
  
  quad.SetNodalPoints(r, pts);

  R2 pt_loc(0.2, 0.4); VectReal_wp phi;
  quad.ComputeValuesPhiNodalRef(pt_loc, phi);
  
  QuadrangleHdivFirstFamily tri;
  
  tri.ConstructFiniteElement(6);

  //MyFiniteElement elt;
  //elt.ConstructFiniteElement(10);

  return FinalizeMontjoie();
}
