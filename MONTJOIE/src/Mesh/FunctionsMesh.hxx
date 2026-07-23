#ifndef MONTJOIE_FILE_FUNCTIONS_MESH_HXX

namespace Montjoie
{
  
#ifdef MONTJOIE_WITH_TWO_DIM
  void TranslateMesh(Mesh<Dimension2>& mesh, const R2& vec_u);
  
  void RotateMesh(Mesh<Dimension2>& mesh, const R2& center, const Real_wp& teta);
  
  void ScaleMesh(Mesh<Dimension2>& mesh, const R2& vec_u);
#endif
  
#ifdef MONTJOIE_WITH_THREE_DIM
  void TranslateMesh(Mesh<Dimension3>& mesh, const R3& vec_u);
  
  void RotateMesh(Mesh<Dimension3>& mesh, const R3& center, const R3& axis_, const Real_wp& teta);
  
  void ScaleMesh(Mesh<Dimension3>& mesh, const R3& vec_u);

  void ScaleVariableMesh(Mesh<Dimension3>& mesh, ScalingMeshVirtualFunction3D& fct);
#endif
  
  template<class T>
  void WriteElementMesh(const Mesh<Dimension2, T>& mesh,
			const ElementGeomReference<Dimension2>& Fb,
			const SetPoints<Dimension2>& PointsElem,
                        const string& file_name, int i);

  template<class T>
  void WriteElementMesh(const Mesh<Dimension3, T>& mesh,
			const ElementGeomReference<Dimension3>& Fb,
			const SetPoints<Dimension3>& PointsElem,
                        const string& file_name, int i);
  
  class BoundaryConditionEnum
  {
  public :
    
    //! boundary conditions
    enum {LINE_INSIDE, LINE_DIRICHLET, LINE_NEUMANN, LINE_ABSORBING,
	  LINE_HIGH_CONDUCTIVITY, LINE_IMPEDANCE, LINE_THIN_SLOT, LINE_NEIGHBOR,
	  LINE_TRANSMISSION, LINE_DTN, LINE_SUPPORTED, LINE_TRANSPARENT};

    // type of periodicity
    enum {PERIODIC_CTE, PERIODIC_THETA, PERIODIC_X, PERIODIC_Y, PERIODIC_Z,
          PERIODIC_XY, PERIODIC_XZ, PERIODIC_YZ, PERIODIC_XYZ, PERIODIC_ZTHETA};

    static int GetInteger(string);

    static int GetCompositionPeriodicity(int t1, int t2);
    static int GetCompositionPeriodicity(int t1, int t2, int t3);
    
  };
  
}

#define MONTJOIE_FILE_FUNCTIONS_MESH_HXX
#endif
