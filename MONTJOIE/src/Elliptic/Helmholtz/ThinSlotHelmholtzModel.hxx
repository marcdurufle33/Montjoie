#ifndef MONTJOIE_FILE_THIN_SLOT_HELMHOLTZ_MODEL_HXX

namespace Montjoie
{
  
  //! empty class, overloaded
  template<class Dimension>
  class SlotModelParameters
  {
  };

#ifdef MONTJOIE_WITH_TWO_DIM
  //! class with slot parameters in 2-D
  template<>
  class SlotModelParameters<Dimension2>
  {
  public :
    //! available models
    /*!
      MESH1D : thin slot model with discretization of the 1-D slot 
      DTN : thin slot model with DtN operator
    */
    enum {MODEL_MESH1D, MODEL_DTN};
    
    //! type of model
    static int type_model;

    //! two extremities
    R2 ptA, ptB;
    //! for each extremity, edge [ptC, ptD] define the section of the slot
    TinyVector<R2, 2> ptC, ptD;
    //! thickness of the slot
    Real_wp epsilon;
    //! length of the slot
    Real_wp length;
    //! 1-D order of approximation
    int order; int offset_nodl; //!< offset number for dofs of the slot
    //! number of points for the 1-D mesh
    int nb_points;
    //! 1-D problem
    EllipticProblem<HelmholtzEquation1D> var1D;
    //! 1-D matrix
    Matrix<Complex_wp, Symmetric, ArrayRowSymSparse> mat_direct1D;
    //! list of edges intersecting the two sections of the slot
    IVect Slot_NumEdge;
    //! local coordinates defining the part of the edge in common with the slot
    VectReal_wp Slot_LocalPos, Slot_SecondPos;
    //! is the edge on the left or right extremity of the slot ?
    IVect Slot_Extremity;
    
    
    SlotModelParameters();
    
    void ComputeVar1D();
    
  };
  
  ostream& operator <<(ostream& out, const SlotModelParameters<Dimension2>& x);
#endif
  
  
} // namespace Montjoie

#define MONTJOIE_FILE_THIN_SLOT_HELMHOLTZ_MODEL_HXX
#endif
