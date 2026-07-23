#ifndef MONTJOIE_FILE_WEDGE_DG_CLASSICAL_CXX

namespace Montjoie
{
  //! default constructor
  WedgeDgClassical::WedgeDgClassical() : WedgeClassical()
  {
    this->discontinuous_element = true;
    
    this->optimized_mass_matrix = true;
    this->sparse_mass_matrix = true;
  }
  
  
  //! constructing finite element
  void WedgeDgClassical::ConstructFiniteElement(int r, int rgeom, int rquad, int type_quad,
                                                int rtri, int r_quad,
                                                int type_surf_tri, int type_surf_quad, int gauss_z)
  {    
    WedgeClassical::ConstructFiniteElement(r, rgeom, rquad, type_quad,
                                           rtri, r_quad, type_surf_tri, type_surf_quad,
					   Globatto<Real_wp>::QUADRATURE_GAUSS);
    
  }
  
  
  //! computation of mass matrix knowing weighted jacobian on each quadrature point
  /*!
    \param[out] A computed mass matrix
    \param[in] coef values of jacobian |DF| on quadrature points
   */
  void WedgeDgClassical::IntegrateMassMatrix(Matrix<Real_wp, Symmetric, RowSymPacked> & A,
                                             const VectReal_wp & coef) const
  {
    A.Reallocate(nb_dof_loc, nb_dof_loc);
    A.Fill(0);
    int Ntri = this->nb_dof_tri;
    Real_wp vloc;
    int Nquad = this->Points2D_tri().GetM();
    for (int i = 0; i < Ntri; i++)
      for (int j = 0; j < Ntri; j++)
        for (int k = 0; k <= order; k++)
          {
            int irow = NumDofTri(i, k);
            int icol = NumDofTri(j, k);
            vloc = 0;
            for (int m = 0; m < Nquad; m++)
              vloc += ChTri(i, m)*ChTri(j, m)*coef(NumQuadTri(m, k));
            
            A(irow, icol) = vloc;
          }
  }
  
}

#define MONTJOIE_FILE_WEDGE_DG_CLASSICAL_CXX
#endif
