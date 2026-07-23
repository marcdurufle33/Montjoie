#ifndef MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX

namespace Montjoie
{
  //! classe pour la propagation du son dans l'air
  class FluideExterieur
  {
  public:
    //! density of the fluid
    Real_wp rho;
    //! c^2, where c is the velocity of sound waves in the fluid
    Real_wp c2;
    //! energy of the fluid
    Real_wp energie;
    //! if IsDofOnDomain(i) is true, the pressure dof i is not on the plate
    //! if IsDofVOnDomain(i) is true, the velocity dof i is not on the plate
    Vector<bool> IsDofOnDomain, IsDofVOnDomain;
    //! main object for the vibro-acoustics (coupling between a plate and a fluid)
    VibroAcousticSession<TimeReissnerMindlinEquation, AcousticEquation<Dimension3> > var_vibro;
    //! chronometre
    MontjoieTimer chrono;
    //! vecteurs d'iterations en temps
    VectReal_wp Pn, Pnm1, KhPn, KhPnm1, Vh_n, Vh_np1, BhVh_np1, Vh_nm1;
    //! vecteurs d'iterations en temps
    VectReal_wp KhPn_unassembled, KhPnm1_unassembled, Mh_unassembled;
    Matrix<Real_wp, General, RowSparse> Cbarre; //!< matrice de couplage Ch
    VectReal_wp Mh, Sh; //!< matrice de masse (diagonale) et matrice de damping.
    VectReal_wp InvMhPlusDtSh; //!< inverse de (Mh + dt/2 Sh)
    
    //! matrice de rigidite
    FemMatrixFreeClass<Real_wp, LaplaceEquation<Dimension3> > Kh;
    
    FluideExterieur();
    
    void SetInputData(const string & keyword, const VectString & param);
    void ConstructAcoustics();
    void CalculeGvol(const VectReal_wp& MhSurfMinus, const VectReal_wp& MhSurfPlus,
                     const VectReal_wp& ChLambda_Minus, const VectReal_wp& ChLambda_Plus,
		     VectReal_wp & Gvol_minus, VectReal_wp& Gvol_plus, Real_wp temps);
    
    Real_wp GetEnergy();
    void Advance(const VectReal_wp & Lambda, const Real_wp& temps,
                 VectReal_wp&, VectReal_wp& );
    
    void ComputeCbarre();
    void InitIndices(int nplate, int nvol);
    
  };
  
}

#define MONTJOIE_FILE_FLUIDE_EXTERIEUR_HXX
#endif
