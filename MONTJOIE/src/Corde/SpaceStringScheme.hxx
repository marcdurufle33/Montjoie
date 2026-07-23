#ifndef MONTJOIE_FILE_SPACE_STRING_SCHEME_HXX

namespace Montjoie
{
  
  class Hammer;
  
  //! Classe pour la resolution d'equations de corde en 1-D avec des elements finis d'ordre eleve
  template<int nb_base, int nb_quad, class TypeEquation>
  class SpaceStringScheme
  {
  public :
    int Nx; // nombre de points principaux (donc Nx-1 elements)
    int Nx_interp; // nombre de points ou on veut calculer la solution (par interpolation)
    Real_wp point_observation; // point ou on veut calculer un sismogramme
    Real_wp Deltax; // pas d'espace
    Real_wp invDeltax; // inverse du pas d'espace
    // si vrai, on utilise les points de Gauss-Lobatto pour l'integration
    bool use_lobatto_quadrature;
    
    // buffer pour les sismogrammes (output_string : sismo sur point_observation,
    //              output_last_point : sismo sur le dernier point de la corde)
    WriteOnTheGoWithTinyBuffer<Real_wp,TypeEquation::dimension+1> output_string, output_last_point;
    // sismo pour l'ecrasement du marteau
    WriteOnTheGoWithTinyBuffer<Real_wp, 3> output_crush;
    // taille des buffers utilises au-dessus
    int size_buffer;
    string file_output_sismo_string, file_output_last_point, file_output_crush, file_output_interp;
    
    Globatto<Real_wp> base;//fonctions de base ( quadrature lobatto )
    Globatto<Real_wp> quad;//fonctions de quadrature (quadrature quad ou lobatto)
    Mesh<Dimension1> mesh;// maillage, points etc...
    
    // poids d'integration sur l'intervalle [0, 1]
    TinyVector<Real_wp, nb_quad> Weights;
    
    // La matrice des \int_0^1 phi_i phi_j
    TinyMatrix<Real_wp, Symmetric, nb_base, nb_base> MassMatrixBase; 
    // La matrice des Phi^L_i(Xi^G_j) sur le segment [0 1]
    TinyMatrix<Real_wp, General, nb_base, nb_quad> valPhi; 
    TinyVector<Real_wp, nb_base> valPhiLastPoint; // vecteur contenant Phi^L_i(1.0)
    // La matrice des d/dx(Phi^L_i)(Xi^G_j) sur le segment [0 1]
    TinyMatrix<Real_wp, General, nb_base, nb_quad> valdPhi;
    TinyVector<Real_wp, nb_base> valdPhiLastPoint; // vecteur contenant d/dx(Phi^L_i)(1.0)
    // La matrice des omega_^G_j Phi^L_i(Xi^G_j)
    TinyMatrix<Real_wp, General, nb_base, nb_quad> valPhiWeight;
    // La matrice des omega_^G_j d/dx(Phi^L_i)(Xi^G_j)
    TinyMatrix<Real_wp, General, nb_base, nb_quad> valdPhiWeight;
    
    // omega_^G_j d/dx(Phi^L_i)(Xi^G_j) d/dx(Phi^L_k)(Xi^G_j) / (\Delta x)
    TinyArray3D<Real_wp, nb_base, nb_quad, nb_base> valDoubledPhiWeight;
    // omega_^G_j Phi^L_i(Xi^G_j) d/dx(Phi^L_k)(Xi^G_j)
    TinyArray3D<Real_wp, nb_base, nb_quad, nb_base> valMixteWeight; 
    // omega_^G_j Phi^L_i(Xi^G_j) Phi^L_k(Xi^G_j) \Delta x
    TinyArray3D<Real_wp, nb_base, nb_quad, nb_base> valDoublePhiWeight; 
    
    // Le marteau est distribue en espace
    // (fonction delta(x-x_0) avec x_0 le point d'impact du marteau)
    // repartition_marteau(j) = integrale de phi_j contre delta(x-x_0)
    Vector<Real_wp> repartition_marteau; 
    Vector<Real_wp, VectSparse> repartition_marteau_sparse; // le meme en creux
    int elt_obs; // l' element contenant observation (celui de droite si observation est a cheval)
    int nb_elt_obs; // nb d'elements contenant observation (1 ou 2)
    Vector<Real_wp> observ; // valeur de chaque fct de base au point d'observation
    
    // grille d'interpolation
    GridInterpolation<Dimension1> grid_interp;
    Matrix<Real_wp, General, ArrayRowSparse> MatInterp; // MatInterp(p,l) = \phi_p (x_interp_l)
    
    
    SpaceStringScheme();
    
    void InitBuffer(const string& DOSSIER, int size, int numero_corde, bool remove_file = true,
                    bool double_prec = false, bool binary = false);
    
    void ConstructMesh(const Real_wp& L);
    
    void ConstruitImpact(Hammer & hammer);
    
    template<class GenericPb>
    void ComputeRepartitionVector(VectReal_wp& repartition,
                                  Vector<Real_wp, VectSparse>& repartition_sparse,
                                  const Real_wp& delta,
                                  const Real_wp& point_impact, GenericPb& hammer);
    
    void ConstruitInterp(const Real_wp & L);    
    void ConstruitObserv();
    
    template<class GenericPb>
    void WriteInterp(int n, GenericPb & var, int numero_corde);
    
    template<class GenericPb>
    void WriteSismo(Real_wp temp, GenericPb & var, int numero_corde);
    
    void CalculeIntPhi(const TinyVector<Real_wp, nb_quad>& valeurs,
                       TinyVector<Real_wp, nb_base> & res);
    
    template<int dimension>
    void CalculeIntPhi(const TinyMatrix<Real_wp, General, dimension, nb_quad>& valeurs,
		       TinyMatrix<Real_wp, General, dimension, nb_base> & res);
    
    void CalculeIntGrad(const TinyVector<Real_wp, nb_quad>& valeurs,
                        TinyVector<Real_wp, nb_base> & res);
    
    template<int dimension>
    void CalculeIntGrad(const TinyMatrix<Real_wp, General, dimension, nb_quad>& valeurs,
			TinyMatrix<Real_wp, General, dimension, nb_base> & res);
    
    template<int dimension, class GenericPb>
    void ComputeTM(int numero_corde, int ne, int instant, Matrix<Real_wp> & val,
		   GenericPb & var, TinyMatrix<Real_wp, General, dimension, nb_quad> & TM );
    
    template<int dimension, class GenericPb>
    void ComputedTM(int numero_corde, int ne, int instant, Matrix<Real_wp> & val,
		    GenericPb & var, TinyMatrix<Real_wp, General, dimension, nb_quad> & dTM );
    
  };
  
}

#define MONTJOIE_FILE_SPACE_STRING_SCHEME_HXX
#endif
