#ifndef MONTJOIE_FILE_TINY_BLOCK_SOLVER_1D_HXX

namespace Seldon
{

  //! Class implementing 1-D finite element matrix
  template<class T, int r>
  class TinyBlockSolver1D
  {
  protected :
    //! number of elements
    int nb_elem;

    //! tridiagonal matrix to solve (after Schur complement)
    TinyBandMatrix<T, 1> mat_tridiag;

    // part of matrix involved in Schur complement
    // case where all the elements have the same mesh size
    // and the material are homogeneous such that each block is the same
    TinyMatrix<T, General, 2, r> a12;    
    TinyMatrix<T, General, r, 2> a21;    
    TinyMatrix<T, General, r, r> inv_a22;

    // case where the matrix is hetereogeneous
    Vector<TinyMatrix<T, General, 2, r> > vec_a12;
    Vector<TinyMatrix<T, General, r, 2> > vec_a21;    
    Vector<TinyMatrix<T, General, r, r> > vec_inv_a22;
    
  public :
    TinyBlockSolver1D();
    
    void Factorize(const Matrix<T, General, BandedCol>& mat, bool hg = false);

    template<class T0>
    void Solve(Vector<T0>& X);

  };

}
#define MONTJOIE_FILE_TINY_BLOCK_SOLVER_1D_HXX
#endif
