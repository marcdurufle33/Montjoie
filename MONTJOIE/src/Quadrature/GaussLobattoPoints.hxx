#ifndef MONTJOIE_FILE_GAUSS_LOBATTO_POINTS_HXX

namespace Montjoie
{
  template<class Vector>
  void ComputeGaussFormulas(int nb_points_gauss, Vector& points,
			    Vector& weights, bool gauss_lobatto);

  template<class real, class Vector>
  void ComputeFourierCoef(const int n, real& cz, Vector& cp,
			  Vector& dcp, Vector& ddcp);
  
  template<class real, class Vector>
  void ComputeLegendrePol_and_Derivative(int n, real theta, const real& cz,
					 const Vector& cp, const Vector& dcp,
					 const Vector& ddcp, real& pb,
					 real& dpb, real& ddpb);
  
  template<class T>
  void ComputeGaussBlendedFormulas(int nb_points, Vector<T>& points,
				   Vector<T>& weights, const T& tau);
  
} // end namespace

#define MONTJOIE_FILE_GAUSS_LOBATTO_POINTS_HXX
#endif
