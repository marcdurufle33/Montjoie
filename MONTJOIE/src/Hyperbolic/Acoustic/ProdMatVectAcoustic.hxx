#ifndef MONTJOIE_FILE_PROD_MAT_VECT_ACOUSTIC_HXX

namespace Montjoie
{
  
  template<class TypeEquationTime,
           class Vector1, class Vector2, class Vector3, class Vector4, class Vector5>
  void GetExtrapolationAcoustic(HyperbolicProblem<TypeEquationTime>& var_time,
				const Vector1& U, const Vector1& V, int level,
				Vector2& extrapolU, Vector3& extrapolV,
				Vector2& Uneighbor, Vector3& Vneighbor,
				Vector4& Un_quad, Vector5& Vn_quad, bool treatU, bool treatV);

  template<class TypeEquation>
  void MltAddStiffnessVectorial_AcousticDG(const Real_wp& alpha, int level,
                                           HyperbolicProblem<TypeEquation>& var,
                                           const VectReal_wp& B, const Real_wp& beta,
                                           VectReal_wp& C);

  template<class TypeEquation>
  void MltAddStiffnessScalar_AcousticDG(const Real_wp& alpha, int level,
                                        HyperbolicProblem<TypeEquation>& var,
					const VectReal_wp& B, const Real_wp& beta, VectReal_wp& C);

  template<class TypeEquation, class Vector1, class Vector2, class Vector3, class Vector4>
  void GetExtrapolationAcousticHDG(HyperbolicProblem<TypeEquation>& var, int level,
				   const VectReal_wp& U, const VectReal_wp& V,
				   Vector1& Un_quad, Vector2& Vn_quad,
				   Vector3& extrapolU, Vector4& extrapolV, VectReal_wp& Lambda);

  template<class TypeEquation, class Vector1, class Vector2, class Vector3, class Vector4>
  void MltAddStiffnessAcousticHDG(const Real_wp& alpha, int level, int level_vol,
				  HyperbolicProblem<TypeEquation>& var,
				  const VectReal_wp& U, const VectReal_wp& V, const VectReal_wp& Lambda,
				  const Vector1& Un_quad, const Vector2& Vn_quad,
				  const Vector3& extrapolU, const Vector4& extrapolV,
				  const Real_wp& beta, VectReal_wp& ProdU, VectReal_wp& ProdV);

  template<class TypeEquation, class T>
  void MltAddLambdaAcousticHDG(const T& alpha, int level,
			       HyperbolicProblem<TypeEquation>& var,
			       const Vector<T>& Lambda, Vector<T>& ProdUV);
  
}

#define MONTJOIE_FILE_PROD_MAT_VECT_ACOUSTIC_HXX
#endif
