#ifndef MONTJOIE_FILE_DONNEES_STRING_HXX

#include <map>

namespace Montjoie
{
  class ParameterMultistring
  {
  protected :
    // liste des notes
    map<string, int> liste_note;
    
    // parameters of strings 
    // length L, section area A, density rho, Young's modulus E, Moment of inertia I
    // Shear modulus G, Timoshenko coefficient Kprime, string tension T0
    VectReal_wp param_L, param_A, param_rho, param_E, param_I, param_G, param_Kprime, param_T0;
    // points where the string are attached on the bridge
    VectR2 points_attache;
    // if is_hammer_damped is true, R follows a law otherwise R is null
    // if is_string_damped is true, b1, b3 follow a law otherwise they are null
    bool is_hammer_damped, is_string_damped;
    
  public :
    ParameterMultistring();
    
    void InitHammerDamping(const string& verite);
    void InitStringDamping(const string& verite);

    void InitPlanCordes(istream& input_stream);
    void InitPlanCordes(const string& file_name);

    int GetNoteNumber(string note);
    string GetNoteString(int num);
    void GetParamHammer(string note, Real_wp& masse, Real_wp& exposant,
                        Real_wp& Ks, Real_wp& Rs);

    void GetParamNote(string note, Real_wp& L, Real_wp& A, Real_wp& rho, Real_wp& T0,
                      Real_wp& E, Real_wp& I, Real_wp& G, Real_wp& k_prime,
                      R3& amo, R3& amoB2);
            
    void GetParamAttache(string note, Real_wp& x0, Real_wp& y0);
    void GetParamAttache(int num, Real_wp& x0, Real_wp& y0);
    
  };
  
}

#define MONTJOIE_FILE_DONNEES_STRING_HXX
#endif
