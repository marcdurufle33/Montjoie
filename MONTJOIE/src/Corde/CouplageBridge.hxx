#ifndef MONTJOIE_FILE_COUPLAGE_BRIDGE_HXX

namespace Montjoie
{
  
  class CouplageBridge : public InputDataProblem_Base
  {
  public :

    //! oscillator coupled with the strings
    //~ Oscillator_Base oscillatorV;
    //~ Oscillator_Base oscillatorH;
    
	Oscillator_Base oscillator[5];
	Oscillator_Base oscillatorTheta;
    
    //! strings
#ifdef PIANO_STIFF_STRING
    MultiString<5, 5, WaveEquationStiffString> note;
#elif defined(PIANO_LINEAR_STIFF_STRING)
    MultiString<5, 5, WaveEquationStiffLINString> note;
#elif defined(PIANO_NONLINEAR_STIFF_STRING)
    MultiString<5, 5, WaveEquationStiffNLString> note;
#elif defined(PIANO_LINEAR_STRING)
    MultiString<5, 5, WaveEquationString> note;
#elif defined(PIANO_LINEAR_STRING_TL)
    MultiString<5, 5, WaveEquationStringLINTL> note;
#elif defined(PIANO_LINEAR_STRING_2T)
    MultiString<5, 5, WaveEquationStringLIN2T> note;
#elif defined(PIANO_NONLINEAR_STRING_2T)
    MultiString<5, 5, WaveEquationStringNL2T> note;   
#elif defined(PIANO_NONLINEAR_STIFF_STRING_2T)
    MultiString<5, 5, WaveEquationStringStiffNL2T> note;     
#elif defined(PIANO_LINEAR_STIFF_STRING_2T)
    MultiString<5, 5, WaveEquationStringStiffLIN2T> note;      
#else
    MultiString<5, 5, WaveEquationNLString> note;
#endif

    //! buffer pour les sorties de l'energie
    WriteOnTheGoWithTinyBuffer<Real_wp, 6> output_ener_all_but_note;
    string file_output_energy;
    
  public:
    CouplageBridge(const string&);
    void SetInputData(const string & keyword, const Vector<string> &param);
    void ConstructAll(const string&);
    void Display();
    
    void AdvanceStringBridge(int nt, const Real_wp& t, bool save_snapshot = true);
    
    Real_wp WriteEnergy();
    
  };

}

#define MONTJOIE_FILE_COUPLAGE_BRIDGE_HXX
#endif
