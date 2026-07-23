#ifndef MONTJOIE_FILE_DONNEES_STRING_CXX

namespace Montjoie
{
  
  //! default constructor
  ParameterMultistring::ParameterMultistring()
  {
    // attention les numeros commencent a 0 et vont jusqu'a 87
    liste_note.insert(pair<string,int>("A0",0));
    liste_note.insert(pair<string,int>("Ad0",1));
    liste_note.insert(pair<string,int>("Bb0",1));
    liste_note.insert(pair<string,int>("B0",2));
    
    liste_note.insert(pair<string,int>("C1",3));
    liste_note.insert(pair<string,int>("Cd1",4));
    liste_note.insert(pair<string,int>("Db1",4));
    liste_note.insert(pair<string,int>("D1",5));
    liste_note.insert(pair<string,int>("Dd1",6));
    liste_note.insert(pair<string,int>("Eb1",6));
    liste_note.insert(pair<string,int>("E1",7));
    liste_note.insert(pair<string,int>("F1",8));
    liste_note.insert(pair<string,int>("Fd1",9));
    liste_note.insert(pair<string,int>("Gb1",9));
    liste_note.insert(pair<string,int>("G1",10));
    liste_note.insert(pair<string,int>("Gd1",11));
    liste_note.insert(pair<string,int>("Ab1",11));
    liste_note.insert(pair<string,int>("A1",12));
    liste_note.insert(pair<string,int>("Ad1",13));
    liste_note.insert(pair<string,int>("Bb1",13));
    liste_note.insert(pair<string,int>("B1",14));
    
    liste_note.insert(pair<string,int>("C2",15));
    liste_note.insert(pair<string,int>("Cd2",16));
    liste_note.insert(pair<string,int>("Db2",16));
    liste_note.insert(pair<string,int>("D2",17));
    liste_note.insert(pair<string,int>("Dd2",18));
    liste_note.insert(pair<string,int>("Eb2",18));
    liste_note.insert(pair<string,int>("E2",19));
    liste_note.insert(pair<string,int>("F2",20));
    liste_note.insert(pair<string,int>("Fd2",21));
    liste_note.insert(pair<string,int>("Gb2",21));
    liste_note.insert(pair<string,int>("G2",22));
    liste_note.insert(pair<string,int>("Gd2",23));
    liste_note.insert(pair<string,int>("Ab2",23));
    liste_note.insert(pair<string,int>("A2",24));
    liste_note.insert(pair<string,int>("Ad2",25));
    liste_note.insert(pair<string,int>("Bb2",25));
    liste_note.insert(pair<string,int>("B2",26));
    
    liste_note.insert(pair<string,int>("C3",27));
    liste_note.insert(pair<string,int>("Cd3",28));
    liste_note.insert(pair<string,int>("Db3",28));
    liste_note.insert(pair<string,int>("D3",29));
    liste_note.insert(pair<string,int>("Dd3",30));
    liste_note.insert(pair<string,int>("Eb3",30));
    liste_note.insert(pair<string,int>("E3",31));
    liste_note.insert(pair<string,int>("F3",32));
    liste_note.insert(pair<string,int>("Fd3",33));
    liste_note.insert(pair<string,int>("Gb3",33));
    liste_note.insert(pair<string,int>("G3",34));
    liste_note.insert(pair<string,int>("Gd3",35));
    liste_note.insert(pair<string,int>("Ab3",35));
    liste_note.insert(pair<string,int>("A3",36));
    liste_note.insert(pair<string,int>("Ad3",37));
    liste_note.insert(pair<string,int>("Bb3",37));
    liste_note.insert(pair<string,int>("B3",38));
    
    liste_note.insert(pair<string,int>("C4",39));
    liste_note.insert(pair<string,int>("Cd4",40));
    liste_note.insert(pair<string,int>("Db4",40));
    liste_note.insert(pair<string,int>("D4",41));
    liste_note.insert(pair<string,int>("Dd4",42));
    liste_note.insert(pair<string,int>("Eb4",42));
    liste_note.insert(pair<string,int>("E4",43));
    liste_note.insert(pair<string,int>("F4",44));
    liste_note.insert(pair<string,int>("Fd4",45));
    liste_note.insert(pair<string,int>("Gb4",45));
    liste_note.insert(pair<string,int>("G4",46));
    liste_note.insert(pair<string,int>("Gd4",47));
    liste_note.insert(pair<string,int>("Ab4",47));
    liste_note.insert(pair<string,int>("A4",48));
    liste_note.insert(pair<string,int>("Ad4",49));
    liste_note.insert(pair<string,int>("Bb4",49));
    liste_note.insert(pair<string,int>("B4",50));
    
    liste_note.insert(pair<string,int>("C5",51));
    liste_note.insert(pair<string,int>("Cd5",52));
    liste_note.insert(pair<string,int>("Db5",52));
    liste_note.insert(pair<string,int>("D5",53));
    liste_note.insert(pair<string,int>("Dd5",54));
    liste_note.insert(pair<string,int>("Eb5",54));
    liste_note.insert(pair<string,int>("E5",55));
    liste_note.insert(pair<string,int>("F5",56));
    liste_note.insert(pair<string,int>("Fd5",57));
    liste_note.insert(pair<string,int>("Gb5",57));
    liste_note.insert(pair<string,int>("G5",58));
    liste_note.insert(pair<string,int>("Gd5",59));
    liste_note.insert(pair<string,int>("Ab5",59));
    liste_note.insert(pair<string,int>("A5",60));
    liste_note.insert(pair<string,int>("Ad5",61));
    liste_note.insert(pair<string,int>("Bb5",61));
    liste_note.insert(pair<string,int>("B5",62));
    
    liste_note.insert(pair<string,int>("C6",63));
    liste_note.insert(pair<string,int>("Cd6",64));
    liste_note.insert(pair<string,int>("Db6",64));
    liste_note.insert(pair<string,int>("D6",65));
    liste_note.insert(pair<string,int>("Dd6",66));
    liste_note.insert(pair<string,int>("Eb6",66));
    liste_note.insert(pair<string,int>("E6",67));
    liste_note.insert(pair<string,int>("F6",68));
    liste_note.insert(pair<string,int>("Fd6",69));
    liste_note.insert(pair<string,int>("Gb6",69));
    liste_note.insert(pair<string,int>("G6",70));
    liste_note.insert(pair<string,int>("Gd6",71));
    liste_note.insert(pair<string,int>("Ab6",71));
    liste_note.insert(pair<string,int>("A6",72));
    liste_note.insert(pair<string,int>("Ad6",73));
    liste_note.insert(pair<string,int>("Bb6",73));
    liste_note.insert(pair<string,int>("B6",74));
      
    liste_note.insert(pair<string,int>("C7",75));
    liste_note.insert(pair<string,int>("Cd7",76));
    liste_note.insert(pair<string,int>("Db7",76));
    liste_note.insert(pair<string,int>("D7",77));
    liste_note.insert(pair<string,int>("Dd7",78));
    liste_note.insert(pair<string,int>("Eb7",78));
    liste_note.insert(pair<string,int>("E7",79));
    liste_note.insert(pair<string,int>("F7",80));
    liste_note.insert(pair<string,int>("Fd7",81));
    liste_note.insert(pair<string,int>("Gb7",81));
    liste_note.insert(pair<string,int>("G7",82));
    liste_note.insert(pair<string,int>("Gd7",83));
    liste_note.insert(pair<string,int>("Ab7",83));
    liste_note.insert(pair<string,int>("A7",84));
    liste_note.insert(pair<string,int>("Ad7",85));
    liste_note.insert(pair<string,int>("Bb7",85));
    liste_note.insert(pair<string,int>("B7",86));
    
    liste_note.insert(pair<string,int>("C8",87));
    
    // points d'attache
    string data = string("0.36000000000000 1.72000000000000\n")
      +string("0.37098846698873 1.71317769037297\n")
      +string("0.38590792914654 1.70342273007615\n")
      +string("0.40369800000000 1.69031000000000\n")
      +string("0.42329829307568 1.67341438103496\n")
      +string("0.44364842190016 1.65231075407148\n")
      +string("0.46368800000000 1.62657400000000\n")
      +string("0.48254168115942 1.59597610829750\n")
      +string("0.50007428019324 1.56107750278707\n")
      +string("0.51633565217391 1.52263571587832\n")
      +string("0.53137565217391 1.48140827998089\n")
      +string("0.54524413526570 1.43815272750438\n")
      +string("0.55799095652174 1.39362659085842\n")
      +string("0.56966597101449 1.34858740245262\n")
      +string("0.58031903381643 1.30379269469661\n")
      +string("0.59000000000000 1.26000000000000\n")
      +string("0.59875872463768 1.21796685077242\n")
      +string("0.60664506280193 1.17845077942348\n")
      +string("0.61370886956522 1.14220931836280\n")
      +string("0.62000000000000 1.11000000000000\n")
      +string("0.22797100000000 1.58124000000000\n")
      +string("0.24883835553447 1.49468000783069\n")
      +string("0.27166176902954 1.41275341385488\n")
      +string("0.29598900000000 1.33523100000000\n")
      +string("0.32136780796061 1.26188354819349\n")
      +string("0.34734595242614 1.19248184036280\n")
      +string("0.37347119291137 1.12679665843536\n")
      +string("0.39929128893107 1.06459878433861\n")
      +string("0.42435400000000 1.00565900000000\n")
      +string("0.44831570954518 0.94976468295990\n")
      +string("0.47126729664258 0.89676959321044\n")
      +string("0.49340826428041 0.84654408635669\n")
      +string("0.51493811544685 0.79895851800370\n")
      +string("0.53605635313013 0.75388324375655\n")
      +string("0.55696248031844 0.71118861922029\n")
      +string("0.57785600000000 0.67074500000000\n")
      +string("0.59885623843985 0.63242052226221\n")
      +string("0.61976181501047 0.59607444441939\n")
      +string("0.64029117236115 0.56156380544545\n")
      +string("0.66016275314123 0.52874564431435\n")
      +string("0.67909500000000 0.49747700000000\n")
      +string("0.69691634952969 0.46765648341350\n")
      +string("0.71389521409411 0.43934899321455\n")
      +string("0.73041000000000 0.41266100000000\n")
      +string("0.74678752787950 0.38765431172728\n")
      +string("0.76314827566647 0.36421208579616\n")
      +string("0.77956113562021 0.34217281696695\n")
      +string("0.79609500000000 0.32137500000000\n")
      +string("0.81278859612622 0.30168596957795\n")
      +string("0.82955999156355 0.28308842007271\n")
      +string("0.84629708893778 0.26559388577849\n")
      +string("0.86288779087467 0.24921390098951\n")
      +string("0.87922000000000 0.23396000000000\n")
      +string("0.89522711165197 0.21979985996452\n")
      +string("0.91102449201852 0.20652572947907\n")
      +string("0.92677300000000 0.19388600000000\n")
      +string("0.94260259627673 0.18167065796374\n")
      +string("0.95851964864889 0.16983606972707\n")
      +string("0.97449962669660 0.15838019662687\n")
      +string("0.99051800000000 0.14730100000000\n")
      +string("1.00655184708503 0.13659582328200\n")
      +string("1.02258468226085 0.12625953830300\n")
      +string("1.03860162878245 0.11628639899182\n")
      +string("1.05458780990480 0.10667065927726\n")
      +string("1.07052834888286 0.09740657308811\n")
      +string("1.08640836897162 0.08848839435319\n")
      +string("1.10221299342605 0.07991037700129\n")
      +string("1.11792734550113 0.07166677496122\n")
      +string("1.13353654845183 0.06375184216178\n")
      +string("1.14902572553313 0.05615983253177\n")
      +string("1.16438000000000 0.04888500000000\n")
      +string("1.17958449510742 0.04192159849527\n")
      +string("1.19462433411035 0.03526388194638\n")
      +string("1.20948464026379 0.02890610428213\n")
      +string("1.22415053682269 0.02284251943133\n")
      +string("1.23860714704204 0.01706738132279\n")
      +string("1.25283959417682 0.01157494388529\n")
      +string("1.26683300148199 0.00635946104765\n")
      +string("1.28057249221253 0.00141518673867\n")
      +string("1.29404318962341 -0.00326362511284\n")
      +string("1.30723021696962 -0.00768272057809\n")
      +string("1.32011869750613 -0.01184784572828\n")
      +string("1.33269375448791 -0.01576474663459\n")
      +string("1.34494051116993 -0.01943916936823\n")
      +string("1.35684409080718 -0.02287686000040\n")
      +string("1.36838961665462 -0.02608356460228\n")
      +string("1.37956221196724 -0.02906502924508\n")
      +string("1.39034700000000 -0.03182700000000\n");
    
    istringstream input_data(data);
    points_attache.ReadText(input_data);
    
    is_hammer_damped = true;
    is_string_damped = true;
  }

  
  //! si verite = NO ou verite = FALSE, le marteau est non-amorti, sinon il est amorti
  void ParameterMultistring::InitHammerDamping(const string& verite)
  {
    if(!verite.compare("FALSE") || !verite.compare("NO"))
      {
	is_hammer_damped = false;
      }
  }
  

  //! si verite = NO ou verite = FALSE, la corde est non-amortie, sinon elle est amortie
  void ParameterMultistring::InitStringDamping(const string& verite)
  {
    // DISP(verite);
    if(!verite.compare("FALSE") || !verite.compare("NO"))
      {
	is_string_damped = false;
	cout << "string damped unactivated" << endl;
      }
  }
  
  
  //! on lit le plan de cordes a partir du flux input_stream
  void ParameterMultistring::InitPlanCordes(istream& input_stream)
  {
    string chaine, char_equal;
    int nb_cordes = 88;
    param_L.Reallocate(nb_cordes);
    param_A.Reallocate(nb_cordes);
    param_rho.Reallocate(nb_cordes);
    param_T0.Reallocate(nb_cordes);
    param_E.Reallocate(nb_cordes);
    param_I.Reallocate(nb_cordes);
    param_G.Reallocate(nb_cordes);
    param_Kprime.Reallocate(nb_cordes);
    
    // on lit chaque ligne String = L A rho T0 E I G Kprime
    for (int i = 0; i < nb_cordes; i++)
      {
	input_stream >> chaine >> char_equal >> param_L(i) >> param_A(i) >>
	  param_rho(i) >> param_T0(i) >> param_E(i) >> param_I(i) >> param_G(i) >> param_Kprime(i); 
      }
  }
  
  
  //! on lit le plan de cordes a partir d'un fichier
  void ParameterMultistring::InitPlanCordes(const string& file_name)
  {
    if (file_name == "AUTO")
      {
	// on met en dur les valeurs des parametres
	string data = string("String = 5.4728 9.518551e-08 7850 6.800859e+01 2.02e+11 7.209943e-16 8.000000e+10 0.8500\n")+
	  string("String = 5.1822 7.259579e-07 7850 5.220112e+02 2.02e+11 4.193851e-14 8.000000e+10 0.8500\n")+
	  string("String = 4.9071 1.077710e-06 7850 7.799284e+02 2.02e+11 9.242596e-14 8.000000e+10 0.8500\n")+
	  string("String = 4.6466 1.242815e-06 7850 9.052167e+02 2.02e+11 1.229144e-13 8.000000e+10 0.8500\n")+
	  string("String = 4.3999 1.288709e-06 7850 9.447204e+02 2.02e+11 1.321599e-13 8.000000e+10 0.8500\n")+
	  string("String = 4.1665 1.263688e-06 7850 9.323949e+02 2.02e+11 1.270779e-13 8.000000e+10 0.8500\n")+
	  string("String = 3.9454 1.201500e-06 7850 8.922867e+02 2.02e+11 1.148783e-13 8.000000e+10 0.8500\n")+
	  string("String = 3.7361 1.124962e-06 7850 8.409071e+02 2.02e+11 1.007084e-13 8.000000e+10 0.8500\n")+
	  string("String = 3.5380 1.048790e-06 7850 7.891102e+02 2.02e+11 8.753201e-14 8.000000e+10 0.8500\n")+
	  string("String = 3.3504 9.818027e-07 7850 7.435699e+02 2.02e+11 7.670764e-14 8.000000e+10 0.8500\n")+
	  string("String = 3.1728 9.286160e-07 7850 7.079318e+02 2.02e+11 6.862186e-14 8.000000e+10 0.8500\n")+
	  string("String = 3.0046 8.909337e-07 7850 6.837034e+02 2.02e+11 6.316564e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.8454 8.685240e-07 7850 6.709356e+02 2.02e+11 6.002799e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.6946 8.599436e-07 7850 6.687367e+02 2.02e+11 5.884778e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.5519 8.630661e-07 7850 6.756550e+02 2.02e+11 5.927592e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.4167 8.754596e-07 7850 6.899578e+02 2.02e+11 6.099052e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.2887 8.946475e-07 7850 7.098288e+02 2.02e+11 6.369335e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.1676 9.182820e-07 7850 7.335044e+02 2.02e+11 6.710305e-14 8.000000e+10 0.8500\n")+
	  string("String = 2.0528 9.442502e-07 7850 7.593617e+02 2.02e+11 7.095195e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.9442 9.707336e-07 7850 7.859711e+02 2.02e+11 7.498774e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.8413 9.962310e-07 7850 8.121223e+02 2.02e+11 7.897874e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.7439 1.019558e-06 7850 8.368320e+02 2.02e+11 8.272066e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.6517 1.039830e-06 7850 8.593378e+02 2.02e+11 8.604288e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.5644 1.056436e-06 7850 8.790830e+02 2.02e+11 8.881293e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.4817 1.069002e-06 7850 8.956970e+02 2.02e+11 9.093840e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.4033 1.077362e-06 7850 9.089715e+02 2.02e+11 9.236629e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.3292 1.081517e-06 7850 9.188365e+02 2.02e+11 9.308004e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.2590 1.081602e-06 7850 9.253363e+02 2.02e+11 9.309478e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.1925 1.077860e-06 7850 9.286065e+02 2.02e+11 9.245161e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.1295 1.070605e-06 7850 9.288529e+02 2.02e+11 9.121126e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.0699 1.060206e-06 7850 9.263324e+02 2.02e+11 8.944794e-14 8.000000e+10 0.8500\n")+
	  string("String = 1.0134 1.047060e-06 7850 9.213365e+02 2.02e+11 8.724360e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.9599 1.031580e-06 7850 9.141764e+02 2.02e+11 8.468291e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.9093 1.014173e-06 7850 9.051715e+02 2.02e+11 8.184917e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.8613 9.952367e-07 7850 8.946391e+02 2.02e+11 7.882117e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.8159 9.751455e-07 7850 8.828862e+02 2.02e+11 7.567092e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.7729 9.542466e-07 7850 8.702036e+02 2.02e+11 7.246217e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.7321 9.328546e-07 7850 8.568613e+02 2.02e+11 6.924972e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.6935 9.112498e-07 7850 8.431049e+02 2.02e+11 6.607924e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.6570 8.896764e-07 7850 8.291543e+02 2.02e+11 6.298749e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.6224 8.683429e-07 7850 8.152018e+02 2.02e+11 6.000295e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.5896 8.474227e-07 7850 8.014127e+02 2.02e+11 5.714659e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.5586 8.270561e-07 7850 7.879253e+02 2.02e+11 5.443272e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.5292 8.073522e-07 7850 7.748521e+02 2.02e+11 5.186999e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.5013 7.883915e-07 7850 7.622814e+02 2.02e+11 4.946226e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.4750 7.702286e-07 7850 7.502787e+02 2.02e+11 4.720951e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.4500 7.528953e-07 7850 7.388890e+02 2.02e+11 4.510860e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.4263 7.364031e-07 7850 7.281385e+02 2.02e+11 4.315403e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.4039 7.207462e-07 7850 7.180369e+02 2.02e+11 4.133852e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.3827 7.059045e-07 7850 7.085796e+02 2.02e+11 3.965354e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.3626 6.918456e-07 7850 6.997496e+02 2.02e+11 3.808979e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.3435 6.785280e-07 7850 6.915196e+02 2.02e+11 3.663749e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.3255 6.659025e-07 7850 6.838537e+02 2.02e+11 3.528673e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.3084 6.539148e-07 7850 6.767096e+02 2.02e+11 3.402769e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2922 6.425071e-07 7850 6.700396e+02 2.02e+11 3.285080e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2769 6.316194e-07 7850 6.637924e+02 2.02e+11 3.174688e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2624 6.211915e-07 7850 6.579147e+02 2.02e+11 3.070726e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2486 6.111634e-07 7850 6.523518e+02 2.02e+11 2.972383e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2356 6.014771e-07 7850 6.470490e+02 2.02e+11 2.878912e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2232 5.920767e-07 7850 6.419526e+02 2.02e+11 2.789627e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2115 5.829094e-07 7850 6.370102e+02 2.02e+11 2.703910e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.2004 5.739260e-07 7850 6.321718e+02 2.02e+11 2.621211e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1899 5.650812e-07 7850 6.273899e+02 2.02e+11 2.541042e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1800 5.563338e-07 7850 6.226204e+02 2.02e+11 2.462981e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1706 5.476469e-07 7850 6.178223e+02 2.02e+11 2.386664e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1616 5.389879e-07 7850 6.129584e+02 2.02e+11 2.311789e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1532 5.303287e-07 7850 6.079953e+02 2.02e+11 2.238105e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1452 5.216454e-07 7850 6.029034e+02 2.02e+11 2.165414e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1376 5.129181e-07 7850 5.976570e+02 2.02e+11 2.093564e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1304 5.041310e-07 7850 5.922340e+02 2.02e+11 2.022446e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1236 4.952719e-07 7850 5.866163e+02 2.02e+11 1.951989e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1171 4.863320e-07 7850 5.807892e+02 2.02e+11 1.882157e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1110 4.773061e-07 7850 5.747415e+02 2.02e+11 1.812943e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.1052 4.681914e-07 7850 5.684653e+02 2.02e+11 1.744364e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0997 4.589883e-07 7850 5.619556e+02 2.02e+11 1.676461e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0945 4.496993e-07 7850 5.552105e+02 2.02e+11 1.609291e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0896 4.403291e-07 7850 5.482304e+02 2.02e+11 1.542925e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0849 4.308842e-07 7850 5.410184e+02 2.02e+11 1.477445e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0805 4.213728e-07 7850 5.335796e+02 2.02e+11 1.412938e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0763 4.118043e-07 7850 5.259211e+02 2.02e+11 1.349497e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0723 4.021895e-07 7850 5.180517e+02 2.02e+11 1.287216e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0686 3.925397e-07 7850 5.099817e+02 2.02e+11 1.226189e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0650 3.828671e-07 7850 5.017229e+02 2.02e+11 1.166504e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0616 3.731845e-07 7850 4.932878e+02 2.02e+11 1.108249e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0584 3.635048e-07 7850 4.846901e+02 2.02e+11 1.051503e-14 8.000000e+10 0.8500\n")+
	  string("String = 0.0554 3.538410e-07 7850 4.759442e+02 2.02e+11 9.963376e-15 8.000000e+10 0.8500\n")+
	  string("String = 0.0525 3.442065e-07 7850 4.670650e+02 2.02e+11 9.428188e-15 8.000000e+10 0.8500\n")+
	  string("String = 0.0498 3.346141e-07 7850 4.580679e+02 2.02e+11 8.910020e-15 8.000000e+10 0.8500\n");
	
	istringstream input_stream(data);
	InitPlanCordes(input_stream);
      }
    else
      {
	// les valeurs sont lues dans un fichier
	ifstream file_in(file_name.data());
	if (file_in.is_open())
	  {
	    InitPlanCordes(file_in);
	    file_in.close();
	  }
	else
	  {
	    cout << "Erreur de lecture du fichier " << file_name << endl;
	    abort();
	  }
      }
  }
  
  
  //! retourne un numero de note a partir de son identifiant
  /*!
    \param[in] note identifiant de la note (par exemple A0, C2, Dd5, Eb4, etc)
  */
  int ParameterMultistring::GetNoteNumber(string note)
  {
    return liste_note[note];
  }
  
  
  //! retourne la chaine correspondant a une note
  /*!
    \param[in] num numero de la note 
    retourne l'identifiant de la note (par exemple A0, C2, Dd5, Eb4, etc)
  */
  string ParameterMultistring::GetNoteString(int num)
  {
    map<string, int>::iterator it;
    for (it = liste_note.begin(); it != liste_note.end(); it++)
      if ( (*it).second == num)
	return (*it).first;
    
    return string("A0");
  }
  
  
  //! calcul des parametres du marteau pour une note donnee
  /*!
    \param[in] note identifiant de la note (par exemple A0, C2, Dd5, Eb4, etc)
    \param[out] masse masse du marteau
    \param[out] exposant exposant du marteau (p)
    \param[out] Ks rigidite du marteau
    \param[out] Rs amortissement du marteau
     */
  void ParameterMultistring::GetParamHammer(string note, Real_wp& masse, Real_wp& exposant,
					    Real_wp& Ks, Real_wp& Rs)
  {
    // on recupere le numero de la note
    int num = GetNoteNumber(note) + 1;
    
    // utilisation de lois polynomiales :
    masse = -6.2348e-5*num + 0.0112;
    exposant = 2.4295e-4*num*num - 0.007703*num + 2.337;
    Ks = pow(10, 5.3097e-2*num + 7.6425);
    cout << "Hammer Damping : " << is_hammer_damped << endl;
    if ( is_hammer_damped )
      {
	Rs = Ks*pow(10, -0.04366*num - 2.294);
      }
    else
      {
	Rs = 0.0;
      }
  }
  
  
  //! calcul des parametres des cordes pour une note donnee
  /*!
    \param[in] note identifiant de la note (par exemple A0, C2, Dd5, Eb4, etc)
    \param[out] L longueur de la corde
    \param[out] A aire de la section de la corde
    \param[out] rho masse volumique de la corde
    \param[out] T0 tension de la corde
    \param[out] E module de Young
    \param[out] I moment d'inertie
    \param[out] G module de cisaillement
    \param[out] k_prime coefficient de Timoshenko
    \param[out] amo amortissement constant b1
    \param[out] amoB2 amortissement quadratique b3
  */
  void ParameterMultistring
  ::GetParamNote(string note, Real_wp& L, Real_wp& A, Real_wp& rho, Real_wp& T0,
		 Real_wp& E, Real_wp& I, Real_wp& G, Real_wp& k_prime,
		 R3& amo, R3& amoB2)
  {
    // on recupere le numero de la note
    int num = GetNoteNumber(note);
    if (num < 3)
      {
	cout << "Les plans de corde etablis fonctionnent a partir de C1" << endl;
	abort();
      }
    
    L = param_L(num);
    A = param_A(num);
    rho = param_rho(num);
    T0 = param_T0(num);
    E = param_E(num);
    I = param_I(num);
    G = param_G(num);
    k_prime = param_Kprime(num);
    num++;
    
    //DISP(is_string_damped);
    if (is_string_damped)
      {          
	// lois polynomiales pour l'amortissement
	amoB2(0) = 2.78e-11 * num + 1.5274e-09;
	amoB2(1) = 1e-9;
	amoB2(2) = amoB2(0);
        
	amo(0) = 5e-3*num - 0.015;
	if (amo(0) < 0)
	  amo(0) = 0;
	
	amo(1) = 0.5;
	amo(2) = amo(0);
      }
    else 
      {
	amo.Fill(0.0);
	amoB2.Fill(0.0);
      }
    
  }
  
  
  //! calcul du point d'attache pour une note donnee
  void ParameterMultistring::GetParamAttache(string note, Real_wp& x0, Real_wp& y0)
  {
    // on recupere le numero de la note
    int num = GetNoteNumber(note);
    
    GetParamAttache(num, x0, y0);
  }
  
  
  //! calcul du point d'attache pour une note donnee
  void ParameterMultistring::GetParamAttache(int num, Real_wp& x0, Real_wp& y0)
  {
    // on utilise une spline basee sur des points isoles
    x0 = points_attache(num)(0);
    y0 = points_attache(num)(1);
    
  }

}

#define MONTJOIE_FILE_DONNEES_STRING_CXX
#endif
