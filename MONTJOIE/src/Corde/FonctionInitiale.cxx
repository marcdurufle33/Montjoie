namespace Montjoie{

  Real_wp FonctionInitiale(Real_wp x, Real_wp n , Real_wp x1, Real_wp U,Real_wp V,int k, int dimension, int fonctioninittype){
	Real_wp res=0;
	//if(k==0 ){
			//|| k==0
		//	if ( (x>(V/2-V/5)) && (x<(V/2+V/5))){
//			res=U*exp(-(x-V/2)*(x-V/2)/((V/n)*(V/n)));
//			}
	         switch (fonctioninittype){
			case 0: {
// Un sinus sur u :
			  if(k==0){res=U*sin(n*pi_wp*x/V);}
			  else {res=0.0;}
			  break;
			}
// Une gaussienne sur u :
		        case 1:{
				res=U*0.1*exp(-((x-0.45)/0.1)*((x-0.45)/0.1));
				//cout<<"je fais une gaussienne" <<endl;	
				break;
			}
		        case 2: {
// Une fonction qui excite u et v
			  switch (k) {
			  case 0:
			    {
			      res=U*exp(-((x-0.45)/0.1)*((x-0.45)/0.1))*cos(pi_wp/4);
			      //cout<<"gaussienne sur u "<<"k = "<<k << endl;
			      break;
			    }
			  case 1:
			    {
			      res=U*exp(-((x-0.45)/0.1)*((x-0.45)/0.1))*sin(pi_wp/4);
			      //cout<<"gaussienne sur v "<<"k = "<<k << endl;
			      break;
			    }		
			  }
			  break;
			}
		 }
				
  
		//if(x==V){
//			res=0;
//		}
	
	//else if(k==dimension-1){
//		res=U*2*pi_wp/V*cos(2*pi_wp*x/V)*(V-x)*x;
//	}
	//DISP(x1);
	//Real_wp x2 = x1+Delta;
	
		// if (x<x1){ // partie gauche : tout nul
		// 		if(k==0){
		// 		// // **** le sinus ****
		// 		// 				res = U*sin(2*pi_wp*x);
		// 		// 		// *************************	
		// 		}
		// 		else{
		// 		res=0.0;
		// 	}
		// 	}
		// 	if(x>=x1 & x<x2){
		// 	//	if(k==dimension-1){ // on ne remplit que pour v
		// 			if(k==0){ // on ne remplit que pour u
		// 		
		// 		// **** l'exponentielle ****
		// 		Real_wp phi = 2.0*((x-x1)/Delta -1.0/2.0);
		// 		res = U*exp(pow(phi,2)/(pow(phi,2)-1));
		// 		// *************************
		// 		// **** le sinus ****
		// 		//		res = U*sin(2*pi_wp*x);
		// 		// *************************
		// 		}
		// 		else{
		// 		res=0.0;
		// 	}
		// 	}
		// 	
		// 	if(x>x2){// partie droite : tout nul
		// 		if(k==0){
		// 		// **** le sinus ****
		// 			//	res = U*sin(2*pi_wp*x);
		// 		// *************************	
		// 		}
		// 		else{
		// 		res=0.0;
		// 	}
		// 	}
		
	return(res);
  }
}
