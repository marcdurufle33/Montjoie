#include "mpreal.h"

namespace mpfr
{
//////////////////////////////////////////////////////////////////////////
// Operators - Assignment
mpreal& mpreal::operator=(const mpreal& v)
{
    if (this != &v)
    {
		mp_prec_t tp = mpfr_get_prec(  mpfr_srcptr());
		mp_prec_t vp = mpfr_get_prec(v.mpfr_srcptr());

		if(tp != vp){
			clear(mpfr_ptr());
			mpfr_init2(mpfr_ptr(), vp);
		}

        mpfr_set(mpfr_ptr(), v.mpfr_srcptr(), mpreal::get_default_rnd());

        MPREAL_MSVC_DEBUGVIEW_CODE;
    }
    return *this;
}

mpreal& mpreal::operator=(const double v)                
{   
#if (MPREAL_DOUBLE_BITS_OVERFLOW > -1)
	if(fits_in_bits(v, MPREAL_DOUBLE_BITS_OVERFLOW))
	{
		mpfr_set_d(mpfr_ptr(),v,mpreal::get_default_rnd());
	}else
		throw conversion_overflow();
#else
	mpfr_set_d(mpfr_ptr(),v,mpreal::get_default_rnd());
#endif

	MPREAL_MSVC_DEBUGVIEW_CODE;
    return *this;
}

mpreal& mpreal::operator=(const char* s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpreal(const char* s,        mp_prec_t prec, int base, mp_rnd_t mode)
    //        mpreal(const std::string& s,mp_prec_t prec, int base, mp_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpfr_t t;

    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));

    if(0 == mpfr_set_str(t, s, 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpfr_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}

mpreal& mpreal::operator=(const std::string& s)
{
    // Use other converters for more precise control on base & precision & rounding:
    //
    //        mpreal(const char* s,        mp_prec_t prec, int base, mp_rnd_t mode)
    //        mpreal(const std::string& s,mp_prec_t prec, int base, mp_rnd_t mode)
    //
    // Here we assume base = 10 and we use precision of target variable.

    mpfr_t t;

    mpfr_init2(t, mpfr_get_prec(mpfr_srcptr()));

    if(0 == mpfr_set_str(t, s.c_str(), 10, mpreal::get_default_rnd()))
    {
        mpfr_set(mpfr_ptr(), t, mpreal::get_default_rnd()); 
        MPREAL_MSVC_DEBUGVIEW_CODE;
    }

    clear(t);
    return *this;
}


#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))

std::string mpreal::toString(const std::string& format) const
{
    char *s = NULL;
    std::string out;

    if( !format.empty() )
    {
        if(!(mpfr_asprintf(&s, format.c_str(), mpfr_srcptr()) < 0))
        {
            out = std::string(s);

            mpfr_free_str(s);
        }
    }

    return out;
}

#endif

std::string mpreal::toString(int n, int b, mp_rnd_t mode) const
{
    // TODO: Add extended format specification (f, e, rounding mode) as it done in output operator
    (void)b;
    (void)mode;

#if (MPFR_VERSION >= MPFR_VERSION_NUM(2,4,0))

    std::ostringstream format;

    int digits = (n >= 0) ? n : 1 + bits2digits(mpfr_get_prec(mpfr_srcptr()));
    
    format << "%." << digits << "RNg";

    return toString(format.str());

#else

    char *s, *ns = NULL; 
    size_t slen, nslen;
    mp_exp_t exp;
    std::string out;

    if(mpfr_inf_p(mp))
    { 
        if(mpfr_sgn(mp)>0) return "+Inf";
        else               return "-Inf";
    }

    if(mpfr_zero_p(mp)) return "0";
    if(mpfr_nan_p(mp))  return "NaN";

    s  = mpfr_get_str(NULL, &exp, b, 0, mp, mode);
    ns = mpfr_get_str(NULL, &exp, b, (std::max)(0,n), mp, mode);

    if(s!=NULL && ns!=NULL)
    {
        slen  = strlen(s);
        nslen = strlen(ns);
        if(nslen<=slen) 
        {
            mpfr_free_str(s);
            s = ns;
            slen = nslen;
        }
        else {
            mpfr_free_str(ns);
        }

        // Make human eye-friendly formatting if possible
        if (exp>0 && static_cast<size_t>(exp)<slen)
        {
            if(s[0]=='-')
            {
                // Remove zeros starting from right end
                char* ptr = s+slen-1;
                while (*ptr=='0' && ptr>s+exp) ptr--; 

                if(ptr==s+exp) out = std::string(s,exp+1);
                else           out = std::string(s,exp+1)+'.'+std::string(s+exp+1,ptr-(s+exp+1)+1);

                //out = string(s,exp+1)+'.'+string(s+exp+1);
            }
            else
            {
                // Remove zeros starting from right end
                char* ptr = s+slen-1;
                while (*ptr=='0' && ptr>s+exp-1) ptr--; 

                if(ptr==s+exp-1) out = std::string(s,exp);
                else             out = std::string(s,exp)+'.'+std::string(s+exp,ptr-(s+exp)+1);

                //out = string(s,exp)+'.'+string(s+exp);
            }

        }else{ // exp<0 || exp>slen
            if(s[0]=='-')
            {
                // Remove zeros starting from right end
                char* ptr = s+slen-1;
                while (*ptr=='0' && ptr>s+1) ptr--; 

                if(ptr==s+1) out = std::string(s,2);
                else         out = std::string(s,2)+'.'+std::string(s+2,ptr-(s+2)+1);

                //out = string(s,2)+'.'+string(s+2);
            }
            else
            {
                // Remove zeros starting from right end
                char* ptr = s+slen-1;
                while (*ptr=='0' && ptr>s) ptr--; 

                if(ptr==s) out = std::string(s,1);
                else       out = std::string(s,1)+'.'+std::string(s+1,ptr-(s+1)+1);

                //out = string(s,1)+'.'+string(s+1);
            }

            // Make final string
            if(--exp)
            {
                if(exp>0) out += "e+"+mpfr::toString<mp_exp_t>(exp,std::dec);
                else       out += "e"+mpfr::toString<mp_exp_t>(exp,std::dec);
            }
        }

        mpfr_free_str(s);
        return out;
    }else{
        return "conversion error!";
    }
#endif
}


//////////////////////////////////////////////////////////////////////////
// I/O
std::ostream& mpreal::output(std::ostream& os) const 
{
    std::ostringstream format;
    const std::ios::fmtflags flags = os.flags();

    format << ((flags & std::ios::showpos) ? "%+" : "%");
    if (os.precision() >= 0)
        format << '.' << os.precision() << "R*"
               << ((flags & std::ios::floatfield) == std::ios::fixed ? 'f' :
                   (flags & std::ios::floatfield) == std::ios::scientific ? 'e' :
                   'g');
    else
        format << "R*e";

    char *s = NULL;
    if(!(mpfr_asprintf(&s, format.str().c_str(),
                        mpfr::mpreal::get_default_rnd(),
                        mpfr_srcptr())
        < 0))
    {
        os << std::string(s);
        mpfr_free_str(s);
    }
    return os;
}

}

