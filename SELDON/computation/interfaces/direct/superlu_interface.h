#ifndef SELDON_FILE_SUPERLU_INTERFACE_H /* allow multiple inclusions */
#define SELDON_FILE_SUPERLU_INTERFACE_H

#ifdef SELDON_WITH_SUPERLU_DIST

namespace superlu
{

#ifdef SELDON_WITH_SUPERLU_DOUBLE

#include "superlu_ddefs.h"
#undef Reduce

#else

#include "superlu_zdefs.h"
#undef Reduce

#endif
}

#else

namespace superlu
{
#ifdef SELDON_WITH_SUPERLU_MT
#include "slu_mt_zdefs.h"
#else
#include "slu_zdefs.h"
#endif
}

#endif

using superlu::int_t;
using superlu::colperm_t;

// since slu_zdefs.h and slu_ddefs.h can not be included
// at the same time, needed functions of slu_ddefs are copied here
// version : SuperLU 5.0, SuperLU_MT 3.0, SuperLU_DIST 4.1
extern "C"
{

namespace superlu
{  

#ifdef SELDON_WITH_SUPERLU_DIST
  
  typedef struct {
    int_t   **Lrowind_bc_ptr; /* size ceil(NSUPERS/Pc)                 */
    double **Lnzval_bc_ptr;  /* size ceil(NSUPERS/Pc)                 */
    double **Linv_bc_ptr;  /* size ceil(NSUPERS/Pc)                 */
    int_t   **Lindval_loc_bc_ptr; /* size ceil(NSUPERS/Pc)  pointers to locations in Lrowind_bc_ptr and Lnzval_bc_ptr */
    int_t   *Unnz; /* number of nonzeros per block column in U*/
	int_t   **Lrowind_bc_2_lsum; /* size ceil(NSUPERS/Pc)  map indices of Lrowind_bc_ptr to indices of lsum  */
    double  **Uinv_bc_ptr;  /* size ceil(NSUPERS/Pc)     	*/
    int_t   **Ufstnz_br_ptr;  /* size ceil(NSUPERS/Pr)                 */
    double  **Unzval_br_ptr;  /* size ceil(NSUPERS/Pr)                 */
        /*-- Data structures used for broadcast and reduction trees. --*/
    BcTree  *LBtree_ptr;       /* size ceil(NSUPERS/Pc)                */
    RdTree  *LRtree_ptr;       /* size ceil(NSUPERS/Pr)                */
    BcTree  *UBtree_ptr;       /* size ceil(NSUPERS/Pc)                */
    RdTree  *URtree_ptr;       /* size ceil(NSUPERS/Pr)			*/
#if 0
    int_t   *Lsub_buf;        /* Buffer for the remote subscripts of L */
    double  *Lval_buf;        /* Buffer for the remote nonzeros of L   */
    int_t   *Usub_buf;        /* Buffer for the remote subscripts of U */
    double  *Uval_buf;        /* Buffer for the remote nonzeros of U   */
#endif
    int_t   *Lsub_buf_2[MAX_LOOKAHEADS];   /* Buffers for the remote subscripts of L*/
    double  *Lval_buf_2[MAX_LOOKAHEADS];   /* Buffers for the remote nonzeros of L  */
    int_t   *Usub_buf_2[MAX_LOOKAHEADS];   /* Buffer for the remote subscripts of U */
    double  *Uval_buf_2[MAX_LOOKAHEADS];   /* Buffer for the remote nonzeros of U   */
    double  *ujrow;           /* used in panel factorization.          */
    int_t   bufmax[NBUFFERS]; /* Maximum buffer size across all MPI ranks:
			       *  0 : maximum size of Lsub_buf[]
			       *  1 : maximum size of Lval_buf[]
			       *  2 : maximum size of Usub_buf[]
			       *  3 : maximum size of Uval_buf[]
			       *  4 : maximum size of tempv[LDA]
			       */

    /*-- Record communication schedule for factorization. --*/
    int   *ToRecv;          /* Recv from no one (0), left (1), and up (2).*/
    int   *ToSendD;         /* Whether need to send down block row.       */
    int   **ToSendR;        /* List of processes to send right block col. */

    /*-- Record communication schedule for forward/back solves. --*/
    int_t   *fmod;            /* Modification count for L-solve            */
    int_t   **fsendx_plist;   /* Column process list to send down Xk       */
    int_t   *frecv;           /* Modifications to be recv'd in proc row    */
    int_t   nfrecvx;          /* Number of Xk I will receive in L-solve    */
    int_t   nfsendx;          /* Number of Xk I will send in L-solve       */
    int_t   *bmod;            /* Modification count for U-solve            */
    int_t   **bsendx_plist;   /* Column process list to send down Xk       */
    int_t   *brecv;           /* Modifications to be recv'd in proc row    */
    int_t   nbrecvx;          /* Number of Xk I will receive in U-solve    */
    int_t   nbsendx;          /* Number of Xk I will send in U-solve       */
    int_t   *mod_bit;         /* Flag contribution from each row blocks    */

    /*-- Auxiliary arrays used for forward/back solves. --*/
    int_t   *ilsum;           /* Starting position of each supernode in lsum
				 (local)  */
    int_t   ldalsum;          /* LDA of lsum (local) */
    int_t   SolveMsgSent;     /* Number of actual messages sent in LU-solve */
    int_t   SolveMsgVol;      /* Volume of messages sent in the solve phase */


    /*********************/
    /* The following variables are used in the hybrid solver */

    /*-- Counts to be used in U^{-T} triangular solve. -- */
    int_t UT_SOLVE;
    int_t L_SOLVE;
    int_t FRECV;
    int_t ut_ldalsum;        /* LDA of lsum (local) */
    int_t *ut_ilsum;         /* ilsum in column-wise                        */
    int_t *utmod;            /* Modification count for Ut-solve.            */
    int_t **ut_sendx_plist;  /* Row process list to send down Xk            */
    int_t *utrecv;           /* Modifications to be recev'd in proc column. */
    int_t n_utsendx;         /* Number of Xk I will receive                 */
    int_t n_utrecvx;         /* Number of Xk I will send                    */
    int_t n_utrecvmod;
    int_t nroot;
    int_t *ut_modbit;
    int_t *Urbs;
    Ucb_indptr_t **Ucb_indptr;/* Vertical linked list pointing to Uindex[] */
    int_t  **Ucb_valptr;      /* Vertical linked list pointing to Unzval[] */

    /* some additional counters for L solve */
    int_t n;
    int_t nleaf;
    int_t nfrecvmod;
    int_t inv; /* whether the diagonal block is inverted*/
} dLocalLU_t;
  
  typedef struct {
    int_t *etree;
    Glu_persist_t *Glu_persist;
    dLocalLU_t *Llu;
    char dt;
} dLUstruct_t;

  /*-- Data structure for communication during matrix-vector multiplication. */
typedef struct {
    int_t *extern_start;
    int_t *ind_tosend;    /* X indeices to be sent to other processes */
    int_t *ind_torecv;    /* X indeices to be received from other processes */
    int_t *ptr_ind_tosend;/* Printers to ind_tosend[] (Size procs)
			     (also point to val_torecv) */
    int_t *ptr_ind_torecv;/* Printers to ind_torecv[] (Size procs)
			     (also point to val_tosend) */
    int   *SendCounts;    /* Numbers of X indices to be sent
			     (also numbers of X values to be received) */
    int   *RecvCounts;    /* Numbers of X indices to be received
			     (also numbers of X values to be sent) */
    double *val_tosend;   /* X values to be sent to other processes */
    double *val_torecv;   /* X values to be received from other processes */
    int_t TotalIndSend;   /* Total number of indices to be sent
			     (also total number of values to be received) */
    int_t TotalValSend;   /* Total number of values to be sent.
			     (also total number of indices to be received) */
} pdgsmv_comm_t;

  typedef struct {
    int_t *row_to_proc;
    int_t *inv_perm_c;
    int_t num_diag_procs, *diag_procs, *diag_len;
    pdgsmv_comm_t *gsmv_comm; /* communication metadata for SpMV,
         	       		      required by IterRefine.          */
    pxgstrs_comm_t *gstrs_comm;  /* communication metadata for SpTRSV. */
    int_t *A_colind_gsmv; /* After pdgsmv_init(), the global column
                             indices of A are translated into the relative
                             positions in the gathered x-vector.
                             This is re-used in repeated calls to pdgsmv() */
    int_t *xrow_to_proc; /* used by PDSLin */
} dSOLVEstruct_t;

  typedef struct {
    DiagScale_t DiagScale;
    double *R;
    double *C; 
    int_t  *perm_r;
    int_t  *perm_c;
} dScalePermstruct_t;

  int_t dQuerySpace_dist(int_t, dLUstruct_t *, gridinfo_t *,
                         SuperLUStat_t *, superlu_dist_mem_usage_t*);
  
  void dCreate_CompRowLoc_Matrix_dist(SuperMatrix *, int_t, int_t, int_t, int_t,
                                      int_t, double *, int_t *, int_t *,
                                      Stype_t, Dtype_t, Mtype_t);
  
  void  pdgssvx(superlu_dist_options_t *, SuperMatrix *, 
                dScalePermstruct_t *, double *,
                int, int, gridinfo_t *, dLUstruct_t *,
                dSOLVEstruct_t *, double *, SuperLUStat_t *, int *);

  void dLUstructFree(dLUstruct_t *);
  void   dScalePermstructFree(dScalePermstruct_t *);
  void dDestroy_LU(int_t, gridinfo_t *, dLUstruct_t *); 
#else
  
  void
  dCreate_CompCol_Matrix(SuperMatrix *, int_t, int_t, int_t, double *,
                         int_t *, int_t *, Stype_t, Dtype_t, Mtype_t);
  
  void
  dCreate_Dense_Matrix(SuperMatrix *, int_t, int_t, double *, int_t,
                       Stype_t, Dtype_t, Mtype_t);
    
#ifdef SELDON_WITH_SUPERLU_MT
  int_t  superlu_dQuerySpace (int_t, SuperMatrix *, SuperMatrix *, int_t, 
                              superlu_memusage_t *);
  
  void pdgstrf (superlumt_options_t *, SuperMatrix *, int_t *, 
                SuperMatrix *, SuperMatrix *, Gstat_t *, int_t *);

  void pdgstrf_init (int_t, fact_t, trans_t, yes_no_t, int_t, int_t, double, yes_no_t, double,
                     int_t *, int_t *, void *, int_t, SuperMatrix *,
                     SuperMatrix *, superlumt_options_t *, Gstat_t *);
  
  void dgstrs (trans_t, SuperMatrix *, SuperMatrix*, 
               int_t*, int_t*, SuperMatrix*, Gstat_t *, int_t *);  
#else
  void    dgstrf (superlu_options_t*, SuperMatrix*,
                  int, int, int*, void *, int, int *, int *, 
                  SuperMatrix *, SuperMatrix *, GlobalLU_t *,
                  SuperLUStat_t*, int *);
  
  void    dgstrs (trans_t, SuperMatrix *, SuperMatrix *, int *, int *,
                  SuperMatrix *, SuperLUStat_t*, int *);
  
  int     dQuerySpace (SuperMatrix *, SuperMatrix *, mem_usage_t *);
#endif

#endif
}
  
}

#endif /* __SUPERLU_INTERFACE */
