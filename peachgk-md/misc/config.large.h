/*******************************
 *  config.h Ver.3.0           *
 *      for peachgk_md.f       *
 *            by G.Kikugawa    *
 *******************************/
// Time-stamp: <>

/***** configuration of the dimension of arrays *****/

#define _PGKINIVER_ "5.8"
#define _PGKINIDATE_ "18.10.03"

#define _MAXNSTAGE_ 10

#define _MAXNPROC_ 1000

#define _MAXNINTTYP_ 20

#define _MAXNATOM_ 100000
#define _MAXNMOL_  500
#define _MAXNBOND_ 100000
#define _MAXNANGL_ 100000
#define _MAXNTORS_ 100000
#define _MAXNCONST_ 100000
#define _MAXNLIST_ 1500
#define _MAXMORLIST_ 1500
#define _MAXSHLIST_ 500
#define _MAXRFHLIST_ 500
#define _MAXDOULIST_ 500
#define _MAXRPVWLIST_ 1500
#define _MAXCSTMNBLIST_ 500
#define _MAXNPOLYTYP_ 20
#define _MAXNMATYP_ 20

#define _MAXNVDWTYP_ 500
#define _MAXNATMTYP_ 500
#define _MAXNBONDTYP_ 200
#define _MAXNANGLTYP_ 200
#define _MAXNTORSTYP_ 100
#define _MAXNMORTYP_ 100
#define _MAXNSHTYP_ 5
#define _MAXNRFHTYP_ 15
#define _MAXNDOUTYP_ 10
#define _MAXNRPVWTYP_ 10
#define _MAXNCSTMNBTYP_ 100
#define _MAXNCONSTTYP_ 100

#define _MAXWAVE_ 4631              // =(kmax*2+1)**3/2

#define _MAXNHC_ 10                 // max num. of Nose-Hoover chain
#define _MAXNYOSH_ 5                // max num. of Yoshida-Suzuki method

#define _MAXCELLX_ 20               // max num. of cell-index cell (x)
#define _MAXCELLY_ 20               // max num. of cell-index cell (y)
#define _MAXCELLZ_ 25               // max num. of cell-index cell (z)
#define _MAX_ATOM_CELL_ 1000        // maximum number of atoms for each cell

#define _MAX_INTERPOL_ 10000        // max num. of spline interpolation points
#define _MAX_SPL_ORDER_ 3           // maximum order of spline function

#define _MAXNTCREGION_ 20           // maximum num. of regions for temp. control

#define _MAXNWORD_ 20

#define _MAXNHFREGION_ 30

#define _MAXNINTPBIAS_ 500          // maximum num. of sampling according to
                                    // reaction coordinate with bias potential

#define _MAXNVELREGION_ 500         // max num. of slabs for streaming velocity

/*** macro for spme method ***/
#define _MAXFFT1_ 200
#define _MAXFFT2_ 200
#define _MAXFFT3_ 250
#define _MAX_GRID_ 250               // this macro should be max(MAXFFT1,2,3)

#define _MAXORDER_ 10               // max order of B-spline

/*** macro for MDGRAPE-3 library ***/
#define _MAX_TYPE_MDG_ 32

#define _MAX_CELL_X_VDW_ 5           // max cell index for x-direction (vdW)
#define _MAX_CELL_Y_VDW_ 5           // max cell index for y-direction (vdW)
#define _MAX_CELL_Z_VDW_ 5           // max cell index for z-direction (vdW)

#define _MAX_ATOM_CELL_VDW_ 1000     // max number of atoms per cell (vdW)

#define _MAX_CELL_ALLATOM_VDW_ 200000 // max number of all atoms of j particle

#define _MAX_CELL_X_ELE_ 5           // max cell index for x-direction (ele)
#define _MAX_CELL_Y_ELE_ 5           // max cell index for y-direction (ele)
#define _MAX_CELL_Z_ELE_ 5           // max cell index for z-direction (ele)

#define _MAX_ATOM_CELL_ELE_ 1000     // max number of atoms per cell (ele)

#define _MAX_CELL_ALLATOM_ELE_ 200000 // max number of all atoms of j particle
