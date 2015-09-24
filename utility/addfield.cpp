//////////////////////////////////////////////////////////////////////////////
// addfield.C: process semtex/NEKTON-type field files, computing and
// adding vorticity and divergence, rate of strain magnitude, velocity
// gradient discriminant, etc.
//
// Copyright (c) 1998 <--> $Date: 2015/08/13 06:16:20 $, 
//   Hugh Blackburn, Murray Rudman
//
// NB: the input field file is assumed to contain ONLY velocity and
// pressure data.
//
// Usage:
// -----
// addfield [options] -s session session.fld
//   options:
//   -h        ... print this message
//   -q        ... add kinetic energy per unit mass 0.5(u.u) (default)
//   -d        ... add divergence div(u)
//   -v        ... add vorticity w=curl(u)
//   -e        ... add enstrophy 0.5(w.w)
//   -H        ... add helicity 0.5(u.w) (3D only)
//   -L        ... add divergence of Lamb vector div(uxw) (3D only)
//   -g        ... add strain rate magnitude sqrt(2SijSji)
//   -D        ... add discriminant of velocity gradient tensor
//                 NB: divergence is assumed to be zero. (3D only)
//   -J        ... add vortex core measure of Jeong & Hussain. (3D only)
//   -a        ... add all fields derived from velocity (above)
//   -f <func> ... add a computed function <func> of x, y, z, t, etc.
//
// Reserved field names used/assumed:
// ---------------------------------
//
// u -- x velocity component (cylindrical: axial)
// v -- y velocity component (cylindrical: radial)
// w -- z velocity component (cylindrical: azimuthal)
// p -- pressure/density
// c -- scalar
//
// The following are reserved names, not used by addfield.
// ------------------------------------------------------
//
// A -- uu covariance
// B -- uv covariance
// C -- vv covariance
// D -- uw covariance
// E -- vw covariance
// F -- ww covariance
// 
// Computed variables:
// ------------------
//
// d -- divergence
// e -- enstrophy 0.5*(r^2 + s^2 + t^2) = 0.5 (omega . omega)
// f -- a computed function of spatial variables
// g -- strain rate magnitude sqrt(2SijSij)
// q -- kinetic energy per unit mass 0.5*(u^2 + v^2 + w^2) = 0.5 (u . u)
// r -- x component vorticity
// s -- y component vorticity
// t -- z component vorticity
// H -- helicity  0.5*(u*r + v*s + w*t) = 0.5 (u . omega) 3D only.
// J -- vortex core identification measure, see [2]. 3D only.
// D -- discriminant of velocity gradient tensor, see [1]. 3D only.
// L -- divergence of Lamb vector = div (omega x u), see [3] 3D only.
//
// NB: product terms -- such as are used to calculate enstrophy,
// helicity, the invariants and discriminant of the velocity gradient
// tensor, and the strain rate magnitude, all computed in physical
// space -- are not dealiased.  Therefore it is advisable to project
// the original field to a greater number of planes (3/2 rule) before
// these terms are calculated, otherwise the products can be quite
// different from what would be expected (especially if N_Z is small,
// say 4). If this is done you need to edit a matching session file
// with the appropriate value of N_Z.
//
// References
//-----------
//
// [1] Chong, Perry & Cantwell (1990) A general classification of
// three-dimensional flow fields, PF(A) 2:765--777; see also Blackburn
// et al. (1996) JFM 310:269--292
//
// [2] Jeong & Hussain (1995) On the identification of a vortex, JFM
// 285:69--94
//
// [3] Hamman, Klewicki & Kirby (2008) On the Lamb vector divergence
// in Navier--Stokes flows JFM 610:261--284
//
// --
// This file is part of Semtex.
// 
// Semtex is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
// 
// Semtex is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Semtex (see the file COPYING); if not, write to the Free
// Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA
// 02110-1301 USA
//////////////////////////////////////////////////////////////////////////////

static char RCS[] = "$Id: addfield.cpp,v 8.2 2015/08/13 06:16:20 hmb Exp $";

#include <sem.h>
#include <tensorcalcs.h>

#define FLDS_MAX 64 // -- More than we'll ever want.
#define FLAG_MAX 10 // -- NB: FLAG_MAX should tally with the following enum:
enum {
  ENERGY      ,     // -- NB: the placing of ENERGY and FUNCTION in the first
  FUNCTION    ,	    //    two positions is significant: don't break this.
  DIVERGENCE  ,
  ENSTROPHY   ,
  DISCRIMINANT,
  HELICITY    ,
  DIVLAMB  ,
  STRAINRATE  ,
  VORTICITY   ,
  VORTEXCORE
};

static char  prog[] = "addfield";
static void  getargs  (int, char**, char*&, char*&, char*&, bool[]);
static bool  getDump  (Domain*, ifstream&);
static void  putDump  (Domain*, vector<AuxField*>&, int_t, ostream&);


int main (int    argc,
	  char** argv)
// ---------------------------------------------------------------------------
// Driver.
// ---------------------------------------------------------------------------
{
  Geometry::CoordSys         system;
  char                       *session, *dump, *func, fields[StrMax];
  int_t                      i, j, k, p, q;
  int_t                      np, nz, nel, allocSize, NCOM, NDIM;
  int_t                      iAdd = 0;
  bool                       add[FLAG_MAX], need[FLAG_MAX], gradient;
  ifstream                   file;
  FEML*                      F;
  Mesh*                      M;
  BCmgr*                     B;
  Domain*                    D;
  vector<Element*>           elmt;
  AuxField                   *Ens, *Hel, *Div, *Disc, *Strain;
  AuxField                   *Func, *Vtx, *DivL, *Nrg, *work;
  vector<AuxField*>          velocity, vorticity, lamb, addField(FLDS_MAX);
  vector<vector<AuxField*> > Vij;     // -- Usually computed, for internal use.
  vector<vector<real_t*> >   VijData; // -- For pointwise access in Vij.

  vector<real_t*> VorData; // -- Ditto in vorticity.
  vector<real_t*> LamData; // -- Ditto in Lamb vector.

  real_t          *DisData, *DivData, *StrData, *VtxData, *HelData, *EnsData;
  real_t          vel[3], vort[3], tensor[9];

  Femlib::initialize (&argc, &argv);
  for (i = 0; i < FLAG_MAX; i++) add [i] = need [i] = false;

  // -- Read command line.

  getargs (argc, argv, session, dump, func, add);

  file.open (dump);
  if (!file) message (prog, "can't open input field file", ERROR);

  // -- Set up domain.

  F      = new FEML (session);
  M      = new Mesh (F);
  nel    = M -> nEl ();  
  np     =  Femlib::ivalue ("N_P");
  nz     =  Femlib::ivalue ("N_Z");
  system = (Femlib::ivalue ("CYLINDRICAL") ) ?
                       Geometry::Cylindrical : Geometry::Cartesian;
  Geometry::set (np, nz, nel, system);

  allocSize = Geometry::nTotal();

  elmt.resize (nel);
  for (i = 0; i < nel; i++) elmt[i] = new Element (i, np, M);

  B = new BCmgr  (F, elmt);
  D = new Domain (F, elmt, B);

  if      (strstr (D -> field, "uvw")) NCOM = 3;
  else if (strstr (D -> field, "uv"))  NCOM = 2;
  else message (prog, "lacking velocity components: is session valid?", ERROR);
  NDIM = Geometry::nDim();

  velocity.resize   (NCOM);
  vorticity.resize ((NDIM == 2) ? 1 : 3);
  for (i = 0; i < NCOM; i++) velocity[i] = D -> u[i];

  // -- From the requested fields, flag dependencies.

  // -- First, only allow the "coherent structures" measures for flows
  //    that are 3D.

  if  (NDIM == 2)
    add[HELICITY]=add[DISCRIMINANT]=add[VORTEXCORE]=add[DIVLAMB] = false;

  for (p = 0, i = 0; i < FLAG_MAX; i++) p += (add[i]) ? 1 : 0;
  if  (p == 0) message (prog, "nothing to be done", ERROR);

  // -- Check if we just have the (first two) cases not requiring derivatives.

  for (p = 0, i = 0; i < FLAG_MAX; i++) p += (add[i]) ? (i + 1) : 0;
  if (p <= 2) gradient = false; else gradient = true;

  for (i = 0; i < FLAG_MAX; i++) need[i] = add[i];

  // -- If any gradients need to be computed, we make all of them.
  //    Vij = du_j/dx_i.  So j labels velocity component and i labels
  //    spatial dimension.

  // -- Despite the fact that it's overkill we will take the
  //    simple-minded approach and always make Vij as a 3x3 tensor.
  //    That's a wasteful for anything which is 2D but presumably for
  //    those cases the memory can be afforded.

  // -- The velocity gradient tensor is initially just a naive matrix
  //    of derivatives.  In cylindrical coordinates we need to make
  //    modifications.
  
  if (gradient) {
    Vij    .resize (3);
    VijData.resize (3);
    for (i = 0; i < 3; i++) {
      Vij    [i].resize (3);
      VijData[i].resize (3);
      for (j = 0; j < 3; j++) {
	VijData[i][j] = new real_t [allocSize];
	Vij    [i][j] = new AuxField (VijData[i][j], nz, elmt);
	*Vij   [i][j] = 0.0;
      }
    }
  }
  
  // -- Fields without dependants.

  if (need[ENERGY]) {
    Nrg = new AuxField (new real_t[allocSize], nz, elmt, 'q');
    addField[iAdd++] = Nrg;
  }

  if (need[FUNCTION]) {
    Func = new AuxField (new real_t[allocSize], nz, elmt, 'f');
    addField[iAdd++] = Func;
  }

  if (need[DIVERGENCE]) {
    DivData = new real_t [allocSize];
    *(Div  =  new AuxField (DivData, nz, elmt, 'd')) = 0.0;
    addField[iAdd++] = Div;
  }

  if (need[DISCRIMINANT]) {
    DisData = new real_t [allocSize];
    *(Disc = new AuxField (DisData, nz, elmt, 'D')) = 0.0;
    addField[iAdd++] = Disc;
  }

  if (need[STRAINRATE]) {
    StrData = new real_t [allocSize];
    Strain  = new AuxField (StrData, nz, elmt, 'g');
    addField[iAdd++] = Strain;
  }

  if (need[VORTEXCORE]) {
    VtxData = new real_t [allocSize];
    Vtx = new AuxField (VtxData, nz, elmt, 'J');
    addField[iAdd++] = Vtx;
  }

  if (need[VORTICITY])
    if (NDIM == 2) {
      vorticity.resize (1);
      VorData  .resize (1);
      VorData[0]       = new real_t [allocSize];
      vorticity[0]     = new AuxField (VorData[0], nz, elmt, 't');
      addField[iAdd++] = vorticity[0];
    } else {
      vorticity.resize (3);
      VorData  .resize (3);
      for (i = 0; i < 3; i++) {
	VorData[i]       = new real_t [allocSize];
	vorticity[i]     = new AuxField (VorData[i], nz, elmt, 'r' + i);
	addField[iAdd++] = vorticity[i];
      }
    }

  if (add[DIVLAMB]) { 		// -- Know also NDIM == 3.
    lamb   .resize (3);
    LamData.resize (3);
    for (i = 0; i < 3; i++) {
      LamData[i] = new real_t [allocSize];
      lamb[i]    = new AuxField (LamData[i], nz, elmt, 'L');
    }
    addField[iAdd++] = lamb[0];	// -- Where divergence will get stored.
  }

  if (need[ENSTROPHY]) {
    EnsData = new real_t[allocSize];
    Ens = new AuxField (EnsData, nz, elmt, 'e');
    if (add[ENSTROPHY]) addField[iAdd++] = Ens;
  }
  
  if (need[HELICITY]) {
    HelData = new real_t [allocSize];
    Hel = new AuxField (HelData, nz, elmt, 'H');
    if (add[HELICITY]) addField[iAdd++] = Hel;
  }

  // -- Cycle through field dump, first (if required) making the only
  //    two things which just need the velocity and no gradients, then
  //    if needed computing the full VG tensor and work forward from
  //    there.  The order of computation is determined by
  //    dependencies.  Then write requested output, listed in
  //    addField.
  
  while (getDump (D, file)) {
        
    if (need[FUNCTION]) *Func = func;

    if (need[ENERGY]) ((*Nrg) . innerProduct (velocity, velocity)) *= 0.5;

    if (gradient) {		// -- All other things.

      // -- First make all VG components.

      for (i = 0; i < NDIM ; i++)
	for (j = 0; j < NCOM ; j++) {
	  *Vij[i][j] = *velocity[j];
	  if (i == 2) Vij[i][j] -> transform (FORWARD);
	  Vij[i][j] -> gradient (i);
	  if (i == 2) Vij[i][j] -> transform (INVERSE);
	}

      if (Geometry::cylindrical()) {
	work = new AuxField (new real_t[allocSize],  nz, elmt);
	if (NDIM == 3) for (j = 0; j < NCOM; j++) Vij[2][j] -> divY();
	(*work = *velocity[1]) . divY(); *Vij[2][2] += *work;
#if 1
	if (NCOM == 3) { (*work = *velocity[2]) . divY(); *Vij[1][2] += *work; }
#else
	if (NCOM == 3) { (*work = *velocity[2]) . divY(); *Vij[1][2] -= *work; }
#endif
      }

#if 1
      // -- Loop over every point in the mesh and compute everything
      //    from Vij.  Quite likely this could be made more efficient
      //    but for now simplicity is the aim.

      for (i = 0; i < allocSize; i++) {

	for (k = 0, p = 0; p < 3; p++) {
	  for (q = 0; q < 3; q++, k++)
	    tensor [k] = VijData [p][q][i];
	}

	// -- These operations produce a simple scalar result from Vij.

	if (need[DIVERGENCE])   DivData[i] = tensor3::trace      (tensor);
	if (need[ENSTROPHY])    EnsData[i] = tensor3::enstrophy  (tensor);
	if (need[DISCRIMINANT]) DisData[i] = tensor3::discrimi   (tensor);
	if (need[STRAINRATE])   StrData[i] = tensor3::strainrate (tensor);
	if (need[VORTEXCORE])   VtxData[i] = tensor3::lambda2    (tensor);

	// -- Vorticity could be considered scalar in 2D.

	if (need[VORTICITY]) {
	  tensor3::vorticity (tensor, vort);
	  if (NDIM == 2) 
	    VorData[0][i] = vort[2];
	  else { 
	    VorData[0][i] = vort[0]; 
	    VorData[1][i] = vort[1]; 
	    VorData[2][i] = vort[2];
	  }
	}

	if (!(need[HELICITY] || need[DIVLAMB])) continue;

	// -- Last two measures need velocity too, only made for NDIM = 3.

	vel[0] = velocity[0] -> data()[i];
	vel[1] = velocity[1] -> data()[i];
	vel[2] = velocity[2] -> data()[i];

	if (need[HELICITY]) HelData[i] = tensor3::helicity (tensor, vel);


	if (need[DIVLAMB]) {
	  tensor3::lambvector (tensor, vel, vort); // -- vort is a dummy.
	  LamData[0][i] = vort[0]; 
	  LamData[1][i] = vort[1]; 
	  LamData[2][i] = vort[2];
	}
      }

      // -- For this case, we still need to compute a divergence:
      
      if (need[DIVLAMB]) {
	lamb[0] -> gradient (0);
	(*lamb[2]) . transform (FORWARD) . gradient (2) . transform(INVERSE);
	if (Geometry::cylindrical()) lamb[2] -> divY();
	*lamb[0] += *lamb[2];
	if (Geometry::cylindrical()) *lamb[0] += (*lamb[2] = *lamb[1]) . divY();
	*lamb[0] += (*lamb[1]) . gradient (1);
      }
    }
#else

      if (need[DISCRIMINANT]) {

	// -- 2nd invariant (Q from Chong et al.).
	
	InvQ -> times      (*Vij[0][0], *Vij[1][1]);
	InvQ -> timesMinus (*Vij[0][1], *Vij[1][0]);
	InvQ -> timesPlus  (*Vij[0][0], *Vij[2][2]);

	InvQ -> timesMinus (*Vij[0][2], *Vij[2][0]);
	InvQ -> timesPlus  (*Vij[1][1], *Vij[2][2]);
	InvQ -> timesMinus (*Vij[1][2], *Vij[2][1]);

	// -- 3rd invariant: determinant of Vij (R from Chong et al.).

	work -> times      (*Vij[1][1], *Vij[2][2]);
	work -> timesMinus (*Vij[2][1], *Vij[1][2]);
	InvR -> times      (*work,      *Vij[0][0]);

	work -> times      (*Vij[1][2], *Vij[2][0]);
	work -> timesMinus (*Vij[2][2], *Vij[1][0]);
	InvR -> timesPlus  (*work,      *Vij[0][1]);

	work -> times      (*Vij[2][1], *Vij[1][0]);
	work -> timesMinus (*Vij[1][1], *Vij[2][0]);
	InvR -> timesPlus  (*work,      *Vij[0][2]);

	// -- Discriminant L of Vij.
	//    NB: DIVERGENCE (P from Chong et al.) ASSUMED = 0.

	work -> times (*InvQ, *InvQ);
	Disc -> times (*work, *InvQ);
	work -> times (*InvR, *InvR);
	*work *= 6.75;
	*Disc += *work;
      }

      if (need[DIVERGENCE])
	for (i = 0; i < nComponent; i++) *Div += *Vij[i][i];

    
      if (need[VORTICITY]) {
	if (nComponent == 2) {
	  *vorticity[0]  = *Vij[1][0];
	  *vorticity[0] -= *Vij[0][1];
	} else {
	  *vorticity[0]  = *Vij[2][1];
	  *vorticity[0] -= *Vij[1][2];
	  *vorticity[1]  = *Vij[0][2];
	  *vorticity[1] -= *Vij[2][0];
	  *vorticity[2]  = *Vij[1][0];
	  *vorticity[2] -= *Vij[0][1];
	}
      }

      if (need[DIVLAMB]) {
	if (nComponent == 2) {
	  *lamb[0] = 0.0;
	  lamb[0] -> timesMinus (*D -> u[1], *vorticity[0]);
	  lamb[1] -> times      (*D -> u[0], *vorticity[0]);
	  *DivL  = (*work = *lamb[0]) . gradient (0);
	  *DivL += (*work = *lamb[1]) . gradient (1);
	  if (Geometry::cylindrical())
	    *DivL += (*work = *lamb[1]) . divY();
	} else {
	  lamb[0] -> times      (*D -> u[2], *vorticity[1]);
	  lamb[0] -> timesMinus (*D -> u[1], *vorticity[2]);
	  lamb[1] -> times      (*D -> u[0], *vorticity[2]);
	  lamb[1] -> timesMinus (*D -> u[2], *vorticity[0]);
	  lamb[2] -> times      (*D -> u[1], *vorticity[0]);
	  lamb[2] -> timesMinus (*D -> u[0], *vorticity[1]);
	  *DivL  = (*work = *lamb[0]) . gradient (0);
	  *DivL += (*work = *lamb[1]) . gradient (1);
	  (*work = *lamb[2]).transform(FORWARD).gradient(2).transform(INVERSE);
	  if (Geometry::cylindrical()) {
	    *DivL += work -> divY();
	    *DivL += (*work = *lamb[1]) . divY();
	  } else
	    *DivL += *work;
	}
      }
    
      if (need[ENSTROPHY])
	Ens -> innerProduct (vorticity, vorticity) *= 0.5;

      if (need[HELICITY])
	Hel -> innerProduct (vorticity, velocity)  *= 0.5;

      if (need[STRAINTENSOR]) {
	for (i = 0; i < nComponent; i++)
	  for (j = i; j < nComponent; j++) {
	    *Sij[i][j]  = *Vij[i][j];
	    *Sij[i][j] += *Vij[j][i];
	    *Sij[i][j] *= 0.5;
	  }
      }

      if (need[STRAINRATE]) {
	*Strain = 0.0;
	for (i = 0; i < nComponent; i++)
	  for (j = 0; j < nComponent; j++)
	    Strain -> timesPlus (*Sij[i][j], *Sij[j][i]);
	(*Strain *= 2.0) . sqroot();
      }

      if (need[VORTEXCORE])	// -- Only done for 3-component fields.
	for (i = 0; i < allocSize; i++) {
	  for (k = 0, p = 0; p < 3; p++)
	    for (q = 0; q < 3; q++, k++)
	      tensor [k] = VijData[p][q][i];
	  VtxData[i] = lambda2 (tensor);
	}
    }
#endif
    // -- Finally, add mass-projection smoothing on everything.

    for (i = 0; i < iAdd; i++) D -> u[0] -> smooth (addField[i]);

    putDump (D, addField, iAdd, cout);
  }
  
  file.close();
  Femlib::finalize();
  return EXIT_SUCCESS;
}


static void getargs (int    argc   ,
		     char** argv   ,
		     char*& session,
		     char*& dump   ,
		     char*& func   ,
		     bool*  flag   )
// ---------------------------------------------------------------------------
// Deal with command-line arguments.
// ---------------------------------------------------------------------------
{
  char usage[] =
    "Usage: %s [options] -s session dump.fld\n"
    "options:\n"
    "  -h        ... print this message \n"
    "  -q        ... add kinetic energy per unit mass 0.5(u.u) (default)\n"
    "  -d        ... add divergence div(u)\n"
    "  -v        ... add vorticity w=curl(u)\n"
    "  -e        ... add enstrophy 0.5(w.w)\n"
    "  -H        ... add helicity 0.5(u.w) (3D only)\n"
    "  -L        ... add divergence of Lamb vector, div(uxw)\n"
    "  -g        ... add strain rate magnitude sqrt(2SijSji)\n"
    "  -D        ... add discriminant of velocity gradient tensor\n"
    "                NB: divergence is assumed to be zero. (3D only)\n"
    "  -J        ... add vortex core measure of Jeong & Hussain (3D only)\n"
    "  -a        ... add all fields derived from velocity (above)\n"
    "  -f <func> ... add a computed function <func> of x, y, z, t, etc.\n";
              
  int_t i, sum = 0;
  char  buf[StrMax];
 
  while (--argc  && **++argv == '-')
    switch (*++argv[0]) {
    case 'h':
      sprintf (buf, usage, prog);
      cout << buf;
      exit (EXIT_SUCCESS);
      break;
    case 's':
      if (*++argv[0]) session = *argv; else { --argc; session = *++argv; }
      break;
    case 'v': flag[VORTICITY]    = true; break;
    case 'e': flag[ENSTROPHY]    = true; break;
    case 'H': flag[HELICITY]     = true; break;
    case 'd': flag[DIVERGENCE]   = true; break;
    case 'g': flag[STRAINRATE]   = true; break;
    case 'D': flag[DISCRIMINANT] = true; break;
    case 'J': flag[VORTEXCORE]   = true; break;
    case 'L': flag[DIVLAMB]      = true; break;
    case 'q': flag[ENERGY]       = true; break;
    case 'a': flag[0]=true; for (i = 2; i < FLAG_MAX ; i++) flag[i]=true; break;
    case 'f':
      if (*++argv[0]) func = *argv; else { --argc; func = *++argv; }
      flag[FUNCTION] = true; break;
    default: sprintf (buf, usage, prog); cout<<buf; exit(EXIT_FAILURE); break;
    }

  for (i = 0; i < FLAG_MAX; i++) sum += (flag[i]) ? 1 : 0;
  if (!sum) flag[ENERGY] = true;

  if   (!session)  message (prog, "no session file", ERROR);
  if   (argc != 1) message (prog, "no field file",   ERROR);
  else             dump = *argv;
}


static bool getDump (Domain*   D   ,
		     ifstream& dump)
// ---------------------------------------------------------------------------
// Read next set of field dumps from file.
// ---------------------------------------------------------------------------
{
  dump >> *D;
  return dump.good ();
}


static void putDump  (Domain*            D       ,
		      vector<AuxField*>& addField,
		      int_t              nOut    ,
		      ostream&           strm    )
// ---------------------------------------------------------------------------
// This is a version of the normal Domain dump that adds extra AuxFields.
// ---------------------------------------------------------------------------
{
  const char *hdr_fmt[] = { 
    "%-25s "    "Session\n",
    "%-25s "    "Created\n",
    "%-25s "    "Nr, Ns, Nz, Elements\n",
    "%-25d "    "Step\n",
    "%-25.6g "  "Time\n",
    "%-25.6g "  "Time step\n",
    "%-25.6g "  "Kinvis\n",
    "%-25.6g "  "Beta\n",
    "%-25s "    "Fields written\n",
    "%-25s "    "Format\n"
  };

  int_t  i, nComponent;
  char   routine[] = "putDump";
  char   s1[StrMax], s2[StrMax];
  time_t tp (::time (0));

  if (D -> nField() == 1)	// -- Scalar original field.
    nComponent = 1;
  else				// -- Original field was vector.
    nComponent = (D -> nField() == 3) ? 2 : 3;

  sprintf (s1, hdr_fmt[0], D -> name);
  strm << s1;

  strftime (s2, 25, "%a %b %d %H:%M:%S %Y", localtime (&tp));
  sprintf  (s1, hdr_fmt[1], s2);
  strm << s1;

  D -> u[0] -> describe (s2);
  sprintf (s1, hdr_fmt[2], s2);
  strm << s1;

  sprintf (s1, hdr_fmt[3], D -> step);
  strm << s1;

  sprintf (s1, hdr_fmt[4], D -> time);
  strm << s1;

  sprintf (s1, hdr_fmt[5], Femlib::value ("D_T"));
  strm << s1;

  sprintf (s1, hdr_fmt[6], Femlib::value ("KINVIS"));
  strm << s1;

  sprintf (s1, hdr_fmt[7], Femlib::value ("BETA"));
  strm << s1;

  for (i = 0; i <= nComponent; i++) s2[i] = D -> u[i] -> name();
  for (i = 0; i <  nOut; i++)
    s2[nComponent + i + 1] = addField[i] -> name();
  s2[nComponent + nOut + 1] = '\0';

  sprintf (s1, hdr_fmt[8], s2);
  strm << s1;

  sprintf (s2, "binary ");
  Veclib::describeFormat (s2 + strlen (s2));
  sprintf (s1, hdr_fmt[9], s2);
  strm << s1;

  for (i = 0; i <= nComponent; i++) strm << *D -> u[i];
  for (i = 0; i < nOut; i++) strm << *addField[i];

  if (!strm) message (routine, "failed writing field file", ERROR);
  strm << flush;
}
