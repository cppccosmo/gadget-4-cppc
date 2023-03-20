////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//    Copyright 2022 Amol Upadhye                                             //
//                                                                            //
//    This file is part of MuFLR.                                             //
//                                                                            //
//    MuFLR is free software: you can redistribute it and/or modify           //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation, either version 3 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    MuFLR is distributed in the hope that it will be useful,                //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with MuFLR.  If not, see <http://www.gnu.org/licenses/>.          //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//USER-DEFINED GRID PARAMETERS
#define NK (256)
#define KMIN (1e-4)
#define KMAX (20)

////////////////////////////////////////////////////////////////////////////////
//k grid and windowing: NK is number of points in "real" k grid,                 
//I extend this by a factor of four for fast-pt.  The extended grid: 
//   -3NK/2  <= i < -17NK/16 : P=0 
//  -17NK/16 <= i < -NK      : P tapered smoothly from 0 to P[0] 
//     -NK   <= i < 0        : Power spectrum extrapolated to left 
//       0   <= i < NK       : Power spectrum region of interest 
//      NK   <= i < 2NK      : Power spectrum extrapolated to right 
//     2NK   <= i < 33NK/16  : P tapered smoothly from P[NK-1] to 0 
//  33NK/16  <= i < 5NK/2    : P=0

#define NKP (4*NK)
#define NSHIFT ((NKP-NK)/2)
#define DLNK (log(KMAX/KMIN)/(NK-1.0+1e-300))
#define LNKMIN (log(KMIN))

//split up interval between zero-padding and tapering:
//these are measured in units of NK/16, and must add to (NKP/NK-1)*16
#define S_PADL (7)
#define S_TAPL (1)
#define S_EXTL (16)
#define S_EXTR (16)
#define S_TAPR (1)
#define S_PATR (7)

#define LNK_PAD_MIN   (LNKMIN        - DLNK*NSHIFT)
#define LNK_PAD_WINLO (LNK_PAD_MIN   + DLNK*NK*S_PADL/16)
#define LNK_PAD_WINLI (LNK_PAD_WINLO + DLNK*NK*S_TAPL/16)
#define LNK_PAD_WINRI (LNK_PAD_WINLI + DLNK*(NK*(16+S_EXTL+S_EXTR)/16 -1))
#define LNK_PAD_WINRO (LNK_PAD_WINRI + DLNK*NK*S_TAPR/16)
#define LNK_PAD_MAX   (LNK_PAD_WINRO + DLNK*NK*S_PADR/16)

