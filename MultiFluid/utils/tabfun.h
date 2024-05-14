//    Copyright 2023 Amol Upadhye
//
//    This file is part of FlowsForTheMasses.
//
//    FlowsForTheMasses is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    FlowsForTheMasses is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with FlowsForTheMasses.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <fstream>

const int TF_DEBUG_ALLOC = 0;

///////////////////////////////////////////////////////////////////////////////

class tabulated_function{
 public:
  tabulated_function(const char *filename, int nCol, 
		     const int nvars, const int nrx, const int nrf){

    //initialize x and f tables
    tableSize = nCol;
    tableSizeX = 0;
    tableSizeY = 0;
    xTable = new double[tableSize];
    fTable = new double[tableSize];

    if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated_function: constructor0: allocated x and f "
                << "arrays of size " << tableSize << std::endl;

    //read from file
    int n=0;
    double temp_in[nvars];

    std::ifstream tableFile(filename, std::ios::in);
    int status = 1;
    do{
      tabulated_function::discard_comments(&tableFile);
      for(int i=0; i<nvars; i++) status = status && (tableFile >> temp_in[i]); 
      xTable[n] = temp_in[nrx];
      fTable[n] = temp_in[nrf];
      n++;
    } while(status && n<nCol);

    tableSize = n;

    tableFile.close();
  }

  tabulated_function(int nCol, const double *xArray, const double *fArray){

    //initialize x and f tables
    tableSize = nCol;
    tableSizeX = 0;
    tableSizeY = 0;
    xTable = new double[tableSize];
    fTable = new double[tableSize];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: constructor1: allocated x and f "
	      << "arrays of size " << tableSize << std::endl;

    //input from arrays
    for(int i=0; i<nCol; i++){
      xTable[i] = xArray[i];
      fTable[i] = fArray[i];
    }
  }

  tabulated_function(int nCol){
    
    //initialize x and f tables
    tableSize = nCol;
    tableSizeX = 0;
    tableSizeY = 0;
    xTable = new double[tableSize];
    fTable = new double[tableSize];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: constructor2: allocated x and f "
	      << "arrays of size " << tableSize << std::endl;

    //fill with zeros for now
    for(int i=0; i<nCol; i++){
      xTable[i] = 0;
      fTable[i] = 0;
    }
  }

  tabulated_function(){
    //do nothing yet; initialize later
    tableSize=0;
    tableSizeX=0;
    tableSizeY=0;
  }

  tabulated_function(const char *filename, int nCx, int nCy, const int nvars, 
                     const int nrx, const int nry, const int nrf){

    //initialize x and f tables
    tableSize = 0;
    tableSizeX = nCx;
    tableSizeY = nCy;
    xTable = new double[tableSizeX];
    yTable = new double[tableSizeY];
    fTable = new double[tableSizeX * tableSizeY];

    if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated_function: constructor4: allocated x, y, and f "
                << "arrays of size " << tableSizeX << " " << tableSizeY 
                << std::endl;

    //read from file
    int n=0;
    double temp_in[nvars];

    std::ifstream tableFile(filename, std::ios::in);
    int status = 1;
    do{
      tabulated_function::discard_comments(&tableFile);
      for(int i=0; i<nvars; i++) status = status && (tableFile >> temp_in[i]);
      xTable[n/nCy] = temp_in[nrx];
      yTable[n%nCy] = temp_in[nry];
      fTable[n]     = temp_in[nrf];
      n++;
    } while(status && n<nCx*nCy);

    tableFile.close();
  }

  tabulated_function(int nColx, int nColy, const double *xArray, 
		     const double *yArray, const double *fArray){

    //initialize x, y, and f tables
    tableSizeX = nColx;
    tableSizeY = nColy;
    tableSize = 0;
    xTable = new double[tableSizeX];
    yTable = new double[tableSizeY];
    fTable = new double[tableSizeX * tableSizeY];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: allocated x, y, and f arrays of size "
              << tableSizeX << " " << tableSizeY << std::endl;

    //fill tables from inputs
    for(int i=0; i<tableSizeX; i++) xTable[i] = xArray[i];
    for(int i=0; i<tableSizeY; i++) yTable[i] = yArray[i];
    for(int i=0; i<tableSizeX*tableSizeY; i++) fTable[i] = fArray[i];
  }

  tabulated_function(int nColx, int nColy){

    //initialize x, y, and f tables
    tableSizeX = nColx;
    tableSizeY = nColy;
    tableSize = 0;
    xTable = new double[tableSizeX];
    yTable = new double[tableSizeY];
    fTable = new double[tableSizeX * tableSizeY];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: allocated x, y, and f arrays of size "
              << tableSizeX << " " << tableSizeY << std::endl;
  }

  int initialize(int nCol){ //just allocate memory

    //initialize x and f tables
    tableSize = nCol;
    tableSizeX = 0;
    tableSizeY = 0;
    xTable = new double[tableSize];
    fTable = new double[tableSize];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: constructor2: allocated x and f "
              << "arrays of size " << tableSize << std::endl;

    //fill with zeros for now
    for(int i=0; i<nCol; i++){
      xTable[i] = 0;
      fTable[i] = 0;
    }

    return 0;
  }

  int initialize(int nCol, const double *xArray, const double *fArray){

    //initialize x and f tables  
    tableSize = nCol;
    tableSizeX = 0;
    tableSizeY = 0;
    xTable = new double[tableSize];
    fTable = new double[tableSize];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: initialize: allocated x and f "
	      << "arrays of size " << tableSize << std::endl;

    //fill tables from inputs  
    for(int i=0; i<tableSize; i++) xTable[i] = xArray[i];
    for(int i=0; i<tableSize; i++) fTable[i] = fArray[i];

    return 0;
  }

  int initialize(const char *filename, int nCx, int nCy, const int nvars,
                 const int nrx, const int nry, const int nrf){

    //initialize x and f tables
    tableSize = 0;
    tableSizeX = nCx;
    tableSizeY = nCy;
    xTable = new double[tableSizeX];
    yTable = new double[tableSizeY];
    fTable = new double[tableSizeX * tableSizeY];

    if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated_function: constructor4: allocated x, y, and f "
                << "arrays of size " << tableSizeX << " " << tableSizeY
                << std::endl;

    //read from file
    int n=0;
    double temp_in[nvars];

    std::ifstream tableFile(filename, std::ios::in);
    int status = 1;
    do{
      tabulated_function::discard_comments(&tableFile);
      for(int i=0; i<nvars; i++) status = status && (tableFile >> temp_in[i]);
      xTable[n/nCy] = temp_in[nrx];
      yTable[n%nCy] = temp_in[nry];
      fTable[n]     = temp_in[nrf];
      n++;
    } while(status && n<nCx*nCy);

    tableFile.close();
    return 0;
  }

  int initialize(int nColx, int nColy, const double *xArray, 
		 const double *yArray, const double *fArray){

    //initialize x, y, and f tables
    tableSizeX = nColx;
    tableSizeY = nColy;
    tableSize = 0;
    xTable = new double[tableSizeX];
    yTable = new double[tableSizeY];
    fTable = new double[tableSizeX * tableSizeY];

    if(TF_DEBUG_ALLOC)
    std::cout << "#tabulated_function: allocated x, y, and f arrays of size "
              << tableSizeX << " " << tableSizeY << std::endl;

    //fill tables from inputs
    for(int i=0; i<tableSizeX; i++) xTable[i] = xArray[i];
    for(int i=0; i<tableSizeY; i++) yTable[i] = yArray[i];
    for(int i=0; i<tableSizeX*tableSizeY; i++) fTable[i] = fArray[i];

    return 0;
  }

  int input_arrays(const double *xArray, const double *fArray){
    //input from arrays
    for(int i=0; i<tableSize; i++){
      xTable[i] = xArray[i];
      fTable[i] = fArray[i];
    }
    return 0;
  }

  int input_arrays(const double *xArray, const double *yArray, 
		   const double *fArray){
    //input from arrays
    for(int i=0; i<tableSizeX; i++) xTable[i] = xArray[i];
    for(int j=0; j<tableSizeY; j++) yTable[j] = yArray[j];
    for(int i=0; i<tableSizeX; i++){
      for(int j=0; j<tableSizeY; j++){ 
	int k = nf(i,j);
	fTable[k] = fArray[k]; 
      }
    }
    return 0;
  }

  double f(double x) const {
    int n1 = tabulated_function::findN(x);
    if(n1<=0) //allows for extrapolation to the left
      return tabulated_function::linInterp(&xTable[0],&fTable[0],x);
    else if(n1>=tableSize-2) //allows for extrapolation to the right
      return tabulated_function::linInterp(&xTable[tableSize-2],
					   &fTable[tableSize-2],x);
    else 
      return tabulated_function::cubicInterp(&xTable[n1-1],&fTable[n1-1],x);
  }

  double f(double x, double y) const {
    int nx = tabulated_function::findNx_2d(x);
    int ny = tabulated_function::findNy_2d(y);
    if(nx<0 || nx>=tableSizeX){ 
      std::cout << "ERROR: x=" << x << " out of bounds; nx=" << nx 
		<< ", tableSizes=" << tableSizeX << ' ' << tableSizeY
		<<std::endl;
      std::abort(); }
    if(ny<0 || ny>=tableSizeY){
      std::cout << "ERROR: ny=" << ny << " out of bounds\n"; 
      std::cout << "#      tableSizeX=" << tableSizeX << ", tableSizeY=" 
		<< tableSizeY << std::endl;
      std::abort(); 
    }

    double fTabY[] = {0,0,0,0}; //f[] vs. y[] at interpolated x

    if(nx>0 && nx<tableSizeX-2){ //interior: use cubic interpolation
      double fTabA[] = { fTable[nf(nx-1,ny-1)], fTable[nf(nx,ny-1)],
			 fTable[nf(nx+1,ny-1)], fTable[nf(nx+2,ny-1)] };
      double fTabB[] = { fTable[nf(nx-1,ny)],   fTable[nf(nx,ny)],
			 fTable[nf(nx+1,ny)],   fTable[nf(nx+2,ny)]   };
      double fTabC[] = { fTable[nf(nx-1,ny+1)], fTable[nf(nx,ny+1)],
			 fTable[nf(nx+1,ny+1)], fTable[nf(nx+2,ny+1)] };
      double fTabD[] = { fTable[nf(nx-1,ny+2)], fTable[nf(nx,ny+2)],
			 fTable[nf(nx+1,ny+2)], fTable[nf(nx+2,ny+2)] };
      fTabY[0] = tabulated_function::cubicInterp(&xTable[nx-1],fTabA,x);
      fTabY[1] = tabulated_function::cubicInterp(&xTable[nx-1],fTabB,x);
      fTabY[2] = tabulated_function::cubicInterp(&xTable[nx-1],fTabC,x);
      fTabY[3] = tabulated_function::cubicInterp(&xTable[nx-1],fTabD,x);
      
      if(ny>0 && ny<tableSizeY-2)
	return tabulated_function::cubicInterp(&yTable[ny-1],fTabY,y);
      else if(ny==0)
        return tabulated_function::linInterp(&yTable[ny],&fTabY[1],y);
      else if(ny==tableSizeY-2)
	return tabulated_function::linInterp(&yTable[ny],&fTabY[1],y);
      else{
	std::cout << "#tabulated_function::f: ERROR: extrapolation in y." 
		  << std::endl;
	abort();
      }
    }
    else{ //boundary: use linear interpolation
      double fTabA[] = { fTable[nf(nx,ny-1)], fTable[nf(nx+1,ny-1)] } ;
      double fTabB[] = { fTable[nf(nx,ny)],   fTable[nf(nx+1,ny)]   } ;
      double fTabC[] = { fTable[nf(nx,ny+1)], fTable[nf(nx+1,ny+1)] } ;
      double fTabD[] = { fTable[nf(nx,ny+2)], fTable[nf(nx+1,ny+2)] } ;

      fTabY[0] = tabulated_function::linInterp(&xTable[nx],fTabA,x);
      fTabY[1] = tabulated_function::linInterp(&xTable[nx],fTabB,x);
      fTabY[2] = tabulated_function::linInterp(&xTable[nx],fTabC,x);
      fTabY[3] = tabulated_function::linInterp(&xTable[nx],fTabD,x);
      
      if(ny>0 && ny<tableSizeY-2)
	return tabulated_function::cubicInterp(&yTable[ny-1],fTabY,y);
      else if(ny==0)
	return tabulated_function::linInterp(&yTable[ny],&fTabY[1],y);
      else if(ny==tableSizeY-2)
	return tabulated_function::linInterp(&yTable[ny],&fTabY[1],y);
      else{
	std::cout << "#tabulated_function::f: ERROR: extrapolation in y."
                  << std::endl;
        abort();
      }
    }
  }

  double operator()(double x) const { return f(x); }
  double operator()(double x, double y) const { return f(x,y); }

  double df(double x) const {
    int n1 = tabulated_function::findN(x);
    return tabulated_function::df(n1);
  }
 
  double df(int n1) const {
    if(n1 == tableSize-1) 
      return (fTable[n1]-fTable[n1-1])/(xTable[n1]-xTable[n1-1]);
    else
      return (fTable[n1+1]-fTable[n1])/(xTable[n1+1]-xTable[n1]);
  }

  double df(double x, double y, double theta_dir) const {

    //find nx and ny
    int nx = tabulated_function::findNx_2d(x);
    int ny = tabulated_function::findNy_2d(y);
    if(nx<0 || nx>=tableSizeX){
      std::cout << "ERROR: x out of bounds\n"; std::abort(); }
    if(ny<0 || ny>=tableSizeY){
      std::cout << "ERROR: ny=" << ny << " out of bounds\n"; std::abort(); }

    //linear interpolation to find gradient
    int ixmax = (nx==tableSizeX-1), iymax = (ny==tableSizeY-1);
    double f00 = fTable[(nx-ixmax)*tableSizeY + ny-iymax];
    double f01 = fTable[(nx-ixmax)*tableSizeY + ny+1-iymax];
    double f10 = fTable[(nx+1-ixmax)*tableSizeY + ny-iymax];
    double f11 = fTable[(nx+1-ixmax)*tableSizeY + ny+1-iymax];
    
    double Dx0 = xTable[nx+1-ixmax] - xTable[nx-ixmax];
    double Dy0 = yTable[ny+1-iymax] - yTable[ny-iymax];

    double dfdxTab[2] = {(f10-f00) / Dx0, (f11-f01) / Dx0};
    double dfdyTab[2] = {(f01-f00) / Dy0, (f11-f10) / Dy0};
    double dfdx = linInterp(&yTable[ny-iymax],dfdxTab,y);
    double dfdy = linInterp(&xTable[nx-ixmax],dfdyTab,x);

    //output derivative in theta_dir direction
    return dfdx*cos(theta_dir) + dfdy*sin(theta_dir);
  }

  double d2f(double x) const {
    int n1 = tabulated_function::findN(x);
    return tabulated_function::d2f(n1);
  }

  double d2f(int n1) const {
    return 2.0/(xTable[n1+1]-xTable[n1-1]) 
           * ( tabulated_function::df(n1) - tabulated_function::df(n1-1) );
  }

  int clear(){
    //clean up before re-initializing arrays
    if(tableSize+tableSizeX > 0){
      if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated function: deleting x Table of dim " 
		<< tableSize+tableSizeX << std::endl;
      delete [] xTable;
    }
    if(tableSizeY > 0){
      if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated function: deleting y Table of dim " 
		<< tableSizeY << std::endl;
      delete [] yTable;
    }
    if(tableSize + tableSizeX + tableSizeY > 0){
      if(TF_DEBUG_ALLOC)
      std::cout << "#tabulated function: deleting f Table of dim " 
		<< tableSize + tableSizeX*tableSizeY << std::endl;
      delete [] fTable;
      tableSize = tableSizeX = tableSizeY = 0;
    }
    return 0;
  }

  ~tabulated_function(){
    clear();
  }


 private:
  int tableSize, tableSizeX, tableSizeY;
  double *xTable, *yTable;
  double *fTable;

  inline double fmax(double x, double y) const { return( (x>y)?(x):(y) ); }
  inline double fmin(double x, double y) const { return( (x<y)?(x):(y) ); }
  inline double fmax(double *x, int len) const {
    return( (len==2)?(fmax(x[0],x[1])):(fmax(fmax(x,len-1),x[len-1])) ); }
  inline double fmin(double *x, int len) const {
    return( (len==2)?(fmin(x[0],x[1])):(fmin(fmin(x,len-1),x[len-1])) ); }
  inline void discard_comments(std::ifstream *file) const {
    while(file->peek()=='#' || file->peek()=='\n'){file->ignore(10000,'\n');} }
  
  inline int nf(int nx, int ny) const { return ny + tableSizeY*nx; }

  double linInterp(double *xTab, double *fTab, double xEval) const {
    double F = fTab[0] + (fTab[1]-fTab[0])/(xTab[1]-xTab[0])*(xEval-xTab[0]);
    return F;
  }

  //cubic polynomial interpolation
  double cubicInterp(double *xTab, double *fTab, double xEval) const {
    //Input the function f(x), in the form of a table of 4 values fTab[]
    //given at the four points xTab[].  Output the interpolated value
    //of the function at xEval, using a cubic polynomial interpolation.
    
    //make sure that we're interpolating, not extrapolating!
    if(xEval<fmin(xTab,4) || xEval>fmax(xTab,4)){
      std::cout << "cubicInterp WARNING: xEval=" << xEval
		<< " out of bounds.  You are" 
		<< std::endl
		<< "                     extrapolating, not interpolating!" 
		<< std::endl;
      abort();
    }
    
    double F = (xEval-xTab[1])*(xEval-xTab[2])*(xEval-xTab[3])
      /(xTab[0]-xTab[1])/(xTab[0]-xTab[2])/(xTab[0]-xTab[3])*fTab[0] 
      + (xEval-xTab[0])*(xEval-xTab[2])*(xEval-xTab[3])
      /(xTab[1]-xTab[0])/(xTab[1]-xTab[2])/(xTab[1]-xTab[3])*fTab[1] 
      + (xEval-xTab[0])*(xEval-xTab[1])*(xEval-xTab[3])
      /(xTab[2]-xTab[0])/(xTab[2]-xTab[1])/(xTab[2]-xTab[3])*fTab[2] 
      + (xEval-xTab[0])*(xEval-xTab[1])*(xEval-xTab[2])
      /(xTab[3]-xTab[0])/(xTab[3]-xTab[1])/(xTab[3]-xTab[2])*fTab[3];
    return(F);
  }

  int findN(double x) const { 
    //do something simple for now; can use a spiffier algorithm if 
    //large tables are necessary
    int n = 0;
    while(xTable[n+1]<x && n<tableSize-2){n++;}
    return n;
  }

  int findNx_2d(double x) const { 
    //do something simple for now; can use a spiffier algorithm if 
    //large tables are necessary
    int n = 0;
    while(xTable[n+1]<x && n<tableSizeX-2){n++;}
    return n;
  }

  int findNy_2d(double y) const { 
    //do something simple for now; can use a spiffier algorithm if 
    //large tables are necessary
    int n = 0;
    while(yTable[n+1]<y && n<tableSizeY-2){n++;}
    return n;
  }

};

