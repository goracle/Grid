    /*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid 

    Source file: ./lib/algorithms/iterative/SchurRedBlack.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
    *************************************************************************************/
    /*  END LEGAL */
#ifndef GRID_SCHUR_RED_BLACK_H
#define GRID_SCHUR_RED_BLACK_H

  /*
   * Red black Schur decomposition
   *
   *  M = (Mee Meo) =  (1             0 )   (Mee   0               )  (1 Mee^{-1} Meo)
   *      (Moe Moo)    (Moe Mee^-1    1 )   (0   Moo-Moe Mee^-1 Meo)  (0   1         )
   *                =         L                     D                     U
   *
   * L^-1 = (1              0 )
   *        (-MoeMee^{-1}   1 )   
   * L^{dag} = ( 1       Mee^{-dag} Moe^{dag} )
   *           ( 0       1                    )
   * L^{-d}  = ( 1      -Mee^{-dag} Moe^{dag} )
   *           ( 0       1                    )
   *
   * U^-1 = (1   -Mee^{-1} Meo)
   *        (0    1           )
   * U^{dag} = ( 1                 0)
   *           (Meo^dag Mee^{-dag} 1)
   * U^{-dag} = (  1                 0)
   *            (-Meo^dag Mee^{-dag} 1)
   ***********************
   *     M psi = eta
   ***********************
   *Odd
   * i)                 D_oo psi_o =  L^{-1}  eta_o
   *                        eta_o' = (D_oo)^dag (eta_o - Moe Mee^{-1} eta_e)
   *
   * Wilson:
   *      (D_oo)^{\dag} D_oo psi_o = (D_oo)^dag L^{-1}  eta_o
   * Stag:
   *      D_oo psi_o = L^{-1}  eta =    (eta_o - Moe Mee^{-1} eta_e)
   *
   * L^-1 eta_o= (1              0 ) (e
   *             (-MoeMee^{-1}   1 )   
   *
   *Even
   * ii)  Mee psi_e + Meo psi_o = src_e
   *
   *   => sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
   *
   * 
   * TODO: Other options:
   * 
   * a) change checkerboards for Schur e<->o
   *
   * Left precon by Moo^-1
   * b) Doo^{dag} M_oo^-dag Moo^-1 Doo psi_0 =  (D_oo)^dag M_oo^-dag Moo^-1 L^{-1}  eta_o
   *                              eta_o'     = (D_oo)^dag  M_oo^-dag Moo^-1 (eta_o - Moe Mee^{-1} eta_e)
   *
   * Right precon by Moo^-1
   * c) M_oo^-dag Doo^{dag} Doo Moo^-1 phi_0 = M_oo^-dag (D_oo)^dag L^{-1}  eta_o
   *                              eta_o'     = M_oo^-dag (D_oo)^dag (eta_o - Moe Mee^{-1} eta_e)
   *                              psi_o = M_oo^-1 phi_o
   * TODO: Deflation 
   */
namespace Grid {

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Now make the norm reflect extra factor of Mee
  template<class Field> class SchurRedBlackStaggeredSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackStaggeredSolve(OperatorFunction<Field> &HermitianRBSolver)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
    };

    template<class Matrix>
      void operator() (Matrix & _Matrix,const Field &in, Field &out){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurStaggeredOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);
      
      std::cout << GridLogMessage << " SchurRedBlackStaggeredSolve " <<std::endl;
      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);

      std::cout << GridLogMessage << " SchurRedBlackStaggeredSolve checkerboards picked" <<std::endl;
    
      /////////////////////////////////////////////////////
      // src_o = (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      //src_o = tmp;     assert(src_o.checkerboard ==Odd);
      _Matrix.Mooee(tmp,src_o); // Extra factor of "m" in source from dumb choice of matrix norm.

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver calling the Mpc solver" <<std::endl;
      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver called  the Mpc solver" <<std::endl;

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver reconstructed other CB" <<std::endl;
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );
      std::cout<<GridLogMessage << "SchurRedBlackStaggeredSolver inserted solution" <<std::endl;

      // Verify the unprec residual
      _Matrix.M(out,resid); 
      resid = resid-in;
      RealD ns = norm2(in);
      RealD nr = norm2(resid);
      std::cout<<GridLogMessage << "SchurRedBlackStaggered solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
    }     
  };
  template<class Field> using SchurRedBlackStagSolve = SchurRedBlackStaggeredSolve<Field>;

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagMooeeSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagMooeeSolve(OperatorFunction<Field> &HermitianRBSolver,int cb=0)  :  _HermitianRBSolver(HermitianRBSolver) 
  { 
    CBfactorise=cb;
  };
    template<class Matrix>
      void operator() (Matrix & _Matrix,const Field &in, Field &out){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagMooeeOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      _Matrix.M(out,resid); 
      resid = resid-in;
      RealD ns = norm2(in);
      RealD nr = norm2(resid);

      std::cout<<GridLogMessage << "SchurRedBlackDiagMooee solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
    }     
  };


  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoSolve {
  private:
    OperatorFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoSolve(OperatorFunction<Field> &HermitianRBSolver)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
    };

    template<class Matrix>
      void operator() (Matrix & _Matrix,const Field &in, Field &out){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
//      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      _Matrix.MooeeInv(tmp,sol_o);        assert(  sol_o.checkerboard   ==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      _Matrix.M(out,resid); 
      resid = resid-in;
      RealD ns = norm2(in);
      RealD nr = norm2(resid);

      std::cout<<GridLogMessage << "SchurRedBlackDiagTwo solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
    }     
  };
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class Field> class SchurRedBlackDiagTwoMixed {
  private:
    LinearFunction<Field> & _HermitianRBSolver;
    int CBfactorise;
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
  SchurRedBlackDiagTwoMixed(LinearFunction<Field> &HermitianRBSolver)  :
     _HermitianRBSolver(HermitianRBSolver) 
    { 
      CBfactorise=0;
    };

    template<class Matrix>
      void operator() (Matrix & _Matrix,const Field &in, Field &out){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      Field src_e(grid);
      Field src_o(grid);
      Field sol_e(grid);
      Field sol_o(grid);
      Field   tmp(grid);
      Field  Mtmp(grid);
      Field resid(fgrid);

      pickCheckerboard(Even,src_e,in);
      pickCheckerboard(Odd ,src_o,in);
      pickCheckerboard(Even,sol_e,out);
      pickCheckerboard(Odd ,sol_o,out);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e,tmp);     assert(  tmp.checkerboard ==Even);
      _Matrix.Meooe   (tmp,Mtmp);      assert( Mtmp.checkerboard ==Odd);     
      tmp=src_o-Mtmp;                  assert(  tmp.checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp,src_o);     assert(src_o.checkerboard ==Odd);       

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
//      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
//      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      _HermitianRBSolver(src_o,tmp);  assert(tmp.checkerboard==Odd);
      _Matrix.MooeeInv(tmp,sol_o);        assert(  sol_o.checkerboard   ==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp);        assert(  tmp.checkerboard   ==Even);
      src_e = src_e-tmp;               assert(  src_e.checkerboard ==Even);
      _Matrix.MooeeInv(src_e,sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out,sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out,sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      _Matrix.M(out,resid); 
      resid = resid-in;
      RealD ns = norm2(in);
      RealD nr = norm2(resid);

      std::cout<<GridLogMessage << "SchurRedBlackDiagTwo solver true unprec resid "<< std::sqrt(nr/ns) <<" nr "<< nr <<" ns "<<ns << std::endl;
    }     
  };

  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  // Take a matrix and form a Red Black solver calling a Herm solver
  // Use of RB info prevents making SchurRedBlackSolve conform to standard interface
  ///////////////////////////////////////////////////////////////////////////////////////////////////////
  template<class LinopPolicyD, class LinopPolicyF,typename std::enable_if<getPrecision<typename LinopPolicyD::FermionField>::value == 2 && getPrecision<typename LinopPolicyF::FermionField>::value == 1 , int>::type = 0> class SplitConjugateGradientReliableUpdate;
}
namespace cps{
  template<class Field, typename GridPolicies> class deflateGuessRead;
}
namespace Grid{
  template<class Field, class LinopPolicyD, class LinopPolicyF, typename GridPolicies> class SchurRedBlackDiagTwoSplit {
  private:

    //SplitConjugateGradientReliableUpdate
    Grid::SplitConjugateGradientReliableUpdate<LinopPolicyD, LinopPolicyF> & _HermitianRBSolver;
    int CBfactorise;
    cps::deflateGuessRead<Field, GridPolicies> & _guesser; //assumes single prec. evecs
  public:

    /////////////////////////////////////////////////////
    // Wrap the usual normal equations Schur trick
    /////////////////////////////////////////////////////
    SchurRedBlackDiagTwoSplit(SplitConjugateGradientReliableUpdate<LinopPolicyD, LinopPolicyF> &HermitianRBSolver, cps::deflateGuessRead<Field, GridPolicies> &guesser):
      _HermitianRBSolver(HermitianRBSolver),
      _guesser(guesser)
    { 
      CBfactorise=0;
    };


    //split cg version
    template<class Matrix>
    void operator() (Matrix & _Matrix,const std::vector<Field> &in, std::vector<Field> &out, bool do_defl=false, bool defl_sub=false){

      // FIXME CGdiagonalMee not implemented virtual function
      // FIXME use CBfactorise to control schur decomp
      GridBase *grid = _Matrix.RedBlackGrid();
      GridBase *fgrid= _Matrix.Grid();
      const size_t nsolve = in.size();

      SchurDiagTwoOperator<Matrix,Field> _HermOpEO(_Matrix);
 
      std::vector<Field> src_e(nsolve, Field(grid));
      std::vector<Field> src_o(nsolve, Field(grid));
      std::vector<Field> defl(defl_sub? nsolve : 0, Field(grid));
      Field sol_e(grid);
      Field sol_o(grid);
      std::vector<Field>   tmp(nsolve, Field(grid));
      Field  Mtmp(grid);
      Field resid(fgrid);
      Field resid_rb(grid);

      for(size_t s=0; s<nsolve; s++){
	pickCheckerboard(Even,src_e[s],in[s]);
	pickCheckerboard(Odd ,src_o[s],in[s]);
	pickCheckerboard(Even,sol_e,out[s]);
	pickCheckerboard(Odd ,sol_o,out[s]);
    
      /////////////////////////////////////////////////////
      // src_o = Mdag * (source_o - Moe MeeInv source_e)
      /////////////////////////////////////////////////////
      _Matrix.MooeeInv(src_e[s],tmp[s]);     assert(  tmp[s].checkerboard ==Even);
      _Matrix.Meooe   (tmp[s],Mtmp);      assert( Mtmp.checkerboard ==Odd);
      tmp[s]=src_o[s]-Mtmp;                  assert(  tmp[s].checkerboard ==Odd);     

      // get the right MpcDag
      _HermOpEO.MpcDag(tmp[s],src_o[s]);     assert(src_o[s].checkerboard ==Odd);
      }
      //////////////////////////////////////////////////////////////
      // get the subtraction terms for A2A propagator
      //////////////////////////////////////////////////////////////
      if(do_defl && !defl_sub){
	_guesser(src_o, tmp);
      }else if(do_defl && defl_sub){
	_guesser(src_o, tmp, defl);
      }else if(!do_defl && defl_sub){
	_guesser(src_o, defl);
      }

      //////////////////////////////////////////////////////////////
      // Call the red-black solver
      //////////////////////////////////////////////////////////////
      std::cout<<GridLogMessage << "SchurRedBlack solver calling the MpcDagMp solver" <<std::endl;
//      _HermitianRBSolver(_HermOpEO,src_o,sol_o);  assert(sol_o.checkerboard==Odd);
//      _HermitianRBSolver(_HermOpEO,src_o,tmp);  assert(tmp.checkerboard==Odd);
      _HermitianRBSolver(src_o,tmp);
      for(size_t s=0; s<nsolve; s++){
	assert(tmp[s].checkerboard==Odd);
	//check residual
	{
	  _HermOpEO.HermOp(tmp[s], resid_rb);
	  resid_rb = resid_rb-src_o[s];
	  RealD ns = norm2(src_o[s]);
	  RealD nr = norm2(resid_rb);
	  std::cout<<GridLogMessage << "SchurRedBlackDiagTwoSplit solver true post-solve resid [" << s << "] = " << std::sqrt(nr/ns) <<" nr = "<< nr <<" ns = "<<ns << std::endl;
	}

	//Pull low-mode part out of solution
	if(defl_sub){
	  assert(defl[s].checkerboard==Odd);
	  std::cout<<GridLogMessage << "SchurRedBlackDiagTwoSplit, subtracting low mode [" << s << "]." <<std::endl;
	  axpy(tmp[s], -1.0, defl[s], tmp[s]);
	}
	_Matrix.MooeeInv(tmp[s],sol_o);        assert(  sol_o.checkerboard   ==Odd);

      ///////////////////////////////////////////////////
      // sol_e = M_ee^-1 * ( src_e - Meo sol_o )...
      ///////////////////////////////////////////////////
      _Matrix.Meooe(sol_o,tmp[s]);        assert(  tmp[s].checkerboard   ==Even);
      src_e[s] = src_e[s]-tmp[s];               assert(  src_e[s].checkerboard ==Even);
      _Matrix.MooeeInv(src_e[s],sol_e);   assert(  sol_e.checkerboard ==Even);
     
      setCheckerboard(out[s],sol_e); assert(  sol_e.checkerboard ==Even);
      setCheckerboard(out[s],sol_o); assert(  sol_o.checkerboard ==Odd );

      // Verify the unprec residual
      _Matrix.M(out[s],resid); 
      resid = resid-in[s];
      RealD ns = norm2(in[s]);
      RealD nr = norm2(resid);

      std::cout<<GridLogMessage << "SchurRedBlackDiagTwoSplit solver true unprec resid [" << s << "] = " << std::sqrt(nr/ns) <<" nr = "<< nr <<" ns = "<<ns << std::endl;
      }
    }     


  };

}
#endif
