#!/usr/bin/env python

#---------------------------------------------------------------------------
## pythonFlu - Python wrapping for OpenFOAM C++ API
## Copyright (C) 2010- Alexey Petrov
## Copyright (C) 2009-2010 Pebble Bed Modular Reactor (Pty) Limited (PBMR)
## 
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
## 
## See http://sourceforge.net/projects/pythonflu
##
## Author : Alexey PETROV
##


#---------------------------------------------------------------------------
from Foam import man, ref


#---------------------------------------------------------------------------
def createFields( runTime, mesh):
    
    ref.ext_Info() << "Reading thermophysical properties\n" << ref.nl
  
    pThermo = man.basicPsiThermo.New( mesh )
  
    p = man.volScalarField( pThermo.p(), man.Deps( pThermo ) )
    e = man.volScalarField( pThermo.e(), man.Deps( pThermo ) )
    psi = man.volScalarField( pThermo.psi(), man.Deps( pThermo ) )

    rho = man.volScalarField( man.IOobject( ref.word( "rho" ),
                                            ref.fileName( runTime.timeName() ),
                                            mesh ),
                              man.volScalarField( pThermo.rho(), man.Deps( pThermo ) ) )

    ref.ext_Info() << "Reading field U\n" << ref.nl
    U = man.volVectorField( man.IOobject( ref.word( "U" ),
                                          ref.fileName( runTime.timeName() ),
                                          mesh,
                                          ref.IOobject.MUST_READ,
                                          ref.IOobject.AUTO_WRITE ), mesh );

    phi = man.compressibleCreatePhi( runTime, mesh, rho, U )
 
    ref.ext_Info() << "Creating turbulence model\n" << ref.nl
    turbulence = man.compressible.turbulenceModel.New( rho, U, phi, pThermo )
    
    return pThermo, p, e, psi, rho, U, phi, turbulence


#--------------------------------------------------------------------------------------
def fun_Ueqn( rho, U, phi, turbulence, p ):
    UEqn = man.fvm.ddt( rho, U ) + man.fvm.div( phi, U ) + man( turbulence.divDevRhoReff( U() ), man.Deps( turbulence, U ) )
  
    ref.solve( UEqn == -man.fvc.grad( p ) );
    return UEqn


#--------------------------------------------------------------------------------------
def fun_eEqn( rho, e, phi, turbulence, p, thermo ):
    ref.solve(  ref.fvm.ddt( rho, e ) + ref.fvm.div( phi , e ) - ref.fvm.laplacian( turbulence.alphaEff(), e )  ==
          - p() * ref.fvc.div( phi() / ref.fvc.interpolate( rho ) ) ) # mixed calculation

    thermo.correct()
    pass
    
    
#--------------------------------------------------------------------------------------
def fun_pEqn( mesh, runTime, thermo, rho, p, psi, U, phi, turbulence, UEqn, cumulativeContErr, nNonOrthCorr ):
      
     rho<<thermo.rho()
 
     rAU = 1.0 / UEqn.A()
  
     U<< rAU * UEqn.H()
     
     phid = ref.surfaceScalarField( ref.word( "phid" ),
                                    ref.fvc.interpolate( psi ) * 
                                             ( ( ref.fvc.interpolate( U ) & mesh.Sf() ) + ref.fvc.ddtPhiCorr( rAU, rho, U, phi ) ) )

     for nonOrth in range( nNonOrthCorr + 1 ):
         pEqn = ref.fvm.ddt( psi, p ) + ref.fvm.div( phid, p ) - ref.fvm.laplacian( rho * rAU, p )
         
         pEqn.solve()
         pass

     if nonOrth == nNonOrthCorr:
         phi<<  pEqn.flux()
         pass
         
     ref.rhoEqn( rho, phi )
  
     cumulativeContErr = ref.compressibleContinuityErrs( rho(), thermo, cumulativeContErr ) #it is necessary to force "mixed calculations" implementation

     U -= rAU * ref.fvc.grad( p )
     U.correctBoundaryConditions()
      
     return cumulativeContErr

    
                      
#--------------------------------------------------------------------------------------
def main_standalone( argc, argv ):

    args = ref.setRootCase( argc, argv )

    runTime = man.createTime( args )

    mesh = man.createMesh( runTime )
    
    pThermo, p, e, psi, rho, U, phi, turbulence = createFields( runTime, mesh )

    cumulativeContErr = ref.initContinuityErrs()
    
    ref.ext_Info() << "\nStarting time loop\n" << ref.nl

    while runTime.loop():
        ref.ext_Info() << "Time = " << runTime.timeName() << ref.nl << ref.nl

        piso, nCorr, nNonOrthCorr, momentumPredictor, transonic, nOuterCorr = ref.readPISOControls( mesh )
        
        CoNum, meanCoNum = ref.compressibleCourantNo( mesh, phi, rho, runTime )
        
        ref.rhoEqn( rho, phi );
        
        UEqn = fun_Ueqn( rho, U, phi, turbulence, p )

        fun_eEqn( rho, e, phi, turbulence, p, pThermo )

        for corr in range( nCorr ) :
            cumulativeContErr = fun_pEqn( mesh, runTime, pThermo, rho, p, psi, U, phi, turbulence, UEqn, cumulativeContErr, nNonOrthCorr )
            pass

        turbulence.correct()

        rho << pThermo.rho()

        runTime.write()

        ref.ext_Info() << "ExecutionTime = " << runTime.elapsedCpuTime() << " s" << "  ClockTime = " << runTime.elapsedClockTime() << " s" << ref.nl << ref.nl
        pass

    ref.ext_Info() << "End\n"

    import os
    return os.EX_OK

    
#--------------------------------------------------------------------------------------
from Foam import FOAM_REF_VERSION
if FOAM_REF_VERSION( ">=", "020000" ):
   if __name__ == "__main__" :
      import sys, os
      argv = sys.argv
      os._exit( main_standalone( len( argv ), argv ) )
      pass
   pass
else:
   from Foam.OpenFOAM import ext_Info
   ext_Info()<< "\nTo use this solver, It is necessary to SWIG OpenFoam1.6 \n "     
    


    
#--------------------------------------------------------------------------------------

