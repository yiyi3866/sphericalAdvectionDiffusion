/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    laplacianFoam

Description
    Solves a simple Laplace equation, e.g. for thermal diffusion in a solid.

\*---------------------------------------------------------------------------*/
#include "fvMesh.H"
#include "argList.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "simpleControl.H"

#include "fvmLaplacian.H"
#include <iomanip>  // Required for std::setprecision


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"
   

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;
   
    scalar yearToSec = 31556952; // Conversion factor
   
    scalar tGrowth_sec = tGrowth * yearToSec;
    scalar tDiffuse_sec = tDiffuse * yearToSec;
    runTime.setEndTime(tDiffuse_sec);
   
   
    //scalar deltaT_growth_sec = 10 * yearToSec; //10 years
    //scalar deltaT_diffuse_sec = 40 * yearToSec; //40 years
   
    //scalar writeInterval_growth_sec = 100 * yearToSec; //100 years
    //scalar writeInterval_diffuse_sec = 10000 * yearToSec; //10000 years
       
    scalar deltaT_growth_sec = tStepGrowth*yearToSec;//tGrowth_sec / 10000;//10000
    scalar deltaT_diffuse_sec = tStepDiffuse*yearToSec;//tGrowth_sec / 50000;//250
   
    scalar writeInterval_growth_sec = tGrowth_sec / 100;
    scalar writeInterval_diffuse_sec = tDiffuse_sec / 200;

    //while (simple.loop(runTime))ra
    while (runTime.loop())
    {
        runTime++;
        Info<< "Time = " << runTime.userTimeName() << nl << endl;

        fvModels.correct();
       
        const fvMesh& mesh = T.mesh();
       
        scalar a = Rp * Foam::pow(tGrowth_sec,-1.0/3.0);


    if (runTime.value() <= tGrowth_sec)
    {
        runTime.setDeltaT(deltaT_growth_sec);
        runTime.setWriteInterval(writeInterval_growth_sec);
    }
    else
    {
        runTime.setDeltaT(deltaT_diffuse_sec);
        runTime.setWriteInterval(writeInterval_diffuse_sec);
    }

        // V
        if (runTime.value() <= tGrowth_sec) {
        forAll(mesh.cells(), cellI) {
            scalar z = mesh.C()[cellI].z(); // Get the z-coordinate of the cell center
            scalar t = runTime.value(); // Get the current time
            scalar Vz; // Variable for the z-component of velocity
            scalar Vd = Foam::min(1.0/3.0*a*Foam::pow(t,-2.0/3.0), 4e-08);
            scalar rp = a*Foam::pow(t,1.0/3.0); //instaneous pluton radius
            if (z <= rp) {
                Vz = Vd;
            } else if (z > rp + Ld) {
                Vz = 0.0;
            } else {
                //linear decrease to 0:
                //Vz = Vd / Ld * (Ld - z + a*Foam::pow(t,1.0/3.0));
                Vz = Vd * (rp*rp) / (z*z);//to be checked
            }
            //Vz = 1;
       
            // Accessing internal field of V and modifying the z-component
            V[cellI] = vector(0.0, 0.0, Vz); // Vx = Vy = 0
        }
        } else  {
        forAll(mesh.cells(), cellI) {
            V[cellI] = vector(0.0, 0.0, 0.0);
            }
        }

        V.correctBoundaryConditions();

        Info << "Min/Max V: " << gMin(V) << " / " << gMax(V) << endl; //check if V exists before phi
       
        //Radius at cell centers (distances from the pluton center)
        volScalarField r = mag(mesh.C());
        //Jacobian J = r^2
        volScalarField J = r*r;
        //Face-interpolated J
        surfaceScalarField Jf = fvc::interpolate(J);
        //Face flux: phir = Jf * (v_face.*Sf)
        surfaceScalarField phif = Jf * (fvc::interpolate(V) & mesh.Sf());

        //If needed, print phif inside your solver for debugging:
        Info << "Min/Max phi: " << gMin(phif) << " / " << gMax(phif) << endl;
       
        // DT

        forAll(mesh.cells(), cellI)
        {

            // Access temperature field value at cell center
            double temp = T[cellI];
            double DT_value;
            if (temp < 846.0) {
            DT_value= (567.3 / temp - 0.062)*0.000001;
              } else {
               DT_value = (0.732 - 0.000135 * temp)*0.000001;
               }
             DT[cellI] = DT_value;
           
        }
        DT.correctBoundaryConditions();
       


        while (simple.correctNonOrthogonal()) {

        //volScalarField sourceField = (runTime.value() > tstop.value()) ? S0 : S;

        //Diffusion: compute face diffusive fluxes
        //face gradient of T in  radial direction: use cell gradient then interpolate to faces
        volVectorField gradT = fvc::grad(T);
        surfaceVectorField gradTfn = fvc::interpolate(gradT);//face gradients (vector)
       
        //For radial 1D aligned with z axis, radial direction is unit vector from origin to face center
        surfaceVectorField faceTtr =  mesh.Cf();
        
        surfaceScalarField rf = mag(faceTtr); //face radius

        surfaceVectorField er = faceTtr / rf; //unit radial vector on faces (carefully handle rf~0)
       
        //Radial derivative at face = gradTtn dot er
        surfaceScalarField dTdr_f = (gradTfn & er);
       
        //Interpolate DT to faces
        surfaceScalarField DTf = fvc::interpolate(DT);
       
      
        //Mesh is 1D in z
        //We can use diffFluxScalar = -Jf * DTf * dTdr_f * mag(Sf) * sign(er.*Sf). For robust code, do:
        surfaceScalarField erDotSf = (er & mesh.Sf());
        surfaceScalarField diffFluxScalar = -Jf * DTf * dTdr_f * erDotSf;
       
        //Add diffusion divergence (note sign)
        //TEqn -= fvc::div(diffFluxScalar);
        Info << "DT range: " << min(DT).value() << " to " << max(DT).value() << endl;
        Info << "T range: " << min(T).value() << " to " << max(T).value() << endl;
        fvScalarMatrix TEqn 
        (
            fvm::ddt(J,T) 
            + fvc::div(phif,T) //Advection
            + fvc::div(diffFluxScalar) //Diffusion
            ==
            J * S //RHS source
        );
        

        fvConstraints.constrain(TEqn);
        TEqn.solve();
        fvConstraints.constrain(T);
    }
        Info << "DT range: " << min(DT).value() << " to " << max(DT).value() << endl;
        Info << "T range: " << min(T).value() << " to " << max(T).value() << endl;

       
        //set inner boundary and pluton temperature only before tGrowth
        if (runTime.value() <= tGrowth_sec) {
        scalar zb = a * Foam::pow(runTime.value(),1.0/3.0);
        Info << "Current time: " << runTime.value() << ", zb: " << zb << nl;//check lines
        forAll(mesh.cells(), cellI) {
            scalar z = mesh.C()[cellI].z(); // Get the z-coordinate of the cell center
                if (z < zb) {
                        T[cellI] = Tp; // K
                }
            }
        }
        T.correctBoundaryConditions();

   
        #include "write.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
