Info<< "Reading field T\n" << endl;

volScalarField T
(
    IOobject
    (
        "T",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);



Info<< "Reading field DT\n" << endl;

volScalarField DT
(
    IOobject
    (
        "DT",
        runTime.name(),
        mesh,
        IOobject::NO_READ, //if not changing with T, then MUST_READ
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field V\n" << endl;

volVectorField V
(
    IOobject
    (
        "V",
        runTime.name(),
        mesh,
        IOobject::NO_READ, //if not changing with T, then MUST_READ
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("zero",dimensionSet(0,1,-1,0,0),vector::zero)//(5 dim) is (7 dim) here, because the last two dims can be omitted.
);

Info<< "Reading field S\n" << endl;

volScalarField S
(
    IOobject
    (
        "S",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


Info<< "Reading field S0\n" << endl;

volScalarField S0
(
    IOobject
    (
        "S0",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE
    ),
    mesh
);


IOdictionary physicalProperties
(
    IOobject
    (
        "physicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//Create a tstop to limit heat source S within time "tstop"

Info<< "Reading stop source time tstop\n" << endl;

//dimensionedScalar tstop
//(
//       physicalProperties.lookup("tstop")
      
//);

//scalar tstop = readScalar(physicalProperties.lookup("tstop"));
scalar Tp = readScalar(physicalProperties.lookup("Tp"));
//scalar Tr = readScalar(physicalProperties.lookup("Tr"));
scalar Rp = readScalar(physicalProperties.lookup("Rp"));
scalar Ld = readScalar(physicalProperties.lookup("Ld"));
scalar tGrowth = readScalar(physicalProperties.lookup("tGrowth"));
scalar tDiffuse = readScalar(physicalProperties.lookup("tDiffuse"));
scalar tStepGrowth = readScalar(physicalProperties.lookup("tStepGrowth"));
scalar tStepDiffuse = readScalar(physicalProperties.lookup("tStepDiffuse"));

#include "createFvModels.H"
#include "createFvConstraints.H"
