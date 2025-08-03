/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "dynamicRekEpsilon.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class BasicTurbulenceModel>
void dynamicRekEpsilon<BasicTurbulenceModel>::correctNut()
{
    this->nut_ = Cmu_*sqr(k_)/epsilon_;
    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> dynamicRekEpsilon<BasicTurbulenceModel>::kSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            k_,
            dimVolume*this->rho_.dimensions()*k_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> dynamicRekEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
dynamicRekEpsilon<BasicTurbulenceModel>::dynamicRekEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    Re_  ("Re",10000.0),
    Cmu_ ("Cmu", 0.09),
    C1_ ("C1", 1.44),
    C2_ ("C2", 1.92),
    C3_ ("C3", 0.0),
    sigmak_ ("sigmak", 1.0),
    sigmaEps_ ("sigmaEps", 1.3),
    k_
    (
        IOobject
        (
            IOobject::groupName("k", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", alphaRhoPhi.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    )
{
    if (!this->coeffDict().found("Re"))
    { 
//        InfoErr
          FatalErrorInFunction
            << "No reference Reynolds number specified in dictionary." << nl
            << "Please add the Re field in turbulence dictionary." << nl
            << exit(FatalError);
    }
    else
    {
        Re_.read(this->coeffDict());
        if ((Re_.value() < 5000))
        {
            Re_ = 5000;
            Info << nl << "Warning: Reynolds number out of bounds. The Reynolds number is set to 5000" << endl;
        }
        if ((Re_.value() > 86000))
        {
            Re_ = 86000;
            Info << nl << "Warning: Reynolds number out of bounds. The Reynolds number is set to 86000" << endl;
        }
    }
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    read();
    
    Info << endl << "dynamicRekEpsilon Coeffs:" << endl << "Reynolds number = " << Re_.value() << endl 
                                    << "Cmu = " << Cmu_.value() << endl 
                                    << "C1 = " << C1_.value() << endl  
                                    << "C2 = " << C2_.value() << endl  
                                    << "C3 = " << C3_.value() << endl  
                                    << "Sigmak = " << sigmak_.value() << endl 
                                    << "SigmaEps = " << sigmaEps_.value() << endl << endl;
//    if (type == typeName)
//    {
//        this->printCoeffs(type);
//    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool dynamicRekEpsilon<BasicTurbulenceModel>::read()
{
    if (Re_.value() < 23000) 
    { 
         Cmu_ = -2.27066e-14*sqr(Re_)*Re_ + 8.24265e-10*sqr(Re_) - 7.82289e-06*Re_ + 0.112089;
         C1_ = 1.78878e-13*sqr(Re_)*Re_ - 9.00433e-09*sqr(Re_) + 0.000138488*Re_ + 0.802416;
         C2_ = 4.63161e-13*sqr(Re_)*Re_ - 2.21972e-08*sqr(Re_) + 0.000314514*Re_ + 0.756003;
         sigmak_ = -2.79032e-13*sqr(Re_)*Re_ + 1.48796e-08*sqr(Re_) - 0.000234409*Re_+ 1.57649;
         sigmaEps_ =  -1.46374e-13*sqr(Re_)*Re_ + 8.38289e-09*sqr(Re_) - 0.000168028*Re_ + 2.86225;
    }
    else
    {
         Cmu_ = -6.83529e-12*sqr(Re_) + 9.296e-07*Re_ + 0.0741632;
         C1_ = -1.45115e-07*Re_ +  1.402;
         C2_ = 2.88315e-11*sqr(Re_) - 3.92108e-06*Re_ + 1.95772;
         sigmak_ = 1.32574e-07*Re_ + 0.660281;
         sigmaEps_ = 1.24209e-10*sqr(Re_) - 1.68924e-05*Re_ + 1.97404;
    }
    C3_ = 0.0;
    return true;
}


template<class BasicTurbulenceModel>
void dynamicRekEpsilon<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    const volScalarField::Internal divU
    (
        fvc::div(fvc::absolute(this->phi(), U))().v()
    );

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        IOobject::scopedName(this->type(), "GbyNu"),
        tgradU().v() && devTwoSymm(tgradU().v())
    );
    const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    tgradU.clear();

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        C1_*alpha()*rho()*GbyNu*Cmu_*k_()
      - fvm::SuSp(((2.0/3.0)*C1_ - C3_)*alpha()*rho()*divU, epsilon_)
      - fvm::Sp(C2_*alpha()*rho()*epsilon_()/k_(), epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    epsEqn.ref().relax();
    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha()*rho()*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
