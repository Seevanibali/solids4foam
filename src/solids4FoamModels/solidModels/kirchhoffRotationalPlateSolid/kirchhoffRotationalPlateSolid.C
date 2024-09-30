/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifndef OPENFOAM_ORG

#include "kirchhoffRotationalPlateSolid.H"
#include "fvm.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "addToRunTimeSelectionTable.H"
#include "faCFD.H"
#include "linearElastic.H"

//#include "BlockLduSystem.H"
#include "SparseMatrixTemplate.H"
// #include "sparseMatrix.H"
#include "sparseMatrixTools.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(kirchhoffRotationalPlateSolid, 0);
addToRunTimeSelectionTable(solidModel, kirchhoffRotationalPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool kirchhoffRotationalPlateSolid::converged
(
    const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
    const SolverPerformance<scalar>& solverPerfM,
    const SolverPerformance<scalar>& solverPerfw,
#else
    const lduSolverPerformance& solverPerfM,
    const lduSolverPerformance& solverPerfw,
#endif
    const areaScalarField& M,
    const areaScalarField& w
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    const scalar residualM =
        gMax
        (
            (
                mag(M - M.prevIter())
               /max
                (
                    gMax(mag(M - M.oldTime())()), SMALL
                )
            )()
        );

    const scalar residualw =
        gMax
        (
            (
                mag(w - w.prevIter())
               /max
                (
                    gMax(mag(w - w.oldTime())()), SMALL
                )
            )()
        );

    // Calculate material residual
    const scalar materialResidual = mechanical().residual();

    // If one of the residuals has converged to an order of magnitude
    // less than the tolerance then consider the solution converged
    // force at leaast 1 outer iteration and the material law must be converged
    if (iCorr > 1 && materialResidual < materialTol())
    {
        bool convergedM = false;
        bool convergedw = false;

        if
        (
            (
                solverPerfM.initialResidual() < solutionTol()
             && residualM < solutionTol()
            )
         || solverPerfM.initialResidual() < alternativeTol()
         || residualM < alternativeTol()
        )
        {
            convergedM = true;
        }

        if
        (
            (
                solverPerfw.initialResidual() < solutionTol()
             && residualw < solutionTol()
            )
         || solverPerfw.initialResidual() < alternativeTol()
         || residualw < alternativeTol()
        )
        {
            convergedw = true;
        }


        if (convergedM && convergedw)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, res (M & w), relRes (M & w), matRes, iters (M & w)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfM.initialResidual()
            << ", " << solverPerfw.initialResidual()
            << ", " << residualM
            << ", " << residualw
            << ", " << materialResidual
            << ", " << solverPerfM.nIterations()
            << ", " << solverPerfw.nIterations() << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the M-w loop" << endl;
    }

    return converged;
}


const fvPatch& kirchhoffRotationalPlateSolid::areaPatch() const
{
    if (areaPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaPatchID_];
}


const fvPatch& kirchhoffRotationalPlateSolid::areaShadowPatch() const
{
    if (areaShadowPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaShadowPatchID_];
}


void kirchhoffRotationalPlateSolid::calcAreaPatches() const
{
    // Note: face0PatchID may be -1 if this processor has no faces on the
    // finiteArea patch

    // Check that all areaMesh faces map to the same patch

    const polyMesh& pMesh = mesh();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();
    const labelList& faceLabels = aMesh_.faceLabels();
    const label pMeshNFaces = pMesh.nFaces();

    if (faceLabels.size() > 0)
    {
        const label face0ID = faceLabels[0];

        if (face0ID < pMeshNFaces)
        {
            areaPatchID_ = bm.whichPatch(face0ID);

            // Check all faces map to the same fvMesh patch

            forAll(faceLabels, aFaceI)
            {
                const label faceID = faceLabels[aFaceI];

                // Escape if face is beyond active faces, eg belongs to a face
                // zone
                if (faceID < pMeshNFaces)
                {
                    const label curPatchID = bm.whichPatch(face0ID);

                    if (curPatchID != areaPatchID_)
                    {
                        FatalErrorIn
                        (
                            "void kirchhoffRotationalPlateSolid::calcAreaPatches() const"
                        )   << "The finiteArea patch should correspond to a "
                            << "patch on the boundary of the polyMesh!"
                            << abort(FatalError);
                    }
                }
            }
        }
    }


    // We will now check if the polyMesh has the same number of cells as the
    // number of faces on the areaPatch, as we are assuming the polyMesh to be
    // one cell thick
    if (pMesh.nCells() != pMesh.boundaryMesh()[areaPatchID_].size())
    {
        FatalErrorIn
        (
            "void kirchhoffRotationalPlateSolid::calcAreaPatches() const"
        )   << "The solid polyMesh should be one cell thick, where there is "
            << "the same number of cells as the number of faces on the "
            << "areaPatch" << endl
            << "areaPatchID: " << areaPatchID_
            << abort(FatalError);
    }


    // To find the areaShadowPatch, we will ...
    const unallocLabelList& faceCells =
        pMesh.boundaryMesh()[areaPatchID_].faceCells();

    if (faceCells.size())
    {
        const cellList& cells = pMesh.cells();

        const label face0ID = bm[areaPatchID_].start();
        const label cell0ID = faceCells[0];
        const vector& face0N = bm[areaPatchID_].faceNormals()[0];
        const labelList& curCellFaces = cells[cell0ID];

        scalar mostNegativeDotProduct = GREAT;

        forAll(curCellFaces, fI)
        {
            const label curFaceID = curCellFaces[fI];

            if (curFaceID != face0ID)
            {
                if (!pMesh.isInternalFace(curFaceID))
                {
                    const label otherPatchID = bm.whichPatch(curFaceID);
                    const label curLocalFaceID =
                        curFaceID - bm[otherPatchID].start();

                    const vector& curFaceN =
                        bm[otherPatchID].faceNormals()[curLocalFaceID];

                    const scalar dotProduct = face0N & curFaceN;

                    if (dotProduct < mostNegativeDotProduct)
                    {
                        mostNegativeDotProduct = dotProduct;
                        areaShadowPatchID_ = otherPatchID;
                    }
                }
            }
        }
    }


    // Check if the areaPatch and areaShadowPatch have the same number of faces
    if
    (
        pMesh.boundaryMesh()[areaShadowPatchID_].size()
     != pMesh.boundaryMesh()[areaPatchID_].size()
    )
    {
        FatalErrorIn
        (
            "void kirchhoffRotationalPlateSolid::calcAreaPatches() const"
        )   << "The polyMesh should be one cell thick, where there should be "
            << "two patches opposite each other that have the same number of "
            << "faces" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

kirchhoffRotationalPlateSolid::kirchhoffRotationalPlateSolid
(
    Time& runTime,
    const word& region
)
:
    solidModel(typeName, runTime, region),
    aMesh_(mesh()),
    w_
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    wVf_
    (
        IOobject
        (
            "wVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimLength, 0.0)
    ),
    M_
    (
        IOobject
        (
            "M",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    MVf_
    (
        IOobject
        (
            "MVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure*dimArea, 0.0)
    ),
    gradM_(fac::grad(M_)),
    p_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    pVf_
    (
        IOobject
        (
            "pVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", p_.dimensions(), 0.0)
    ),
    // theta_
    // (
    //     IOobject
    //     (
    //         "theta",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::NO_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     aMesh_,
    //     dimensionedVector("zero", dimless, vector::zero)
    // ),
    thetaX_
    (
        IOobject
        (
            "thetaX",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    thetaXVf_
    (
        IOobject
        (
            "thetaXVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    thetaY_
    (
        IOobject
        (
            "thetaY",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    thetaYVf_
    (
        IOobject
        (
            "thetaYVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    // gradTheta_(fac::grad(theta_)),
    gradThetaX_
    (
        IOobject
        (
            "grad(" + thetaX_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimLength/dimArea, vector::zero)
    ),
    gradThetaY_
    (
        IOobject
        (
            "grad(" + thetaY_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimLength/dimArea, vector::zero)
    ),
    gradThXX_
    (
        IOobject
        (
            "grad(" + thetaX_.name() + ")_X",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh_,
        dimensionedScalar("zero", dimLength/dimArea, 0.0)  
    ),
    gradThXY_("grad(" + thetaX_.name() + ")_Y", gradThXX_),
    gradThYX_("grad(" + thetaY_.name() + ")_X", gradThXX_),
    gradThYY_("grad(" + thetaY_.name() + ")_Y", gradThXX_),
    QThetaXX_
    (
        IOobject
        (
            "Q(" + thetaX_.name() + ")_X",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        aMesh_,
        dimensionedScalar("zero", dimForce, 0.0)
    ),
    QThetaXY_("Q(" + thetaX_.name() + ")_Y", QThetaXX_),
    QThetaYX_("Q(" + thetaY_.name() + ")_X", QThetaXX_),
    QThetaYY_("Q(" + thetaY_.name() + ")_Y", QThetaXX_),
    QThetaX_
    (
        IOobject
        (
            "Q(" + thetaX_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_,
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    QThetaY_("Q(" + thetaY_.name() + ")", QThetaX_),
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidModelDict().lookup("plateThickness")),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    torsionalStiffness_("zero", dimPressure*dimVolume, 0.0),
    areaPatchID_(-1),
    areaShadowPatchID_(-1),
    coupled_(solidModelDict().getOrDefault<bool>("coupled", true)),
    debug_(solidModelDict().getOrDefault<bool>("debug", true))
{
    const PtrList<mechanicalLaw>& mechLaws = mechanical();

    // Only the linearElastic mechanicalLaw is allow and one material
    if (mechLaws.size() != 1)
    {
        FatalErrorIn(type() + "::" + type())
            << " can currently only be used with a single material"
            << abort(FatalError);
    }
    else if (!isA<linearElastic>(mechLaws[0]))
    {
        FatalErrorIn(type() + "::" + type())
            << " can only be used with the linearElastic "
            << "mechanicalLaw" << nl
            << abort(FatalError);
    }

    // Cast the mechanical law to a linearElastic mechanicalLaw
    const linearElastic& mech = refCast<const linearElastic>(mechLaws[0]);

    // Set plate properties
    rho_ = mech.rhoScalar();
    E_ = mech.E();
    nu_ = mech.nu();
    bendingStiffness_ = E_*pow(h_, 3)/(12*(1 - pow(nu_, 2)));
    torsionalStiffness_ = E_*pow(h_, 3)/(12*(1 + nu_));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool kirchhoffRotationalPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<scalar> solverPerfM;
        SolverPerformance<scalar> solverPerfw;
        SolverPerformance<scalar>::debug = 0;
#else
        lduSolverPerformance solverPerfM;
        lduSolverPerformance solverPerfw;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the Kirchhoff plate equation for w and M" << endl;

        // Read the time scheme from the test case 
        // WRONG! It is reading from system/fvScehemes, bu we need to read from system/faSchemes
        // const word ddtSchemeName(mesh().ddtSchemes());
        // Read ddtScheme name from system/faSchemes
        const bool ddtScheme(aMesh_.ddtSchemes().found("ddt(w)"));
        const word ddtSchemeName
        (
            ddtScheme
            ?
            (aMesh_.ddtSchemes().lookup("ddt(w)"))
            :
            (aMesh_.ddtSchemes().lookup("default"))
        );

        // Read d2dt2Scheme name from system/faSchemes
        const bool d2dt2Scheme(aMesh_.d2dt2Schemes().found("d2dt2(w)"));
        const word d2dt2SchemeName
        (
            d2dt2Scheme
            ?
            (aMesh_.d2dt2Schemes().lookup("d2dt2(w)"))
            :
            (aMesh_.d2dt2Schemes().lookup("default"))
        );

        WarningIn("evolve() function in kirchhoff plate rotation-free solid")<< nl 
            << "d2dt2Scheme in system/faSchemes cannot take steadyState as a valid keyword!" << nl 
            << "If you want plate-case to be solved for steady state condition, "
            << "set ddtScheme to be steadyState instead!! " << endl;

        // Algorithm for rotational Kirchhoff plate formulation
        // The M equation is:
        //     rho*h*fac::d2dt2(w) = fam::laplacian(M) + p
        // Have to add more here
        // where
        // M is the moment sum
        // w is the transvere (out of plane) displacement
        // rho is the density
        // h is the plate thickness
        // p is the net transverse pressure
        // D is the bending stiffness
     
        
        // Store the previous iteration values for computing source vector
        gradM_.storePrevIter();
        thetaX_.storePrevIter();
        thetaY_.storePrevIter();
        gradThetaX_.storePrevIter();
        gradThetaY_.storePrevIter();

        Info<< "\nUsing segregated approach to solve for w & M eqns " 
            << "separately and iteratively update them!" << endl;

        do
        { 
            // Store fields for under-relaxation and residual calculation
            M_.storePrevIter();

            // Solve M equation
            // Also, "==" complains so we will move all terms to left
            faScalarMatrix MEqn
            (
                - fam::laplacian(M_) - p_
            );

            // d2dt2 can only take Euler as keyword, but if the user wants it to
            // be steadyState, it cannot happen. Hence check for ddtScheme 
            // and add inertial terms for not steady state!!
            if(ddtSchemeName != "steadyState")
            {
                MEqn += rho_*h_*(fac::d2dt2(w_));
            }

            // Relax the linear system
            MEqn.relax();

            // Solve the linear system
            solverPerfM = MEqn.solve();

            // Relax the field
            M_.relax();

            // Update the gradient of moment sum
            gradM_ = fac::grad(M_);

            // Interpolate the gradient values to edge centres
            const edgeVectorField gradMEdge(fac::interpolate(gradM_));

            // Access to area mesh data 
            const areaVectorField& cellCentres(aMesh_.areaCentres());
            const label nCells(aMesh_.faceCells().size());
            const labelList& edgeOwn(aMesh_.edgeOwner());
            const labelList& edgeNei(aMesh_.edgeNeighbour());
            const edgeVectorField& edgeCentres(aMesh_.edgeCentres());
            const edgeVectorField& Le(aMesh_.Le());

            scalarField& gradThXXI(gradThXX_.ref());
            scalarField& gradThXYI(gradThXY_.ref());
            scalarField& gradThYXI(gradThYX_.ref());
            scalarField& gradThYYI(gradThYY_.ref());

            forAll(gradThXXI, cellI)
            {
                gradThXXI[cellI] = gradThetaX_.internalField()[cellI].component(0);
                gradThXYI[cellI] = gradThetaX_.internalField()[cellI].component(1);
                gradThYXI[cellI] = gradThetaY_.internalField()[cellI].component(0);
                gradThYYI[cellI] = gradThetaY_.internalField()[cellI].component(1);
            }
        
            forAll(gradThXX_.boundaryField(), patchI) 
            {
                scalarField& pGradThXX(gradThXX_.boundaryFieldRef()[patchI]);
                scalarField& pGradThXY(gradThXY_.boundaryFieldRef()[patchI]);
                scalarField& pGradThYX(gradThYX_.boundaryFieldRef()[patchI]);
                scalarField& pGradThYY(gradThYY_.boundaryFieldRef()[patchI]);

                forAll (pGradThXX, faceI) 
                {
                    pGradThXX[faceI] = gradThetaX_.boundaryField()[patchI][faceI].component(0);
                    pGradThXY[faceI] = gradThetaX_.boundaryField()[patchI][faceI].component(1);
                    pGradThYX[faceI] = gradThetaY_.boundaryField()[patchI][faceI].component(0);
                    pGradThYY[faceI] = gradThetaY_.boundaryField()[patchI][faceI].component(1);
                
                }
            }

            QThetaXX_ =
            (
                (bendingStiffness_ - torsionalStiffness_)*gradThXX_
                + (bendingStiffness_*nu_*gradThYY_)
            );
            QThetaXY_ = 
            (
                torsionalStiffness_*gradThYX_ 
            );
            QThetaYX_ = 
            ( 
                torsionalStiffness_*gradThXY_
            );
            QThetaYY_ =
            (
                (bendingStiffness_ - torsionalStiffness_)*gradThYY_
                + (bendingStiffness_*nu_*gradThXX_)
            );
            
            vectorField& QThetaXI = QThetaX_.ref();
            vectorField& QThetaYI = QThetaY_.ref();
            forAll(QThetaXI, cellI)
            {
                QThetaXI[cellI].component(0) = QThetaXX_.internalField()[cellI];
                QThetaXI[cellI].component(1) = QThetaXY_.internalField()[cellI];
                QThetaYI[cellI].component(0) = QThetaYX_.internalField()[cellI];
                QThetaYI[cellI].component(1) = QThetaYY_.internalField()[cellI];
            }
        
            forAll(QThetaX_.boundaryField(), patchI) 
            {
                vectorField& pQThetaX(QThetaX_.boundaryFieldRef()[patchI]);
                vectorField& pQThetaY(QThetaY_.boundaryFieldRef()[patchI]);

                forAll (pQThetaX, faceI) 
                {
                    pQThetaX[faceI].component(0) = QThetaXX_.boundaryField()[patchI][faceI];
                    pQThetaX[faceI].component(1) = QThetaXY_.boundaryField()[patchI][faceI];
                    pQThetaY[faceI].component(0) = QThetaYX_.boundaryField()[patchI][faceI];
                    pQThetaY[faceI].component(1) = QThetaYY_.boundaryField()[patchI][faceI];
                }
            }

            //- Initialise another edgeVectorField to store edgeCentres minus cellCentres
            //- For every edge, subtract position vector of the owner of the edge 
            //- from the edgeCentre position vector because origin assumed at centre of CV  
            edgeVectorField rVec
            (
                IOobject
                (
                    "rVec",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                aMesh_,
                dimensionedVector("zero", dimLength, vector(0, 0, 0)) 
            );

            forAll(edgeCentres.internalField(), edgeI)
            {
                rVec[edgeI] = edgeCentres[edgeI] - cellCentres[edgeOwn[edgeI]];
            }

            forAll(edgeCentres.boundaryField(), patchI)
            {
                faePatchField<vector>& pRVec = rVec.boundaryFieldRef()[patchI];
                forAll(edgeCentres.boundaryField()[patchI], edgeI)
                {
                    const label bI = aMesh_.boundary()[patchI].edgeFaces()[edgeI];
                    pRVec[edgeI] = edgeCentres.boundaryField()[patchI][edgeI] 
                        - cellCentres.internalField()[bI];
                }
            }

            // Solve thetaX equation 
            faScalarMatrix thetaXEqn
            (
                fam::laplacian(torsionalStiffness_, thetaX_)
                    - fac::edgeIntegrate(Le & (mag(rVec.component(0))*gradMEdge))
                    + fac::div(QThetaX_)
            );

        
            Info<< "*-----------------------------------*" << nl
            << "theta X source " << thetaXEqn.source() << nl 
            << "QThetaX div term " << fac::div(QThetaX_) << nl
            << "mesh Le " << aMesh_.Le() << nl 
            // << "rVec " << rVec << nl 
            << "*---------------------------------------*" << nl
            << endl;

            // Solve the linear system
            thetaXEqn.solve();
            thetaX_.relax();

            // gradThetaX_ = fac::grad(thetaX_);

            // Solve thetaY equation
            faScalarMatrix thetaYEqn
            (
                fam::laplacian(torsionalStiffness_, thetaY_) 
                - fac::edgeIntegrate(Le & (mag(rVec.component(1))*gradMEdge))
                + fac::div(QThetaY_) 
            );

            // Solve the linear system
            thetaYEqn.solve();
            thetaY_.relax();

            // For now, update gradThetaX and gradThetaY after solving
            // both thetaX and thetaY; Dont know exactly if this is correct
            // or better
            gradThetaX_ = fac::grad(thetaX_);
            gradThetaY_ = fac::grad(thetaY_);

            // Store fields for under-relaxation and residual calculation
            w_.storePrevIter();
            
            // Solve w equation
            // faScalarMatrix wEqn
            // (
            //     fam::laplacian(bendingStiffness_, w_) + M_
            // );

            /*---------------------------------------------------------------*/
            // const edgeScalarField& thetaXEdge(fac::interpolate(thetaX_));
            // NOTE: The above code compiles but does not display desired results at runtime.
            // Because you are trying to catch the results of a temporary object which deletes
            // itself after creation. Hence Result is zero.

            const edgeScalarField thetaXEdge(fac::interpolate(thetaX_));
            const scalarField& thetaXEdgeI(thetaXEdge.internalField());
            // scalarField thetaXEdgeI(aMesh_.edges().size(), 1.0);
            const edgeScalarField thetaYEdge(fac::interpolate(thetaY_));
            const scalarField& thetaYEdgeI(thetaYEdge.internalField());
            // scalarField thetaYEdgeI(aMesh_.edges().size(), 1.0);

            // Initialise matrix (For w eqn per cell)
            SparseMatrixTemplate<scalar> matrix(nCells);

            // Initialse source vector
            scalarField source(nCells, 0.0);

            // Initioalise solution w field
            scalarField solveW(nCells, 0.0);

            forAll(aMesh_.faces(),faceI)
            {
                forAll(aMesh_.edges(), edgeI)
                {
                    if(aMesh_.isInternalEdge(edgeI))
                    {
                        if (edgeOwn[edgeI] == faceI)
                        {
                            matrix(faceI, edgeNei[edgeI]) = -1;
                            matrix(edgeNei[edgeI], faceI) = -1;
                        }
                        if ((edgeOwn[edgeI] == faceI) || (edgeNei[edgeI] == faceI))
                        {
                            vector deltaR(vector::zero);
                            if(edgeOwn[edgeI] == faceI)
                            {
                                deltaR = cellCentres[edgeNei[edgeI]] 
                                    - cellCentres[edgeOwn[edgeI]];
                            }
                            else
                            {
                                deltaR = cellCentres[edgeOwn[edgeI]] 
                                    - cellCentres[edgeNei[edgeI]];
                            }
                            source[faceI] += deltaR.component(0)*thetaXEdgeI[edgeI]
                                + deltaR.component(1)*thetaYEdgeI[edgeI];
                        }
                    }
                }
            }

            // Populating the diagonal entries of the matrix - negSumDiag 
            forAll (aMesh_.faceCells(), rowI)
            {
                label count = 0;
                forAll (aMesh_.faceCells(), colI)
                {
                    if (matrix(rowI,colI) == -1)
                    {
                        count += 1;
                    }
                }
                matrix(rowI, rowI) = count;
            }
            if(debug_)
            {
                Info<< "Before boundary data " << endl;   
                Info<< "matrix " << nl << matrix.data()
                    << "source " << nl << source << nl 
                    << "****************************" << endl; 
            }
            //- Contribution of the boundary edges to the w eqn
            forAll(aMesh_.boundary(), patchI)
            {
                const vectorField pDelta(aMesh_.boundary()[patchI].delta());
                
                // This seems to be like a dynamic list of integers
                // const UList<int>& bouEdgeLabels(aMesh_.boundary().edgeLabels()[patchI]);
                const scalarField& pThetaX = thetaX_.boundaryField()[patchI];
                const scalarField& pThetaY = thetaY_.boundaryField()[patchI];
                const scalarField& pW = w_.boundaryField()[patchI];

                if(debug_)
                {
                    Info<< "---------------------------------" << nl 
                        <<  "patchI " << aMesh_.boundary()[patchI].name() << nl
                        << "delta " << tab << pDelta << nl 
                        << "pThetaX " << tab << pThetaX << nl 
                        << "pThetaY " << tab << pThetaY << nl 
                        << "pW " << tab << pW <<  endl;
                }
                forAll(aMesh_.boundary()[patchI], faceI)
                {
                    const label bI = aMesh_.boundary()[patchI].edgeFaces()[faceI];
                    matrix(bI,bI) += 1;

                    // Info<< "source term: " << tab <<  pDelta[faceI].component(0)*pThetaX[faceI]
                    //      + pDelta[faceI].component(1)*pThetaY[faceI] + pW[faceI] << tab 
                    //     << " added at cell location " << tab << bI << endl;


                    source[bI] += pDelta[faceI].component(0)*pThetaX[faceI]
                            + pDelta[faceI].component(1)*pThetaY[faceI] + pW[faceI];

                }

            }
            
            if(debug_)
            {
                Info<< "****************************" << endl; 
                Info<< "After adding  boundary data " << endl;
                Info<< "matrix " << nl << matrix.data()
                    << "source " << nl << source << nl 
                    << "****************************" << endl; 
            }
            // Relax the linear system
            // wEqn.relax();

            // Solve the linear system
            // solverPerfw = wEqn.solve();

            // Using Eigen SparseLU direct solver
            sparseMatrixTools::solveLinearSystemEigen
            (
                matrix, source, solveW, false, debug
            );

            // Retrieve solution
            for(label i = 0; i < nCells; ++i)
            {
                w_[i] = solveW[i];
            }

            // Relax the field
            w_.relax();

            // Update the angle of rotation
            // theta_ = -fac::grad(w_);
            // thetaX_(theta.component(0));
            // thetaY_(theta.component(1));

            // Update the gradient of rotation field, used for non-orthogonal
            // correction in clamped boundary conditions
            // gradTheta_ = fac::grad(theta_);
            

            // Info<< "theta " << theta_ << nl << endl;
            // << "gradTheta " << gradTheta_ << endl;         
        }
        while
        (
            // !converged(iCorr, solverPerfM, solverPerfw, M_, w_)
            // && 
            ++iCorr < nCorr()
        );

    
        
        // Map area fields to vol fields
        mapAreaFieldToSingleLayerVolumeField(M_, MVf_);
        mapAreaFieldToSingleLayerVolumeField(w_, wVf_);
        mapAreaFieldToSingleLayerVolumeField(thetaX_, thetaXVf_);
        mapAreaFieldToSingleLayerVolumeField(thetaY_, thetaYVf_);
        mapAreaFieldToSingleLayerVolumeField(p_, pVf_);
        {
            const areaVectorField Ds(w_*aMesh_.faceAreaNormals());
            mapAreaFieldToSingleLayerVolumeField(Ds, D());
        }

        // Interpolate cell displacements to vertices
        mechanical().interpolate(D(), pointD());

        // Increment of displacement
        DD() = D() - D().oldTime();

        // Increment of point displacement
        pointDD() = pointD() - pointD().oldTime();

        // Velocity
        U() = fvc::ddt(D());
    }
    while (mesh().update());

    // Disable writing of the stress fields, as it is not calculated
    sigma().writeOpt() = IOobject::NO_WRITE;

    return true;
}


tmp<vectorField> kirchhoffRotationalPlateSolid::tractionBoundarySnGrad
(
    const vectorField& traction,
    const scalarField& pressure,
    const fvPatch& patch
) const
{
    notImplemented(type() + "::tractionBoundarySnGrad(...)");

    // Keep compiler happy
    return tmp<vectorField>();
}


void kirchhoffRotationalPlateSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Map global field to patch field
    const vectorField patchTraction
    (
        globalPatches()[interfaceI].globalFaceToPatch(faceZoneTraction)
    );

    // Take normal component of the traction field
    // Note: p is the net pressure on the plate (from both sides)
#ifdef OPENFOAM_NOT_EXTEND
    p_.primitiveFieldRef() =
        aMesh_.faceAreaNormals().primitiveField() & patchTraction;
#else
    p_.internalField() =
        aMesh_.faceAreaNormals().internalField() & patchTraction;
#endif
}


void kirchhoffRotationalPlateSolid::writeFields(const Time& runTime)
{
    // Do not call solidModel::writeFields() as we do not want to write the
    // stress and strain fields

    physicsModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif // OPENFOAM_ORG

// ************************************************************************* //
