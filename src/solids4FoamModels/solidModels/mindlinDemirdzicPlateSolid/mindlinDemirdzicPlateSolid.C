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

#include "mindlinDemirdzicPlateSolid.H"
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

defineTypeNameAndDebug(mindlinDemirdzicPlateSolid, 0);
addToRunTimeSelectionTable(solidModel, mindlinDemirdzicPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool mindlinDemirdzicPlateSolid::converged
(
    const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
    const SolverPerformance<scalar>& solverPerfw,
    const SolverPerformance<scalar>& solverPerfThetaX,
    const SolverPerformance<scalar>& solverPerfThetaY,
#else
    const lduSolverPerformance& solverPerfw,
    const lduSolverPerformance& solverPerfThetaX,
    const lduSolverPerformance& solverPerfThetaY,
#endif
    const areaScalarField& w,
    const areaScalarField& thetaX,
    const areaScalarField& thetaY
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    const scalar resThetaX =
        gMax
        (
            (
                mag(thetaX - thetaX.prevIter())
               /max
                (
                    gMax(mag(thetaX - thetaX.oldTime())()), SMALL
                )
            )()
        );

    const scalar resThetaY =
    gMax
    (
        (
            mag(thetaY - thetaY.prevIter())
            /max
            (
                gMax(mag(thetaY - thetaY.oldTime())()), SMALL
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

        bool convergedw = false;
        bool convergedThetaX = false;
        bool convergedThetaY = false;

        if
        (
            (
                // solverPerfw.initialResidual() < solutionTol()
            //  && 
                residualw < solutionTol()
            )
        //  || solverPerfw.initialResidual() < alternativeTol()
        //  || residualw < alternativeTol()
        )
        {
            convergedw = true;
        }

        if
        (
            (
                // solverPerfThetaX.initialResidual() < solutionTol()
            //  && 
                resThetaX < solutionTol()
            )
        //  || solverPerfThetaX.initialResidual() < alternativeTol()
        //  || resThetaX < alternativeTol()
        )
        {
            convergedThetaX = true;
        }

        if
        (
            (
                // solverPerfThetaY.initialResidual() < solutionTol()
            //  && 
                resThetaY < solutionTol()
            )
        //  || solverPerfThetaY.initialResidual() < alternativeTol()
        //  || resThetaY < alternativeTol()
        )
        {
            convergedThetaY = true;
        }

        if (convergedw && (convergedThetaX && convergedThetaY))
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // const bool var = true;
    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, solnRes (w, thetaX, thetaY), relRes (w, thetaX, thetaY), iters (w, thetaX, thetaY)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfw.initialResidual()
            << ", " << solverPerfThetaX.initialResidual()
            << ", " << solverPerfThetaY.initialResidual()
            << ", " << tab << residualw
            << ", " << resThetaX
            << ", " << resThetaY
            << ", " << tab << solverPerfw.nIterations()
            << ", " << solverPerfThetaX.nIterations()
            << ", " << solverPerfThetaY.nIterations()
            << endl;

        if (converged)
        {
            Info<< endl;
        }
    }
    else if (iCorr == nCorr() - 1)
    {
        maxIterReached()++;
        Warning
            << "Max iterations reached within the w-theta loop" << endl;
    }

    return converged;
}


const fvPatch& mindlinDemirdzicPlateSolid::areaPatch() const
{
    if (areaPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaPatchID_];
}


const fvPatch& mindlinDemirdzicPlateSolid::areaShadowPatch() const
{
    if (areaShadowPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaShadowPatchID_];
}


void mindlinDemirdzicPlateSolid::calcAreaPatches() const
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
                            "void mindlinDemirdzicPlateSolid::calcAreaPatches() const"
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
            "void mindlinDemirdzicPlateSolid::calcAreaPatches() const"
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
            "void mindlinDemirdzicPlateSolid::calcAreaPatches() const"
        )   << "The polyMesh should be one cell thick, where there should be "
            << "two patches opposite each other that have the same number of "
            << "faces" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mindlinDemirdzicPlateSolid::mindlinDemirdzicPlateSolid
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
    gradW_
    (
        IOobject
        (
            "grad(" + w_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,   
            IOobject::AUTO_WRITE 
        ),
        aMesh_,
        dimensionedVector("zero", dimless, vector::zero)
    ),
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
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidModelDict().lookup("plateThickness")),
    shearCorrectionFactor_(solidModelDict().getOrDefault<scalar>("shearCorrectionFactor", 0.8333)),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    // torsionalStiffness_("zero", dimPressure*dimVolume, 0.0),
    shearStrainStiffness_("zero",dimForce/dimLength, 0.0),
    areaPatchID_(-1),
    areaShadowPatchID_(-1),
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
    // torsionalStiffness_ = 0.5*(1 - nu_)*bendingStiffness_;
    shearStrainStiffness_ = shearCorrectionFactor_*0.5*E_*h_/(1 + nu_);

    Info<< "Plate mechanical properties\n"
        << "Bending Stiffness\n" << bendingStiffness_ << "\n"
        // << "\nTorsional Stiffness\n" << torsionalStiffness_ << "\n"
        << "\nShear strain stiffness\n" << shearStrainStiffness_
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool mindlinDemirdzicPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Lookup flag for using a compact stencil for the edge normal gradients
    const Switch compactEdgeNormalGrad
    (
        solidModelDict().lookup("compactEdgeNormalGrad")
    );

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        // SolverPerformance<vector> solverPerfShearForce;
        SolverPerformance<scalar> solverPerfw;
        SolverPerformance<scalar> solverPerfThetaX;
        SolverPerformance<scalar> solverPerfThetaY;
        SolverPerformance<scalar>::debug = 0;
#else
        // lduSolverPerformance solverPerfShearForce;
        lduSolverPerformance solverPerfw;
        lduSolverPerformance solverPerfThetaX;
        lduSolverPerformance solverPerfThetaY;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the Mindlin (thick plates) equation for primary variables w, thetaX, thetaY - Demirdzic's Approach" << endl;

        Info<< "\nUsing segregated approach to solve for w, thetaX, and thetaY eqns "
            << "separately and iteratively update them!" << endl;

        // Philip testing
        const scalar alphaW(readScalar(solidModelDict().lookup("alphaW")));
        const scalar alphaTheta(readScalar(solidModelDict().lookup("alphaTheta")));

        // Note: To get in-plane normal unit vectors to an edge, aMesh_.Le()
        // can be used with unit norm
        // Do not use aMesh.unitLe() member from faMesh, since the
        // boundary values of aMesh.unitLe() are set to
        // "calculated (0 0 0)" which is not correct!!
        // Constructing the unit edgeBiNormal instead
        // For 2-D meshes, aMesh_.edgeNormals() gives unit vector in z-direction
        edgeVectorField edgeBiNormal
        (
            IOobject
            (
                "edgeBiNormal",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh_.Le()/aMesh_.magLe()
        );

        // Mesh information required
        const labelList& own(aMesh_.owner());
        const labelList& nei(aMesh_.neighbour());
        // const DimensionedField<scalar, areaMesh>& Sf(aMesh_.S());
        const edgeScalarField& le(aMesh_.magLe());
        const faBoundaryMesh& faBouMesh(aMesh_.boundary());
        const edgeVectorField& edgeCentres(aMesh_.edgeCentres());
        const areaVectorField& cellCentres(aMesh_.areaCentres());
        // const labelList& edgeOwn(aMesh_.edgeOwner());
        // const labelList& edgeNei(aMesh_.edgeNeighbour());

        // NOTE!! - This theta construction was outside the do-loop in
        // previous commit. So, the theta contribution was not getting
        // added to w equation. But now that I add theta contribution
            // the solver is diverging. Do I need more stabilisation terms?
            // CHECK WITH IVAN
        
        do
        {

        // Constructing the theta vector from components
        areaVectorField theta
        (
            IOobject
            (
                "theta",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            // Combine components using unit vectors
            thetaX_ * vector(1, 0, 0)
            + thetaY_ * vector(0, 1, 0)
            + dimensionedScalar("zero", thetaX_.dimensions(), 0.0) * vector(0, 0, 1) 
        );

        // Theta vector at edge centres
        const edgeVectorField thetaEdge(fac::interpolate(theta));
            // Store the previous iteration values for computing source vector
            // The storePrevIter values are also brought inside this do-loop
            // Should not make much difference.

            w_.storePrevIter();
            thetaX_.storePrevIter();
            thetaY_.storePrevIter();
            gradW_.storePrevIter();
            gradThetaX_.storePrevIter();
            gradThetaY_.storePrevIter();

            // Info<< "thetaEdge " << thetaEdge << endl;

            // Solve w equation
            // Also, "==" complains so we will move all terms to left
            // QUESTION - Is fac::div(shearStrainStiffness_*theta)
            // equivalent to physically looping over edges and putting
            // theta contributions interpolated at the edges in the wEqn?
            // I think it is coorect because that is how fac::div
            // is calculated, but need to CHECK WITH IVAN.
            faScalarMatrix wEqn
            (
                fam::laplacian(shearStrainStiffness_, w_)
              - fac::div(shearStrainStiffness_*theta)
              - p_
            );

            // Add stabilisation term
            // laplacian(w) is mathematically the same as div(grad(w)) but
            // numerically different for a given mesh (but they converge to the
            // TODO: we should store gradW to avoid repeatedly calculating it!
            if (alphaW > 0.0)
            {
                wEqn +=
                    alphaW
                   *(
                        fac::div(shearStrainStiffness_*gradW_)
                      - fac::laplacian(shearStrainStiffness_, w_)
                    );
            }

            // // Relax the linear system
            wEqn.relax();

            // Solve the linear system
            solverPerfw = wEqn.solve();

            // Relax the field
            w_.relax();

            // Update the gradient of displacement
            gradW_ = fac::grad(w_);

            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/
            // 2 - thetaX eqn
            // thetaX_.storePrevIter();

            // Initialise thetaX equation with implicit laplacian terms
            faScalarMatrix thetaXEqn
            (
                fam::laplacian(bendingStiffness_, thetaX_)
            );
            
            // Initialise thetaY equation with implicit laplacian terms
            faScalarMatrix thetaYEqn
            (
                fam::laplacian(bendingStiffness_, thetaY_)
            );

            if (alphaTheta > 0.0)
            {
                thetaXEqn +=
                    alphaTheta
                   *(
                        fac::div(bendingStiffness_*gradThetaX_)
                      - fac::laplacian(bendingStiffness_, thetaX_)
                    );

                thetaYEqn +=
                    alphaTheta
                   *(
                        fac::div(bendingStiffness_*gradThetaY_)
                      - fac::laplacian(bendingStiffness_, thetaY_)
                    );
            }

            // NOTE!! - How to get the moment arm (x - x_P) and (y - y_P)?
            // For an orthogonal uniform mesh, the above terms are 
            // half of cell to cell distances and are the equal when we
            // look from the owner and neighbour side. The only difference
            // is that the sign of (x - x_P) is positive for owner and
            // negative for neighbour. CHECK WITH IVAN?

            // Initialise the inverse of delta (cell to cell distance)
            const edgeScalarField invDeltaCoeffs
            (
                IOobject
                (
                    "invDeltaCoeffs",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                1.0/aMesh_.deltaCoeffs()
            );

            const scalarField invDeltaI(invDeltaCoeffs.internalField());

            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/
            // \int_{dSz} (x - x_p) p ds - moment arm due to pressure term


            // Firstly, it does not make sense how to multiply (x_e - x_P) 
            // that is defined at an edge with cell centre pressure field
            // CHECK WITH IVAN..
            
            // Also, according to Torlak, this term is zero if 
            // coordinate system is assumed at cell centre.
            // forAll(thetaXEqn.diag(), cI)
            // {
            //     thetaXEqn.source()[cI] += invDeltaI[cI]*p_[cI]*mag(Sf[cI]);
            //     thetaYEqn.source()[cI] += invDeltaI[cI]*p_[cI]*mag(Sf[cI]);
            // }


            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/

            // SHEAR FORCE MOMENT ARM TERM 
            // The last s_{\phi \l} term of thetaX and thetaY variables
            // in Demirdzic 1997 plate paper.
            // \int_{dl} Gamma (x - x_P) (grad(w) - theta) \cdot n dl
            // \Sum_e Gamma (x_e - x_P) (grad(w)_e - theta_e) \cdot n_e l_e

            // How to get (x_e - x_P) for every edge 'e'? This value is
            // different if looked from the owner side or neighbour side
            // and also contains the SIGN!!!

            // BUT: for orthogonal uniform mesh, it is same. Use that for now!!
            // | (x_e - x_P) | = 0.5*inv(mesh.deltaCoeffs()) for internal edges

            // SIGN: (x_e - x_P) \cdot n_e will always be positive as they
            // cancel each other's signs. CHECK WITH IVAN AGAIN.

            edgeVectorField gradWEdge(fac::interpolate(gradW_));

            // Avoid oscillations in gradient calculations
            if (compactEdgeNormalGrad)
            {
                const edgeScalarField lnGradWEdge(fac::lnGrad(w_));
                gradWEdge +=
                    lnGradWEdge*edgeBiNormal - (sqr(edgeBiNormal) & gradWEdge);
            }

            const vectorField gradWEdgeI(gradWEdge.internalField());

            const scalarField nx(edgeBiNormal.internalField().component(vector::X));
            const scalarField ny(edgeBiNormal.internalField().component(vector::Y));

            const scalarField leI(le.internalField());

            // Shear strain constant after removing dimension (Gamma = G*h)
            const scalar Gamma(shearStrainStiffness_.value());

            /*---------------------------------------------------------------*/
            // APPROACH - 1 (LOOP OVER EDGES)
            // Loop over internal edges
            const Switch loopEdgesShearForce
            (
                solidModelDict().lookup("loopEdgesShearForce")
            );

            if (loopEdgesShearForce)
            {
                forAll(thetaXEqn.upper(), eI)
                {
                    
                    // NOTE (check this again):
                    // Sign for owner and neighbour contribution is same because
                    // (x_e - x_P) \cdot n_e is always positive. The negative sign
                    // due to the normal is cancelled by the sign of (x_e - x_P)
                    // CHECK WITH IVAN AGAIN.

                    // For internal edges uniform orthogonal mesh,
                    // (x_e - x_P) = 0.5*invDeltaCoeffs

                    // Terms are added to RHS source. So, need to subtract from source
                    // Also, the neighbour contribution sign is same as owner.
                    // Hence, it is all minus sign
                    // "source -="

                    // thetaX part
                    thetaXEqn.source()[own[eI]] -=
                        0.5*Gamma*invDeltaI[eI]*leI[eI]
                       *(
                            (
                                gradWEdgeI[eI] - thetaEdge[eI]
                            ) & edgeBiNormal.internalField()[eI]
                        );

                    thetaXEqn.source()[nei[eI]] -=
                        0.5*Gamma*invDeltaI[eI]*leI[eI]
                       *(
                            (
                                gradWEdgeI[eI] - thetaEdge[eI]
                            ) & edgeBiNormal.internalField()[eI]
                        );

                    // thetaY part
                    thetaYEqn.source()[own[eI]] -=
                        0.5*Gamma*invDeltaI[eI]*leI[eI]
                       *(
                            (
                                gradWEdgeI[eI] - thetaEdge[eI]
                            ) & edgeBiNormal.internalField()[eI]
                        );

                    thetaYEqn.source()[nei[eI]] -=
                        0.5*Gamma*invDeltaI[eI]*leI[eI]
                       *(
                            (
                                gradWEdgeI[eI] - thetaEdge[eI]
                            ) & edgeBiNormal.internalField()[eI]
                        );
                }

                // Loop over boundary edges
                forAll(thetaX_.boundaryField(), patchI)
                {
                    const faePatchVectorField pGradW(gradWEdge.boundaryField()[patchI]);
                    const faePatchVectorField pThetaEdge(thetaEdge.boundaryField()[patchI]);

                    forAll(thetaX_.boundaryField()[patchI], pEdge)
                    {
                        // Boundary cell index
                        const label bI = faBouMesh[patchI].edgeFaces()[pEdge];
                        const scalar leB = le.boundaryField()[patchI][pEdge];
                        const scalar pInvDelta = invDeltaCoeffs.boundaryField()[patchI][pEdge];

                        // For boundary edges uniform orthogonal mesh,
                        // (x_e - x_P) = invDeltaCoeffs, no 0.5 coeff needed!!
                        thetaXEqn.source()[bI] -=
                            Gamma*pInvDelta*leB
                           *(
                                (pGradW[pEdge] - pThetaEdge[pEdge])
                              & edgeBiNormal.boundaryField()[patchI][pEdge]
                            );

                        thetaYEqn.source()[bI] -=
                            Gamma*pInvDelta*leB
                           *(
                                (pGradW[pEdge] - pThetaEdge[pEdge])
                              & edgeBiNormal.boundaryField()[patchI][pEdge]
                            );
                    }
                }
            }
            else
            {
                /*---------------------------------------------------------------*/

                // APPROACH 2 (LOOP OVER CELLS and for each cell LOOP over EDGES)
                forAll(theta.internalField(), cellI)
                {
                    forAll(aMesh_.internalEdges(), edgeI)
                    {
                        if (own[edgeI] == cellI)
                        {
                            thetaXEqn.source()[cellI] -=
                                Gamma*leI[edgeI]
                               *(
                                    gradWEdge.internalField()[edgeI].component(vector::X)
                                  - thetaEdge.internalField()[edgeI].component(vector::X)
                                )
                               *(
                                    edgeCentres[edgeI].component(vector::X) 
                                  - cellCentres[cellI].component(vector::X) 
                                )*nx[edgeI];

                            thetaYEqn.source()[cellI] -=
                                Gamma*leI[edgeI]
                                *(
                                    gradWEdge.internalField()[edgeI].component(vector::Y)
                                  - thetaEdge.internalField()[edgeI].component(vector::Y)
                                )
                               *(
                                    edgeCentres[edgeI].component(vector::Y) 
                                  - cellCentres[cellI].component(vector::Y) 
                                )*ny[edgeI];
                        }
                        else if (nei[edgeI] == cellI)
                        {
                            thetaXEqn.source()[cellI] -=
                                Gamma*leI[edgeI]
                               *(
                                    gradWEdge.internalField()[edgeI].component(vector::X)
                                  - thetaEdge.internalField()[edgeI].component(vector::X)
                                )
                               *(
                                    edgeCentres[edgeI].component(vector::X) 
                                  - cellCentres[cellI].component(vector::X) 
                                )*nx[edgeI];

                            thetaYEqn.source()[cellI] -=
                                Gamma*leI[edgeI]
                               *(
                                    gradWEdge.internalField()[edgeI].component(vector::Y)
                                  - thetaEdge.internalField()[edgeI].component(vector::Y)
                                )
                               *(
                                    edgeCentres[edgeI].component(vector::Y) 
                                  - cellCentres[cellI].component(vector::Y) 
                                )*ny[edgeI];
                        }
                    }
                }

                forAll(theta.boundaryField(), patchI)
                {
                    const faePatchVectorField pGradW(gradWEdge.boundaryField()[patchI]);
                    const faePatchVectorField pThetaEdge(thetaEdge.boundaryField()[patchI]);
                    const scalarField nxb
                    (
                        edgeBiNormal.boundaryField()[patchI].component(vector::X)
                    );

                    const scalarField nyb
                    (
                        edgeBiNormal.boundaryField()[patchI].component(vector::Y)
                    );
                    const vectorField pEdgeCentres(edgeCentres.boundaryField()[patchI]);
                    const vectorField pCellCentres(cellCentres.boundaryField()[patchI]);

                    forAll(theta.boundaryField()[patchI], pEdge)
                    {
                        // Boundary cell index
                        const label bI = faBouMesh[patchI].edgeFaces()[pEdge];
                        const scalar leB = le.boundaryField()[patchI][pEdge];

                        thetaXEqn.source()[bI] -=
                            Gamma*leB
                            *(
                                pGradW[pEdge].component(vector::X)
                              - pThetaEdge[pEdge].component(vector::X)
                            )
                           *(
                                pEdgeCentres[pEdge].component(vector::X) 
                              - pCellCentres[bI].component(vector::X) 
                            )*nxb[pEdge];

                        thetaYEqn.source()[bI] -=
                            Gamma*leB
                            *(
                                pGradW[pEdge].component(vector::Y)
                              - pThetaEdge[pEdge].component(vector::Y)
                            )
                           *(
                                pEdgeCentres[pEdge].component(vector::Y) 
                              - pCellCentres[bI].component(vector::Y) 
                            )*nyb[pEdge];
                    }
                }
            }

            Info<< "grad w edge " << gradWEdge << endl;
            Info<< "theta edge " << thetaEdge << endl;
            // Info<< "le " << le << endl;
            // Info<< "invDelta " << invDeltaCoeffs << endl;
            Info<< "thX source " << thetaXEqn.source() << endl;
            Info<< "thY source " << thetaYEqn.source() << endl;
            /*---------------------------------------------------------------*/



            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/
            // GRADIENT THETA TERMS OF THE SOURCE
            // The first three terms of s_{\phi \l} term of thetaX and thetaY variables
            // in Demirdzic 1997 plate paper.
            
            // Interpolate grad(thetaX), grad(thetaY) to edges
            edgeVectorField gradthetaXEdge(fac::interpolate(gradThetaX_));
            edgeVectorField gradthetaYEdge(fac::interpolate(gradThetaY_));

            // Info<< "gradthXE " << gradthetaXEdge << endl; 

            // Avoid oscillations in gradient calculations
            if (compactEdgeNormalGrad)
            {
                const edgeScalarField lnGradThetaXEdge(fac::lnGrad(thetaX_));
                const edgeScalarField lnGradThetaYEdge(fac::lnGrad(thetaY_));

                gradthetaXEdge +=
                    lnGradThetaXEdge*edgeBiNormal - (sqr(edgeBiNormal) & gradthetaXEdge);

                gradthetaYEdge +=
                    lnGradThetaYEdge*edgeBiNormal - (sqr(edgeBiNormal) & gradthetaYEdge);
            }

            // Extracting individual components of gradTheta at edges
            const edgeScalarField gradThXX
            (
                IOobject
                (
                    "gradThXX",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                gradthetaXEdge.component(vector::X)
            );
            const edgeScalarField gradThXY
            (
                IOobject
                (
                    "gradThXY",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                gradthetaXEdge.component(vector::Y)
            );
            const edgeScalarField gradThYX
            (
                IOobject
                (
                    "gradThYX",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                gradthetaYEdge.component(vector::X)
            );
            const edgeScalarField gradThYY
            (
                IOobject
                (
                    "gradThYY",
                    runTime().timeName(),
                    mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                gradthetaYEdge.component(vector::Y)
            );

            // Constants after removing dimensions
            const scalar D(bendingStiffness_.value());
            const scalar nu(nu_.value());

            // Contribution to internal edges
            // The terms when put in source (RHS) will change
            // signs. So the owner contributions are negative
            // and neighbour is positive. // CHECK WITH IVAN??
            forAll(thetaXEqn.upper(), eI)
            {
                // thetaX part
                thetaXEqn.source()[own[eI]] -=
                    D*leI[eI]
                   *(
                        nu*gradThYY[eI]*nx[eI]
                      + 0.5*(1 - nu)*gradThYX[eI]*ny[eI]
                      - 0.5*(1 + nu)*gradThXY[eI]*ny[eI]
                    );
                
                // thetaX part
                thetaXEqn.source()[nei[eI]] +=
                    D*leI[eI]
                   *(
                        nu*gradThYY[eI]*nx[eI]
                      + 0.5*(1 - nu)*gradThYX[eI]*ny[eI]
                      - 0.5*(1 + nu)*gradThXY[eI]*ny[eI]
                    );

                // thetaY part
                thetaYEqn.source()[own[eI]] -=
                    D*leI[eI]
                   *(
                        0.5*(1 - nu)*gradThXY[eI]*nx[eI]
                      + nu*gradThXX[eI]*ny[eI]
                      - 0.5*(1 + nu)*gradThYX[eI]*nx[eI]
                    );
                
                // thetaY part
                thetaYEqn.source()[nei[eI]] +=
                    D*leI[eI]
                   *(
                        0.5*(1 - nu)*gradThXY[eI]*nx[eI]
                      + nu*gradThXX[eI]*ny[eI]
                      - 0.5*(1 + nu)*gradThYX[eI]*nx[eI]
                    );
            }

            // Boundary edges
            forAll(thetaX_.boundaryField(), patchI)
            {
                const faePatchScalarField pGradThXX
                (
                    gradThXX.boundaryField()[patchI]
                );

                const faePatchScalarField pGradThXY
                (
                    gradThXY.boundaryField()[patchI]
                );

                const faePatchScalarField pGradThYX
                (
                    gradThYX.boundaryField()[patchI]
                );

                const faePatchScalarField pGradThYY
                (
                    gradThYY.boundaryField()[patchI]
                );

                const scalarField nxb
                (
                    edgeBiNormal.boundaryField()[patchI].component(vector::X)
                );

                const scalarField nyb
                (
                    edgeBiNormal.boundaryField()[patchI].component(vector::Y)
                );

                forAll(thetaX_.boundaryField()[patchI], pEdge)
                {
                    // Boundary cell index
                    const label bI = faBouMesh[patchI].edgeFaces()[pEdge];
                    const scalar leB = le.boundaryField()[patchI][pEdge];

                    // thetaX part
                    thetaXEqn.source()[bI] -=
                        D*leB
                       *(
                            nu*pGradThYY[pEdge]*nxb[pEdge]
                          + 0.5*(1 - nu)*pGradThYX[pEdge]*nyb[pEdge]
                          - 0.5*(1 + nu)*pGradThXY[pEdge]*nyb[pEdge]
                        );

                    // thetaY part
                    thetaYEqn.source()[bI] -=
                        D*leB
                       *(
                            0.5*(1 - nu)*pGradThXY[pEdge]*nxb[pEdge]
                          + nu*pGradThXX[pEdge]*nyb[pEdge]
                          - 0.5*(1 + nu)*pGradThYX[pEdge]*nxb[pEdge]
                        );
                }
            }
            
            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/
            // Solve the linear system
            solverPerfThetaX = thetaXEqn.solve();

            // Relax thetaX field
            thetaX_.relax();

            // Gradient of thetaX (TO BE CHANGED)
            gradThetaX_ = fac::grad(thetaX_);

            /*---------------------------------------------------------------*/
            /*---------------------------------------------------------------*/
            // 3 - thetaY eqn
            // Solve the linear system
            solverPerfThetaY = thetaYEqn.solve();

            thetaY_.relax();

            // Gradient of thetaY (TO BE CHANGED)
            gradThetaY_ = fac::grad(thetaY_);
        }
        while
        (
            !converged
            (
                iCorr,
                solverPerfw,
                solverPerfThetaX,
                solverPerfThetaY,
                w_,
                thetaX_,
                thetaY_
            )
            &&
            ++iCorr < nCorr()
        );

        // Map area fields to vol fields
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


tmp<vectorField> mindlinDemirdzicPlateSolid::tractionBoundarySnGrad
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


void mindlinDemirdzicPlateSolid::setTraction
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


void mindlinDemirdzicPlateSolid::writeFields(const Time& runTime)
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
