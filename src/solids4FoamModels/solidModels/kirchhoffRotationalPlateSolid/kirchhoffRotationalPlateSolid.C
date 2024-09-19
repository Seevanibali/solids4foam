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
            IOobject::NO_WRITE
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

        // Approach: Block - coupled formulation to Solve M, thetaX, thetaY and w equations
        // equations simulataneously
        // Iterative approach because the source vector values come from the previous iteration 
        // if (coupled_)
        // {
            
        //     Info<< "\nUsing block-coupled approach to solve for w, thetaX, thetaY, and M eqns" << endl;

        //     // Update boundary conditions
        //     M_.correctBoundaryConditions();
        //     w_.correctBoundaryConditions();
        //     thetaX_.correctBoundaryConditions();
        //     thetaY_.correctBoundaryConditions();

        //     const label nCells(aMesh_.faceCells().size());
        //     const labelList& own(aMesh_.owner());
        //     const labelList& nei(aMesh_.neighbour());
        //     const DimensionedField<scalar, areaMesh>& Sf(aMesh_.S());
        //     const edgeScalarField& le(aMesh_.magLe());
        //     const faBoundaryMesh& faBouMesh(aMesh_.boundary());

        //     // Initialise block matrix (2 scalar equations of w and M per cell)
        //     SparseMatrixTemplate<scalar> matrix(2*nCells);
            
        //     matrix.clear();

        //     // Initialise source vector
        //     scalarField source(2*nCells, 0.0);
            
        //     // Initialise solution field
        //     scalarField solveMw(2*nCells, 0.0);

        //     // d2dt2 term
        //     // Note: when running a case, it says available d2dt2 schemes are only Euler,
        //     // When the test case is steadyState, read ddtScheme to be steadyState and inertial
        //     // terms of d2dt2(w) are not added to the matrix.
        //     const faScalarMatrix d2dt2W(rho_*h_*fam::d2dt2(w_));
        //     const scalarField& d2dt2WDiag = d2dt2W.diag();
        //     const scalarField& d2dt2WSource = d2dt2W.source();
            
        //     // Calculate Laplacian discretisation of M (moment sum)
        //     const faScalarMatrix laplacianM(-fam::laplacian(M_));
        //     const scalarField& lapMDiag = laplacianM.diag();
        //     const scalarField& lapMUpper = laplacianM.upper();
        //     const FieldField<Field, scalar>& lapMIntCoeffs = laplacianM.internalCoeffs();
        //     const FieldField<Field, scalar>& lapMBouCoeffs = laplacianM.boundaryCoeffs();

        //     // Calculate Laplacian discretisation of w 
        //     const faScalarMatrix laplacianW(fam::laplacian(bendingStiffness_, w_));
        //     const scalarField& lapWDiag = laplacianW.diag();
        //     const scalarField& lapWUpper = laplacianW.upper();
        //     const FieldField<Field, scalar>& lapWIntCoeffs = laplacianW.internalCoeffs();
        //     const FieldField<Field, scalar>& lapWBouCoeffs = laplacianW.boundaryCoeffs();

        //     // Assembling the diagonal coeffs of MEqn and wEqn into a block matrix
        //     forAll(lapMDiag, i)
        //     {
        //         // Diagonals of the block matrix diagonal 
        //         // Coefficient of M in MEqn   
        //         matrix(2*i, 2*i) = lapMDiag[i];

        //         // Coefficient of w in wEqn
        //         matrix(2*i + 1, 2*i + 1) = lapWDiag[i];

        //         // Off-diagonals of the block diagonal
        //         // Coefficient of M in the wEqn
        //         matrix(2*i + 1, 2*i) = mag(Sf[i]);

        //         // Explicit coeffients of the MEqn go to RHS - source
        //         source[2*i] = p_[i]*mag(Sf[i]);

        //         // Add d2dt2 coeffs
        //         if(ddtSchemeName == "steadyState")
        //         {
        //             // Do not add any inertial contribution for d2dt2 terms
        //         }
        //         else if(ddtSchemeName != "steadyState" && d2dt2SchemeName == "Euler")
        //         {
        //             matrix(2*i, 2*i + 1) = d2dt2WDiag[i];

        //             source[2*i] += d2dt2WSource[i];
        //         }
        //         else
        //         {
        //             FatalError("evolve() function in kirchhoff plate rotattion-free solid") << nl
        //                 << "Incompatible (or not defined) d2dt2Scheme " 
        //                 << d2dt2SchemeName << " is specified! "
        //                 << abort(FatalError);
        //         }

        //     }

        //     // Off-diagonal components of the block matrix
        //     forAll(lapMUpper, faceI)
        //     {
        //         label i = own[faceI];
        //         label j = nei[faceI];
                
        //         // Coefficients of the upper part of the block matrix
        //         // Note: There is no neighbour contribution of w in MEqn and M in wEqn 
        //         // for rotation-free Kirchhoff plate equations
        //         // Neighbour contribution of M in MEqn
        //         matrix(2*i, 2*j) = lapMUpper[i];

        //         // Neighbour contribution of w in wEqn
        //         matrix(2*i + 1, 2*j + 1) = lapWUpper[i];

        //         // Coefficients of the lower part of the block matrix
        //         matrix(2*j, 2*i) = lapMUpper[i];
        //         matrix(2*j + 1, 2*i + 1) = lapWUpper[i];
        //     }

        //     // Loop over boundary patches
        //     forAll(M_.boundaryField(), patchI)
        //     {
        //         const word& patchTypeM(M_.boundaryField()[patchI].type());
        //         const List<scalar>& delta(aMesh_.boundary()[patchI].deltaCoeffs());

        //         // Loop over all faces of boundary patch
        //         forAll(M_.boundaryField()[patchI], faceI)
        //         {
        //             // Boundary cell index
        //             const label bI = faBouMesh[patchI].edgeFaces()[faceI];
        //             const scalar leB = le.boundaryField()[patchI][faceI];
        //             const scalar delB = delta[faceI];

        //             if(patchTypeM == "clampedMoment")
        //             {
        //                 // Coefficients of w in MEqn due to boundary coupling
        //                 // This implict coefficient is only considering an orthogonal mesh
        //                 // For non-orthogonal mesh, correction terms has to be added explicitly
        //                 // So an iterative loop is needed.
        //                 const scalar coeffBouW(-2*leB*bendingStiffness_.value()*pow(delB,3));
        //                 const scalar explicitW(w_.boundaryField()[patchI][faceI]);
                        
        //                 // Diagonal contribution for BC to M_ of the matrix in MEqn
        //                 matrix(2*bI, 2*bI) += lapMIntCoeffs[patchI][faceI];

        //                 // Diagonal contribution for BC to w_ of the matrix in MEqn
        //                 matrix(2*bI, 2*bI + 1) += coeffBouW;

        //                 // Explicit contribution of boundary edges to the source (MEqn)
        //                 source[2*bI] -= coeffBouW*explicitW;                 
        //             }
        //             else
        //             {
        //                 // Contribution of boundary edges to the diagonal of the matrix (MEqn)
        //                 matrix(2*bI, 2*bI) += lapMIntCoeffs[patchI][faceI];

        //                 // Explicit contribution of boundary edges to the source (MEqn)
        //                 source[2*bI] += lapMBouCoeffs[patchI][faceI];

        //                 // Contribution of boundary edges to the diagonal of the matrix (wEqn)
        //                 matrix(2*bI + 1, 2*bI + 1) += lapWIntCoeffs[patchI][faceI];

        //                 // Explicit contribution of boundary edges to the source (wEqn)
        //                 source[2*bI + 1] += lapWBouCoeffs[patchI][faceI];
        //             }   
        //         }
        //     }

        //     if(debug > 1)
        //     {
        //         Info<< "\nBlock Matrix Coefficients: " << matrix.data() << endl;
        //         Info<< "\nSource vector: " << source << endl;
        //     }

        //     // Solve the linear system of equations
        //     // Using Eigen SparseLU direct solver
        //     sparseMatrixTools::solveLinearSystemEigen
        //     (
        //         matrix, source, solveMw, false, debug
        //     );

        //     // Retrieve solution
        //     for(label i = 0; i < nCells; ++i)
        //     {
        //         M_[i] = solveMw[2*i];
        //         w_[i] = solveMw[2*i + 1];
        //     }

        //     // Correct the boundary conditions for M and w
        //     M_.correctBoundaryConditions();
        //     w_.correctBoundaryConditions();
        //     thetaX_.correctBoundaryConditions();
        //     thetaY_.correctBoundaryConditions();

        //     // Update the angle of rotation
        //     vectorField theta(-fac::grad(w_));
        //     // thetaX_(theta.component(0));
        //     // thetaY_(theta.component(1));

        //     // Update the gradient of rotation field, used for non-orthogonal
        //     // correction in clamped boundary conditions
        //     // gradTheta_ = fac::grad(theta_);
        // }
        // else
        // {
            // Approach 2: Using segregated method of solving M and w equations separately 
            // and then iteratively update them until the values fall below a solution tolerance

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

                // Access the cell centre position vectors of the mesh
                const areaVectorField& positionR(aMesh_.areaCentres());

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
                    - (positionR.component(0)*gradM_.component(0))
                );
                QThetaXY_ = 
                (
                    torsionalStiffness_*gradThYX_
                    - (positionR.component(0)*gradM_.component(1))
            
                );
                QThetaYX_ = 
                ( 
                    torsionalStiffness_*gradThXY_
                    - (positionR.component(1)*gradM_.component(0))
                );
                QThetaYY_ =
                (
                    ((bendingStiffness_ - torsionalStiffness_)*gradThYY_)
                    + (bendingStiffness_*nu_*gradThXX_)
                    - (positionR.component(1)*gradM_.component(1))
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

                // Info<< "*-----------------------------------*" << nl
                // // << "Rcell " << nl
                // // << positionRcell.unzip(positionRcell, List<scalar>& rX)<< nl
                // // << "Rbou" << nl <<  positionRbou
                // << nl
                // << "gradThXY " << gradThXY_ << nl
                // << "gradM " << gradM_ << nl 
                // << "QThetaX " << QThetaX_ << nl
                // << "QThetaY " << QThetaY_ << nl
                // // << "gradThXXI " << gradThXXI
                // << "*---------------------------------------*" << nl
                // << endl;

                // Solve thetaX equation 
                faScalarMatrix thetaXEqn
                (
                    fam::laplacian(torsionalStiffness_, thetaX_) + fac::div(QThetaX_)
                );

                // Solve the linear system
                thetaXEqn.solve();
                thetaX_.relax();

                // gradThetaX_ = fac::grad(thetaX_);

                // Solve thetaY equation
                faScalarMatrix thetaYEqn
                (
                    fam::laplacian(torsionalStiffness_, thetaY_) + fac::div(QThetaY_)
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

                const label nCells(aMesh_.faceCells().size());
                const labelList& edgeOwn(aMesh_.edgeOwner());
                const labelList& edgeNei(aMesh_.edgeNeighbour());
                const areaVectorField& posCellCentres(aMesh_.areaCentres());

                // Info<< "edge owners " << nl << edgeOwn << endl;

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
                                    deltaR = posCellCentres[edgeNei[edgeI]] 
                                        - posCellCentres[edgeOwn[edgeI]];
                                }
                                else
                                {
                                    deltaR = posCellCentres[edgeOwn[edgeI]] 
                                        - posCellCentres[edgeNei[edgeI]];
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
