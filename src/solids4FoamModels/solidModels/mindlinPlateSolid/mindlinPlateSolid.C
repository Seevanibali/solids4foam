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

#include "mindlinPlateSolid.H"
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

defineTypeNameAndDebug(mindlinPlateSolid, 0);
addToRunTimeSelectionTable(solidModel, mindlinPlateSolid, dictionary);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

bool mindlinPlateSolid::converged
(
    const int iCorr,
#ifdef OPENFOAM_NOT_EXTEND
    const SolverPerformance<scalar>& solverPerfw,
    const SolverPerformance<vector>& solverPerfTheta,
#else
    const lduSolverPerformance& solverPerfw,
    const lduSolverPerformance& solverPerfTheta,
#endif
    const areaScalarField& w,
    const areaVectorField& theta
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    // Calculate relative residuals
    // const scalar residualM =
    //     gMax
    //     (
    //         (
    //             mag(M - M.prevIter())
    //            /max
    //             (
    //                 gMax(mag(M - M.oldTime())()), SMALL
    //             )
    //         )()
    //     );

    const scalar resTheta =
        gMax
        (
            (
                mag(theta - theta.prevIter())
               /max
                (
                    gMax(mag(theta - theta.oldTime())()), SMALL
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
        bool convergedTheta = false;

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

        if
        (
            (
                mag(solverPerfTheta.initialResidual()) < solutionTol()
             && resTheta < solutionTol()
            )
         || mag(solverPerfTheta.initialResidual()) < alternativeTol()
         || resTheta < alternativeTol()
        )
        {
            convergedTheta = true;
        }



        if (convergedw && convergedTheta)
        {
            Info<< "    The residuals have converged" << endl;
            converged = true;
        }
    }

    // const bool var = true;
    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, solnRes (w & theta), relRes (w  & theta), iters (w & theta)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << solverPerfw.initialResidual()
            << ", " << solverPerfTheta.initialResidual()
            << ", " << tab << residualw
            << ", " << resTheta
            << ", " << tab << solverPerfw.nIterations()
            << ", " << solverPerfTheta.nIterations()
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


const fvPatch& mindlinPlateSolid::areaPatch() const
{
    if (areaPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaPatchID_];
}


const fvPatch& mindlinPlateSolid::areaShadowPatch() const
{
    if (areaShadowPatchID_ == -1)
    {
        calcAreaPatches();
    }

    return mesh().boundary()[areaShadowPatchID_];
}


void mindlinPlateSolid::calcAreaPatches() const
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
                            "void mindlinPlateSolid::calcAreaPatches() const"
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
            "void mindlinPlateSolid::calcAreaPatches() const"
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
            "void mindlinPlateSolid::calcAreaPatches() const"
        )   << "The polyMesh should be one cell thick, where there should be "
            << "two patches opposite each other that have the same number of "
            << "faces" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mindlinPlateSolid::mindlinPlateSolid
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
    // thetaX_
    // (
    //     IOobject
    //     (
    //         "thetaX",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     aMesh_
    // ),
    // thetaY_
    // (
    //     IOobject
    //     (
    //         "thetaY",
    //         runTime.timeName(),
    //         mesh(),
    //         IOobject::MUST_READ,
    //         IOobject::AUTO_WRITE
    //     ),
    //     aMesh_
    // ),
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
    theta_
    (
        IOobject
        (
            "theta",
            runTime.timeName(),
            mesh(),
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh_
    ),
    thetaVf_
    (
        IOobject
        (
            "thetaVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedVector("zero", dimless, vector::zero)
    ),
    gradTheta_
    (
        IOobject
        (
            "grad(" + theta_.name() + ")",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,   
            IOobject::AUTO_WRITE 
        ),
        aMesh_,
        dimensionedTensor("zero", dimLength/dimArea, tensor::zero)
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
        dimensionedVector("zero", dimForce, vector::zero)
    ),
    MsumVf_
    (
        IOobject
        (
            "MsumVf",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimForce, 0.0)
    ),
    rho_("zero", dimDensity, 0.0),
    E_("zero", dimPressure, 0.0),
    nu_("zero", dimless, 0.0),
    h_(solidModelDict().lookup("plateThickness")),
    shearCorrectionFactor_(solidModelDict().getOrDefault<scalar>("shearCorrectionFactor", 0.8333)),
    bendingStiffness_("zero", dimPressure*dimVolume, 0.0),
    torsionalStiffness_("zero", dimPressure*dimVolume, 0.0),
    shearStrainStiffness_("zero",dimForce/dimLength, 0.0),
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

    // Set the area vector field theta_ with x-component thetaX and y-component thetaY
    /*
    Info<< "Setting the theta vector field from two scalar inputs thetaX and thetaY \n" << endl;

    forAll(theta_.internalField(), cellI)
    {
        theta_[cellI] = vector(thetaX_[cellI],thetaY_[cellI], 0.0);
    }

    forAll(theta_.boundaryField(), patchI)
    {
        faPatchVectorField& pTheta = theta_.boundaryFieldRef()[patchI];
        const faPatchScalarField& pThetaX = thetaX_.boundaryField()[patchI];
        const faPatchScalarField& pThetaY = thetaY_.boundaryField()[patchI];

        forAll(pTheta, edgeI)
        {
            pTheta[edgeI] = vector(pThetaX[edgeI], pThetaY[edgeI], 0);
        }
    }

    Info<< "theta vector field is set as below \n" << theta_ << endl;
    */
    // Set plate properties
    rho_ = mech.rhoScalar();
    E_ = mech.E();
    nu_ = mech.nu();
    bendingStiffness_ = E_*pow(h_, 3)/(12*(1 - pow(nu_, 2)));
    torsionalStiffness_ = 0.5*(1 - nu_)*bendingStiffness_;
    shearStrainStiffness_ = shearCorrectionFactor_*0.5*E_*h_/(1 + nu_);

    Info<< "Plate mechanical properties\n"
        << "Bending Stiffness\n" << bendingStiffness_ << "\n"
        << "\nTorsional Stiffness\n" << torsionalStiffness_ << "\n"
        << "\nShear strain stiffness\n" << shearStrainStiffness_
        << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


bool mindlinPlateSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Create volume-to surface mapping object
    volSurfaceMapping vsm(aMesh_);

    // Mesh update loop
    do
    {
        int iCorr = 0;
#ifdef OPENFOAM_NOT_EXTEND
        SolverPerformance<scalar> solverPerfw;
        SolverPerformance<vector> solverPerfTheta;
        SolverPerformance<scalar>::debug = 0;
#else
        lduSolverPerformance solverPerfw;
        lduSolverPerformance solverPerfTheta;
        blockLduMatrix::debug = 0;
#endif

        Info<< "Solving the Mindlin (thick plates) equation for primary variables w, thetaX, thetaY" << endl;

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

        WarningIn("evolve() function in mindlin (thick) plate solid model")<< nl
            << "d2dt2Scheme in system/faSchemes cannot take steadyState as a valid keyword!" << nl
            << "If you want plate-case to be solved for steady state condition, "
            << "set ddtScheme to be steadyState instead!! " << endl;


        // Store the previous iteration values for computing source vector
        // thetaX_.storePrevIter();
        // thetaY_.storePrevIter();
        // gradW_.storePrevIter();

        Info<< "\nUsing segregated approach to solve for w, thetaX, and thetaY eqns "
            << "separately and iteratively update them!" << endl;

        // Philip testing
        const scalar alphaW(readScalar(solidModelDict().lookup("alphaW")));
        const scalar alphaTheta(readScalar(solidModelDict().lookup("alphaTheta")));

        do
        {
            // Store fields for under-relaxation and residual calculation
            w_.storePrevIter();

            // Solve w equation
            // Also, "==" complains so we will move all terms to left
            faScalarMatrix wEqn
            (
                fam::laplacian(shearStrainStiffness_, w_)
              + fac::div(shearStrainStiffness_*theta_)
              + p_
            );

            // Add stabilisation term
            // laplacian(w) is mathematically the same as div(grad(w)) but
            // numerically different for a given mesh (but they converge to the
            // same answer)
            // It may be possible that laplacian(w) is more likely to "lock" in
            // this formulation, while div(grad(w)) may not. As we do not have
            // fam::div(grad(w)), the implicit component must be fam::laplacian
            // but we can change the explicit components such that the laplacian
            // is replaced by div(grad(w)) at convergence
            // Below I add a blending factor alphaW to let us decide how much of
            // div(grad(w)) to use. alphaW = 1 means use all div(grad(w)), while
            // alphaW = 0 means use all laplacian(w). In between values use a
            // blend
            // TODO: we should store gradW to avoid repeatedly calculating it!
            if (alphaW > 0.0)
            {
                wEqn +=
                    alphaW
                   *(
                        fac::div(shearStrainStiffness_*fac::grad(w_))
                      - fac::laplacian(shearStrainStiffness_, w_)
                    );
            }

            // d2dt2 can only take Euler as keyword, but if the user wants it to
            // be steadyState, it cannot happen. Hence check for ddtScheme
            // and add inertial terms for not steady state!!
            if(ddtSchemeName != "steadyState")
            {
                wEqn -= rho_*h_*(fac::d2dt2(w_));
            }

            // Relax the linear system
            wEqn.relax();

            // Solve the linear system
            solverPerfw = wEqn.solve();

            // Relax the field
            w_.relax();

            // Info<< "w field\n" << w_ << endl;
            // Update the gradient of displacement
            // gradW_ = fac::grad(w_);


            // Store fields for under-relaxation and residual calculation
            theta_.storePrevIter();

            // Solve theta equation
            faVectorMatrix thetaEqn
            (
                // div(grad.T) has the same diagonal entry as div(grad) so we
                // could boost our Laplacian diagonal and then explicitly
                // subtract the difference. From some quick tests, it seems
                // that this does not have much effect (this will be true if the
                // bendingStiffness is much less than the torsionalStiffness
                // We have a few ways to write these two terms:
                // Option 1
                fam::laplacian(torsionalStiffness_, theta_)
              + 0.5*bendingStiffness_*(1 + nu_)*fac::grad(fac::div(theta_))
              //   // Option 2
              //   fam::laplacian(torsionalStiffness_ + 0.5*bendingStiffness_*(1 + nu_), theta_)
              // - 0.5*bendingStiffness_*(1 + nu_)*fac::laplacian(theta_)
              // + 0.5*bendingStiffness_*(1 + nu_)*fac::div(T(fac::grad(theta_)))
              //   // Option 3
              //   fam::laplacian(torsionalStiffness_ + 0.5*bendingStiffness_*(1 + nu_), theta_)
              // - 0.5*bendingStiffness_*(1 + nu_)*fac::div(fac::grad(theta_))
              // + 0.5*bendingStiffness_*(1 + nu_)*fac::div(T(fac::grad(theta_)))

                // Other terms
              - shearStrainStiffness_*(fac::grad(w_))
              - fam::Sp(shearStrainStiffness_, theta_)
            );

            // Add stabilisation term
            // Same ideas as alphaW. From initial testing, these terms seem to
            // have a negligible effect, so alphaTheta = 0 is probably fine
            if (alphaTheta > 0.0)
            {
                thetaEqn +=
                    alphaTheta
                   *(
                        fac::div(torsionalStiffness_*fac::grad(theta_))
                      - fac::laplacian(torsionalStiffness_, theta_)
                    );
            }

            // d2dt2 can only take Euler as keyword, but if the user wants it to
            // be steadyState, it cannot happen. Hence check for ddtScheme
            // and add inertial terms for not steady state!!
            if(ddtSchemeName != "steadyState")
            {
                thetaEqn -= rho_*pow(h_,3)*(fac::d2dt2(theta_))/12;
            }

            // Solve the linear system
            solverPerfTheta = thetaEqn.solve();

            // Relax theta fields
            theta_.relax();

            // Gradient of theta
            gradTheta_ = fac::grad(theta_);
        }
        while
        (
            !converged
            (
                iCorr,
                solverPerfw,
                solverPerfTheta,
                w_,
                theta_
            )
            &&
            ++iCorr < nCorr()
        );

        // Calculation of moment vector fields
        const tensorField& gradThetaI = gradTheta_.internalField();

        areaVectorField M
        (
            IOobject
            (
                "M",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh_,
            dimensionedVector("zero", dimForce, vector::zero)
        );

        areaScalarField Msum
        (
            IOobject
            (
                "Msum",
                runTime().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            aMesh_,
            dimensionedScalar("zero", dimForce, 0.0)
        );

        // Evaluate Mxx and Myy for all the internal fields
        forAll(M.internalField(), cellI)
        {
            M[cellI].component(vector::X) =
                bendingStiffness_.value()
                *(
                    gradThetaI[cellI].component(tensor::XX)
                    + nu_.value()*gradThetaI[cellI].component(tensor::YY)
                 );

            M[cellI].component(vector::Y) =
                bendingStiffness_.value()
                *(
                    gradThetaI[cellI].component(tensor::YY)
                    + nu_.value()*gradThetaI[cellI].component(tensor::XX)
                 );

            Msum[cellI] =
                (
                    M[cellI].component(vector::X) + M[cellI].component(vector::Y)
                )/(1 + nu_.value());
        }


        // Evaluate Mxx and Myy for all the boundary patches
        forAll(M.boundaryField(), patchI)
        {
            const faPatchTensorField& pGradTheta = gradTheta_.boundaryField()[patchI];
            faPatchVectorField& pM = M.boundaryFieldRef()[patchI];
            faPatchScalarField& pMsum = Msum.boundaryFieldRef()[patchI];

            forAll(pM, edgeI)
            {
                pM[edgeI].component(vector::X) =
                    bendingStiffness_.value()
                    *(
                        pGradTheta[edgeI].component(tensor::XX)
                        + nu_.value()*pGradTheta[edgeI].component(tensor::YY)
                    );

                pM[edgeI].component(vector::Y) =
                    bendingStiffness_.value()
                    *(
                        pGradTheta[edgeI].component(tensor::YY)
                        + nu_.value()*pGradTheta[edgeI].component(tensor::XX)
                    );

                pMsum[edgeI] =
                (
                    pM[edgeI].component(vector::X) + pM[edgeI].component(vector::Y)
                )/(1 + nu_.value());
            }
        }

        // Map area fields to vol fields
        mapAreaFieldToSingleLayerVolumeField(w_, wVf_);
        // mapAreaFieldToSingleLayerVolumeField(thetaX_, thetaXVf_);
        mapAreaFieldToSingleLayerVolumeField(M, MVf_);
        mapAreaFieldToSingleLayerVolumeField(Msum, MsumVf_);
        mapAreaFieldToSingleLayerVolumeField(theta_, thetaVf_);
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


tmp<vectorField> mindlinPlateSolid::tractionBoundarySnGrad
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


void mindlinPlateSolid::setTraction
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


void mindlinPlateSolid::writeFields(const Time& runTime)
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
