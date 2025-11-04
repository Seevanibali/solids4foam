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
                solverPerfw.initialResidual() < solutionTol()
             &&
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
                solverPerfThetaX.initialResidual() < solutionTol()
             &&
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
                solverPerfThetaY.initialResidual() < solutionTol()
             &&
                resThetaY < solutionTol()
            )
        //  || solverPerfThetaY.initialResidual() < alternativeTol()
        //  || resThetaY < alternativeTol()
        )
        {
            convergedThetaY = true;
        }

        if (convergedw && (convergedThetaX && convergedThetaY))
        // if (convergedThetaX && convergedThetaY)
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
        // Info<< "    Corr, solnRes (thetaX, thetaY), relRes (thetaX, thetaY), iters (thetaX, thetaY)"
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


bool mindlinDemirdzicPlateSolid::blockConverged
(
    const int iCorr,
    const areaScalarField& w,
    const areaScalarField& thetaX,
    const areaScalarField& thetaY
)
{
    // We will check a number of different residuals for convergence
    bool converged = false;

    const scalar resW =
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

    if
    (
        resW < solutionTol()
     && (
            resThetaX < solutionTol()
         && resThetaY < solutionTol()
        )
    )
    {
        converged = true;
    }

    // Print residual information
    if (iCorr == 0)
    {
        Info<< "    Corr, solnRes (w, thetaX, thetaY)"
            << endl;
    }
    else if (iCorr % infoFrequency() == 0 || converged)
    {
        Info<< "    " << iCorr
            << ", " << tab << resW
            << ", " << resThetaX
            << ", " << resThetaY
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


const vectorField mindlinDemirdzicPlateSolid::calculateShearForceContribution
(
    const edgeVectorField& shearForceEdge
) const
{
    // Initialise the source contribution vector
    vectorField sf2source(aMesh_.faceCells().size(), vector::zero);

    // Mesh information required
    edgeVectorField edgeBiNormal("edgeBiNormal", aMesh_.Le()/aMesh_.magLe());
    const labelList& own(aMesh_.owner());
    const labelList& nei(aMesh_.neighbour());
    const edgeScalarField& le(aMesh_.magLe());
    const faBoundaryMesh& faBouMesh(aMesh_.boundary());
    const edgeVectorField& edgeCentres(aMesh_.edgeCentres());
    const areaVectorField& cellCentres(aMesh_.areaCentres());
    const scalarField leI(le.internalField());
    // const labelList& edgeOwn(aMesh_.edgeOwner());

    const vectorField shearForceI(shearForceEdge.internalField());
    const vectorField edgeBiNormalI(edgeBiNormal.internalField());

    // forAll(thetaX_.internalField(), cellI)
    // {
    //     forAll(aMesh_.internalEdges(), edgeI)
    //     {
    //         if (own[edgeI] == cellI)
    //         {
    //             // thetaXEqn.source()[cellI] -=
    //             sf2source[cellI].x() -=
    //                 leI[edgeI]
    //                 *(
    //                     edgeCentres[edgeI].component(vector::X)
    //                   - cellCentres[cellI].component(vector::X)
    //                 )
    //                 *(shearForceI[edgeI] & edgeBiNormalI[edgeI]);

    //             // thetaYEqn.source()[cellI] -=
    //             sf2source[cellI].y() -=
    //                 leI[edgeI]
    //                 *(
    //                     edgeCentres[edgeI].component(vector::Y)
    //                   - cellCentres[cellI].component(vector::Y)
    //                 )
    //                 *(shearForceI[edgeI] & edgeBiNormalI[edgeI]);
    //         }
    //         else if (nei[edgeI] == cellI)
    //         {
    //             // Note: we use "+=" as nx and ny need to be flipped
    //             // thetaXEqn.source()[cellI] +=
    //             sf2source[cellI].x() +=
    //                 leI[edgeI]
    //                 *(
    //                     edgeCentres[edgeI].component(vector::X)
    //                   - cellCentres[cellI].component(vector::X)
    //                 )
    //                 *(shearForceI[edgeI] & edgeBiNormalI[edgeI]);

    //             // Note: we use "+=" as nx and ny need to be flipped
    //             // thetaYEqn.source()[cellI] +=
    //             sf2source[cellI].y() +=
    //                 leI[edgeI]
    //                 *(
    //                     edgeCentres[edgeI].component(vector::Y)
    //                   - cellCentres[cellI].component(vector::Y)
    //                 )
    //                 *(shearForceI[edgeI] & edgeBiNormalI[edgeI]);
    //         }
    //     }
    // }

    // Directly add source contributions to own and nei by looping over edges
    forAll(aMesh_.internalEdges(), edgeI)
    {
        const scalar sfCoeff =
            leI[edgeI]*(shearForceI[edgeI] & edgeBiNormalI[edgeI]);

        const vector dROwn = (edgeCentres[edgeI] - cellCentres[own[edgeI]]);
        const vector dRNei = (edgeCentres[edgeI] - cellCentres[nei[edgeI]]);

        sf2source[own[edgeI]].x() -= sfCoeff*dROwn.x();
        sf2source[own[edgeI]].y() -= sfCoeff*dROwn.y();

        // Note: we use "+=" as nx and ny need to be flipped
        sf2source[nei[edgeI]].x() += sfCoeff*dRNei.x();
        sf2source[nei[edgeI]].y() += sfCoeff*dRNei.y();
    }

    forAll(thetaX_.boundaryField(), patchI)
    {
        const vectorField pEdgeCentres(edgeCentres.boundaryField()[patchI]);
        // const vectorField pCellCentres(cellCentres.boundaryField()[patchI]);
        const faePatchVectorField pShearForce(shearForceEdge.boundaryField()[patchI]);
        const vectorField pEdgeBiNormal(edgeBiNormal.boundaryField()[patchI]);

        forAll(thetaX_.boundaryField()[patchI], pEdge)
        {
            // Boundary cell index
            const label bI(faBouMesh[patchI].edgeFaces()[pEdge]);
            const scalar leB(le.boundaryField()[patchI][pEdge]);

            const scalar pSfCoeff =
                leB*(pShearForce[pEdge] & pEdgeBiNormal[pEdge]);

            const vector pDr = (pEdgeCentres[pEdge] - cellCentres[bI]);

            sf2source[bI].x() -= pSfCoeff*pDr.x();
            sf2source[bI].y() -= pSfCoeff*pDr.y();

            // // thetaXEqn.source()[bI] -=
            // sf2source[bI].x() -=
            //     leB
            //     *(
            //         pEdgeCentres[pEdge].component(vector::X)
            //       - pCellCentres[pEdge].component(vector::X)
            //     )
            // *(pShearForce[pEdge] & pEdgeBiNormal[pEdge]);

            // // thetaYEqn.source()[bI] -=
            // sf2source[bI].y() -=
            //     leB
            //     *(
            //         pEdgeCentres[pEdge].component(vector::Y)
            //       - pCellCentres[pEdge].component(vector::Y)
            //     )
            // *(pShearForce[pEdge] & pEdgeBiNormal[pEdge]);
        }
    }
    return sf2source;
}

const vectorField mindlinDemirdzicPlateSolid::calculateGradientThetaContribution
(
    const edgeVectorField& gradthetaXEdge,
    const edgeVectorField& gradthetaYEdge
) const
{
    // Initialise the source contribution vector
    vectorField gradTh2Source(aMesh_.faceCells().size(), vector::zero);

    // Mesh information required
    edgeVectorField edgeBiNormal("edgeBiNormal", aMesh_.Le()/aMesh_.magLe());
    const labelList& own(aMesh_.owner());
    const labelList& nei(aMesh_.neighbour());
    const edgeScalarField& le(aMesh_.magLe());
    const faBoundaryMesh& faBouMesh(aMesh_.boundary());
    const scalarField leI(le.internalField());

    const scalarField nx(edgeBiNormal.internalField().component(vector::X));
    const scalarField ny(edgeBiNormal.internalField().component(vector::Y));

    // Extracting individual components of gradTheta at edges
    const edgeScalarField gradThXX
    (
        IOobject
        (
            "gradThXX",
            runTime().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
            IOobject::NO_WRITE
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
    forAll(aMesh_.internalEdges(), eI)
    {
        // thetaX part
        // thetaXEqn.source()[own[eI]] -=
        gradTh2Source[own[eI]].x() -=
            D*leI[eI]
        *(
            nu*gradThYY[eI]*nx[eI]
          + 0.5*(1 - nu)*gradThYX[eI]*ny[eI]
          - 0.5*(1 + nu)*gradThXY[eI]*ny[eI]
        );

        // thetaX part
        // thetaXEqn.source()[nei[eI]] +=
        gradTh2Source[nei[eI]].x() +=
            D*leI[eI]
        *(
            nu*gradThYY[eI]*nx[eI]
          + 0.5*(1 - nu)*gradThYX[eI]*ny[eI]
          - 0.5*(1 + nu)*gradThXY[eI]*ny[eI]
        );

        // thetaY part
        // thetaYEqn.source()[own[eI]] -=
        gradTh2Source[own[eI]].y() -=
            D*leI[eI]
        *(
            0.5*(1 - nu)*gradThXY[eI]*nx[eI]
          + nu*gradThXX[eI]*ny[eI]
          - 0.5*(1 + nu)*gradThYX[eI]*nx[eI]
        );

        // thetaY part
        // thetaYEqn.source()[nei[eI]] +=
        gradTh2Source[nei[eI]].y() +=
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
            // thetaXEqn.source()[bI] -=
            gradTh2Source[bI].x() -=
                D*leB
            *(
                nu*pGradThYY[pEdge]*nxb[pEdge]
              + 0.5*(1 - nu)*pGradThYX[pEdge]*nyb[pEdge]
              - 0.5*(1 + nu)*pGradThXY[pEdge]*nyb[pEdge]
            );

            // thetaY part
            // thetaYEqn.source()[bI] -=
            gradTh2Source[bI].y() -=
                D*leB
            *(
                0.5*(1 - nu)*pGradThXY[pEdge]*nxb[pEdge]
              + nu*pGradThXX[pEdge]*nyb[pEdge]
              - 0.5*(1 + nu)*pGradThYX[pEdge]*nxb[pEdge]
            );
        }
    }

    return gradTh2Source;
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
            IOobject::NO_READ,
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
    wAnalytical_
    (
        IOobject
        (
            "w.analytical",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        aMesh_
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
    coupled_(solidModelDict().getOrDefault<bool>("coupled", false)),
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

        Info<< "Solving the Mindlin (thick plates) equation for "
            << "primary variables w, thetaX, thetaY - Demirdzic's Approach"
            << endl;

        // Philip testing
        const scalar alphaW(readScalar(solidModelDict().lookup("alphaW")));
        const scalar alphaTheta(readScalar(solidModelDict().lookup("alphaTheta")));

        // Read the time scheme from the test case
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


        // Note: To get in-plane normal unit vectors to an edge, aMesh_.Le()
        // can be used with unit norm
        // Do not use aMesh.unitLe() member from faMesh, since the
        // boundary values of aMesh.unitLe() are set to
        // "calculated (0 0 0)" which is not correct!!
        // Constructing the unit edgeBiNormal instead
        // For 2-D meshes, aMesh_.edgeNormals() gives unit vector in z-direction
        edgeVectorField edgeBiNormal("edgeBiNormal", aMesh_.Le()/aMesh_.magLe());

        // Mesh information required
        const label nCells(aMesh_.faceCells().size());
        const labelList& own(aMesh_.owner());
        const labelList& nei(aMesh_.neighbour());
        const edgeScalarField& interpWeights(aMesh_.weights());
        const DimensionedField<scalar, areaMesh>& Sf(aMesh_.S());
        const edgeScalarField& le(aMesh_.magLe());
        const faBoundaryMesh& faBouMesh(aMesh_.boundary());
        const edgeVectorField& edgeCentres(aMesh_.edgeCentres());
        const areaVectorField& cellCentres(aMesh_.areaCentres());
        const edgeScalarField& deltaCoeffs(aMesh_.deltaCoeffs());
        // const labelList& edgeOwn(aMesh_.edgeOwner());
        // const labelList& edgeNei(aMesh_.edgeNeighbour());

        // Info<< " aMesh edgeCentres " << edgeCentres << endl;
        // Info<< " aMesh cellCentres " << cellCentres << endl;
        // Info<< " aMesh owners " << own << endl;
        // Info<< " aMesh Nei " << nei << endl;


        const scalarField nx(edgeBiNormal.internalField().component(vector::X));
        const scalarField ny(edgeBiNormal.internalField().component(vector::Y));

        // NOTE!! - This theta construction was outside the do-loop in
        // previous commit. So, the theta contribution was not getting
        // added to w equation. But now that I add theta contribution
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

        // Store the previous iteration values for computing source vector
        // The storePrevIter values are also brought inside this do-loop
        // Should not make much difference.

        // w_.storePrevIter();
        // thetaX_.storePrevIter();
        // thetaY_.storePrevIter();
        // gradW_.storePrevIter();
        // gradThetaX_.storePrevIter();
        // gradThetaY_.storePrevIter();

        // w_.storeOldTime();
        // thetaX_.storeOldTime();
        // thetaY_.storeOldTime();
        // gradW_.storeOldTime();
        // gradThetaX_.storeOldTime();
        // gradThetaY_.storeOldTime();

        // Approach 1: Block - coupled formulation, Solve w, thetaX, and thetaY
        // equations simulataneously
        if (coupled_)
        {

            Info<< "\nUsing implicit block-coupled approach to solve "
                << "Mindlin equations for w, thetaX, and thetaY eqns"
                << endl;
            do
            {

                w_.storePrevIter();
                thetaX_.storePrevIter();
                thetaY_.storePrevIter();
                gradW_.storePrevIter();
                gradThetaX_.storePrevIter();
                gradThetaY_.storePrevIter();

                // Update boundary conditions
                w_.correctBoundaryConditions();
                thetaX_.correctBoundaryConditions();
                thetaY_.correctBoundaryConditions();

                // Theta vector at edge centres
                const edgeVectorField thetaEdge(fac::interpolate(theta));

                // Initialise block matrix
                // (3 scalar equations of w, thetaX, and thetaY per cell)
                SparseMatrixTemplate<scalar> matrix(3*nCells);

                matrix.clear();

                // Initialise source vector
                scalarField source(3*nCells, 0.0);

                // Initialise solution field
                scalarField solveWTheta(3*nCells, 0.0);

                // d2dt2 term
                // Note: when running a case, it says available d2dt2 schemes are only Euler,
                // When the test case is steadyState, read ddtScheme to be steadyState and inertial
                // terms of d2dt2(w) are not added to the matrix.
                const faScalarMatrix d2dt2W(rho_*h_*fam::d2dt2(w_));
                const scalarField& d2dt2WDiag = d2dt2W.diag();
                const scalarField& d2dt2WSource = d2dt2W.source();

                const faScalarMatrix d2dt2ThetaX((1/12)*rho_*pow(h_,3)*fam::d2dt2(thetaX_));
                const scalarField& d2dt2ThetaXDiag = d2dt2ThetaX.diag();
                const scalarField& d2dt2ThetaXSource = d2dt2ThetaX.source();

                const faScalarMatrix d2dt2ThetaY((1/12)*rho_*pow(h_,3)*fam::d2dt2(thetaY_));
                const scalarField& d2dt2ThetaYDiag = d2dt2ThetaY.diag();
                const scalarField& d2dt2ThetaYSource = d2dt2ThetaY.source();


                /*-----------LAPLACIAN TERMS FOR BLOCK DIAGONALS--------------*/
                // Calculate Laplacian discretisation of w
                const faScalarMatrix laplacianW
                (
                    fam::laplacian(shearStrainStiffness_, w_)
                );

                const scalarField& lapWDiag = laplacianW.diag();
                const scalarField& lapWUpper = laplacianW.upper();

                const FieldField<Field, scalar>&
                    lapWIntCoeffs = laplacianW.internalCoeffs();
                const FieldField<Field, scalar>&
                    lapWBouCoeffs = laplacianW.boundaryCoeffs();

                // Calculate Laplacian discretisation of thetaX
                const faScalarMatrix laplacianThetaX
                (
                    fam::laplacian(bendingStiffness_, thetaX_)
                );

                const scalarField& lapThetaXDiag = laplacianThetaX.diag();
                const scalarField& lapThetaXUpper = laplacianThetaX.upper();

                const FieldField<Field, scalar>&
                    lapThetaXIntCoeffs = laplacianThetaX.internalCoeffs();
                const FieldField<Field, scalar>&
                    lapThetaXBouCoeffs = laplacianThetaX.boundaryCoeffs();

                // Calculate Laplacian discretisation of thetaY
                const faScalarMatrix laplacianThetaY
                (
                    fam::laplacian(bendingStiffness_, thetaY_)
                );

                const scalarField& lapThetaYDiag = laplacianThetaY.diag();
                const scalarField& lapThetaYUpper = laplacianThetaY.upper();

                const FieldField<Field, scalar>&
                    lapThetaYIntCoeffs = laplacianThetaY.internalCoeffs();
                const FieldField<Field, scalar>&
                    lapThetaYBouCoeffs = laplacianThetaY.boundaryCoeffs();

                // Explicit shear force calculation
                // const areaVectorField gradW(fac::grad(w_));
                // edgeVectorField gradWEdge(fac::interpolate(gradW));

                // Avoid oscillations in gradient calculations
                // if (compactEdgeNormalGrad)
                // {
                //     const edgeScalarField lnGradWEdge(fac::lnGrad(w_));
                //     gradWEdge +=
                //         lnGradWEdge*edgeBiNormal - (sqr(edgeBiNormal) & gradWEdge);
                // }

                // const edgeVectorField shearForceEdge
                // (
                //     shearStrainStiffness_*(gradWEdge - thetaEdge)
                // );

                // Removed the thetaEdge term because I have added it implicitly now
                // const edgeVectorField shearForceEdge
                // (
                //     shearStrainStiffness_*gradWEdge
                // );

                // This vectorField recieves the contribution to source
                // along with the changed sign when terms are moved to RHS
                // const vectorField shearForceContrib =
                //     calculateShearForceContribution(shearForceEdge);

                // Explicit gradient of theta calculation
                // Interpolate grad(thetaX), grad(thetaY) to edges
                edgeVectorField gradthetaXEdge(fac::interpolate(gradThetaX_));
                edgeVectorField gradthetaYEdge(fac::interpolate(gradThetaY_));

                // Avoid oscillations in gradient calculations
                if (compactEdgeNormalGrad)
                {
                    const edgeScalarField lnGradThetaXEdge(fac::lnGrad(thetaX_));
                    const edgeScalarField lnGradThetaYEdge(fac::lnGrad(thetaY_));

                    gradthetaXEdge +=
                        lnGradThetaXEdge*edgeBiNormal
                      - (sqr(edgeBiNormal) & gradthetaXEdge);

                    gradthetaYEdge +=
                        lnGradThetaYEdge*edgeBiNormal
                      - (sqr(edgeBiNormal) & gradthetaYEdge);
                }

                // Func call to evaluate grad theta contributions to source
                const vectorField gradThetaContrib =
                    calculateGradientThetaContribution
                    (
                        gradthetaXEdge,
                        gradthetaYEdge
                    );

                // Assembling the diagonal coeffs of wEqn, thetaXEqn,
                // and thetaYEqn into a block matrix
                forAll(lapWDiag, i)
                {
                    // Diagonals of the block matrix diagonal
                    // Coefficient of w in the wEqn
                    matrix(3*i, 3*i) = lapWDiag[i];

                    // Coefficient of thetaX in thetaXEqn
                    matrix(3*i + 1, 3*i + 1) = lapThetaXDiag[i];

                    // Coefficient of thetaY in thetaYEqn
                    matrix(3*i + 2, 3*i + 2) = lapThetaYDiag[i];

                    // Explicit coeffients of the wEqn go to RHS - source
                    source[3*i] = p_[i]*mag(Sf[i]);


                    // Source terms - fac::grad(w) and gradient of theta terms
                    // to be added here
                    // The signs are all positive here because they are already
                    // included before
                    // source[3*i + 1] += shearForceContrib[i].component(vector::X);
                    // source[3*i + 2] += shearForceContrib[i].component(vector::Y);
                    source[3*i + 1] += gradThetaContrib[i].component(vector::X);
                    source[3*i + 2] += gradThetaContrib[i].component(vector::Y);

                    // Add d2dt2 coeffs
                    if(ddtSchemeName == "steadyState")
                    {
                        // Do not add any inertial contribution for d2dt2 terms
                    }
                    else if(ddtSchemeName != "steadyState" && d2dt2SchemeName == "Euler")
                    {
                        matrix(3*i, 3*i) += d2dt2WDiag[i];
                        matrix(3*i + 1, 3*i + 1) += d2dt2ThetaXDiag[i];
                        matrix(3*i + 2, 3*i + 2) += d2dt2ThetaYDiag[i];

                        source[3*i] += d2dt2WSource[i];
                        source[3*i + 1] += d2dt2ThetaXSource[i];
                        source[3*i + 2] += d2dt2ThetaYSource[i];
                    }
                    else
                    {
                        FatalError("evolve() function in Mindlin plate solid") << nl
                            << "Incompatible (or not defined) d2dt2Scheme "
                            << d2dt2SchemeName << " is specified! "
                            << abort(FatalError);
                    }
                }

                forAll(lapWUpper, edgeI)
                {
                    const label i = own[edgeI];
                    const label j = nei[edgeI];

                    // 1. Off-diagonal components of the laplacian terms
                    // Upper part of the matrix
                    matrix(3*i, 3*j) = lapWUpper[i];
                    matrix(3*i + 1, 3*j + 1) = lapThetaXUpper[i];
                    matrix(3*i + 2, 3*j + 2) = lapThetaYUpper[i];

                    // Lower part of the matrix
                    matrix(3*j, 3*i) = lapWUpper[i];
                    matrix(3*j + 1, 3*i + 1) = lapThetaXUpper[i];
                    matrix(3*j + 2, 3*i + 2) = lapThetaYUpper[i];

                    // Information reqd to compute fam::div(theta)
                    const scalar wt = interpWeights.internalField()[edgeI];
                    const scalar leI = le.internalField()[edgeI];
                    const scalar nxEdge = nx[edgeI];
                    const scalar nyEdge = ny[edgeI];
                    const scalar Gamma = shearStrainStiffness_.value();
                    const scalar delta = deltaCoeffs.internalField()[edgeI];


                    // 2. fam::div(theta) implicit contributions
                    const scalar thXCoeff = Gamma*nxEdge*leI;
                    const scalar thYCoeff = Gamma*nyEdge*leI;

                    const vector dROwn =
                        edgeCentres[edgeI] - cellCentres[own[edgeI]];
                    const vector dRNei =
                        edgeCentres[edgeI] - cellCentres[nei[edgeI]];

                    // 2. (-) Gamma \int_l fam::div(theta) dl
                    // contribution to wEqn
                    matrix(3*i, 3*i + 1) -= wt*thXCoeff;
                    matrix(3*i, 3*i + 2) -= wt*thYCoeff;

                    matrix(3*i, 3*j + 1) -= (1 - wt)*thXCoeff;
                    matrix(3*i, 3*j + 2) -= (1 - wt)*thYCoeff;

                    matrix(3*j, 3*i + 1) += wt*thXCoeff;
                    matrix(3*j, 3*i + 2) += wt*thYCoeff;

                    matrix(3*j, 3*j + 1) += (1 - wt)*thXCoeff;
                    matrix(3*j, 3*j + 2) += (1 - wt)*thYCoeff;

                    // 3a. (-) Gamma \int_l (x - x_P)*fam::div(theta) dl
                    // contribution to the thetaX Eqn
                    matrix(3*i + 1, 3*i + 1) -= wt*thXCoeff*dROwn.x();
                    matrix(3*i + 1, 3*i + 2) -= wt*thYCoeff*dROwn.x();

                    matrix(3*i + 1, 3*j + 1) -= (1 - wt)*thXCoeff*dROwn.x();
                    matrix(3*i + 1, 3*j + 2) -= (1 - wt)*thYCoeff*dROwn.x();

                    matrix(3*j + 1, 3*i + 1) += wt*thXCoeff*dRNei.x();
                    matrix(3*j + 1, 3*i + 2) += wt*thYCoeff*dRNei.x();

                    matrix(3*j + 1, 3*j + 1) += (1 - wt)*thXCoeff*dRNei.x();
                    matrix(3*j + 1, 3*j + 2) += (1 - wt)*thYCoeff*dRNei.x();

                    // 3b. (-) \int_l Gamma (y - y_P)*fam::div(theta) dl
                    // contribution to the thetaY Eqn
                    matrix(3*i + 2, 3*i + 1) -= wt*thXCoeff*dROwn.y();
                    matrix(3*i + 2, 3*i + 2) -= wt*thYCoeff*dROwn.y();

                    matrix(3*i + 2, 3*j + 1) -= (1 - wt)*thXCoeff*dROwn.y();
                    matrix(3*i + 2, 3*j + 2) -= (1 - wt)*thYCoeff*dROwn.y();

                    matrix(3*j + 2, 3*i + 1) += wt*thXCoeff*dRNei.y();
                    matrix(3*j + 2, 3*i + 2) += wt*thYCoeff*dRNei.y();

                    matrix(3*j + 2, 3*j + 1) += (1 - wt)*thXCoeff*dRNei.y();
                    matrix(3*j + 2, 3*j + 2) += (1 - wt)*thYCoeff*dRNei.y();


                    // 4. \int_l Gamma (x - x_P)*fam::grad(w) \dot n dl
                    //  contribution to thetaX and thetaY
                    matrix(3*i + 1, 3*i) -= Gamma*leI*delta*dROwn.x();
                    matrix(3*i + 2, 3*i) -= Gamma*leI*delta*dROwn.y();

                    matrix(3*i + 1, 3*j) += Gamma*leI*delta*dROwn.x();
                    matrix(3*i + 2, 3*j) += Gamma*leI*delta*dROwn.y();

                    matrix(3*j + 1, 3*i) -= Gamma*leI*delta*dROwn.x();
                    matrix(3*j + 2, 3*i) -= Gamma*leI*delta*dROwn.y();

                    matrix(3*j + 1, 3*j) += Gamma*leI*delta*dROwn.x();
                    matrix(3*j + 2, 3*j) += Gamma*leI*delta*dROwn.y();

                    // 2. (-) \int_l fam::div(theta) dl contribution to wEqn
                    // // Off-diagonal of owner (diagonal) cell block
                    // // thetaX and thetaY contribution in wEqn (fam::div(theta))
                    // matrix(3*i, 3*i + 1) -= Gamma*((wt * nxEdge) * leI);
                    // matrix(3*i, 3*i + 2) -= Gamma*((wt * nyEdge) * leI);

                    // // Off-diagonal of owner cell block
                    // // thetaY contribution in wEqn (fam::div(theta))
                    // matrix(3*i, 3*j + 1) -= Gamma*(((1 - wt) * nxEdge) * leI);
                    // matrix(3*i, 3*j + 2) -= Gamma*(((1 - wt) * nyEdge) * leI);

                    // // Off-diagonal of neighbour cell block
                    // // thetaX contribution in wEqn (fam::div(theta))
                    // matrix(3*j, 3*i + 1) += Gamma*((wt * nxEdge) * leI);
                    // matrix(3*j, 3*i + 2) += Gamma*((wt * nyEdge) * leI);

                    // // Diagonal of neighbour cell block
                    // // thetaX and thetaY contribution in wEqn (fam::div(theta))
                    // matrix(3*j, 3*j + 1) += Gamma*(((1 - wt) * nxEdge) * leI);
                    // matrix(3*j, 3*j + 2) += Gamma*(((1 - wt) * nyEdge) * leI);


                    // // 3. Implicit terms of shear force in thetaX and thetaY Eqns
                    // // Owner cell (thetaX)
                    // 3a. (-) \int_l (x - x_P)*fam::div(theta) dl contribution
                    // to the thetaX Eqn
                    // matrix(3*i + 1, 3*i + 1) -= Gamma*leI*wt*nxEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[own[edgeI]].x());
                    // matrix(3*i + 1, 3*i + 2) -= Gamma*leI*wt*nyEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[own[edgeI]].x());

                    // matrix(3*i + 1, 3*j + 1) -= Gamma*leI*(1 - wt)*nxEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[own[edgeI]].x());
                    // matrix(3*i + 1, 3*j + 2) -= Gamma*leI*(1 - wt)*nyEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[own[edgeI]].x());


                    // // Neighbour cell (thetaX)
                    // 3b. (-) \int_l (y - y_P)*fam::div(theta) dl contribution
                    // to the thetaY Eqn
                    // matrix(3*j + 1, 3*i + 1) += Gamma*leI*wt*nxEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[nei[edgeI]].x());
                    // matrix(3*j + 1, 3*i + 2) += Gamma*leI*wt*nyEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[nei[edgeI]].x());

                    // matrix(3*j + 1, 3*j + 1) += Gamma*leI*(1 - wt)*nxEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[nei[edgeI]].x());
                    // matrix(3*j + 1, 3*j + 2) += Gamma*leI*(1 - wt)*nyEdge
                    //     *(edgeCentres[edgeI].x() - cellCentres[nei[edgeI]].x());

                    // // Owner cell (thetaY)
                    // matrix(3*i + 2, 3*i + 1) -= Gamma*leI*wt*nxEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[own[edgeI]].y());
                    // matrix(3*i + 2, 3*i + 2) -= Gamma*leI*wt*nyEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[own[edgeI]].y());

                    // matrix(3*i + 2, 3*j + 1) -= Gamma*leI*(1 - wt)*nxEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[own[edgeI]].y());
                    // matrix(3*i + 2, 3*j + 2) -= Gamma*leI*(1 - wt)*nyEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[own[edgeI]].y());


                    // // Neighbour cell (thetaY)
                    // matrix(3*j + 2, 3*i + 1) += Gamma*leI*wt*nxEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[nei[edgeI]].y());
                    // matrix(3*j + 2, 3*i + 2) += Gamma*leI*wt*nyEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[nei[edgeI]].y());

                    // matrix(3*j + 2, 3*j + 1) += Gamma*leI*(1 - wt)*nxEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[nei[edgeI]].y());
                    // matrix(3*j + 2, 3*j + 2) += Gamma*leI*(1 - wt)*nyEdge
                    //     *(edgeCentres[edgeI].y() - cellCentres[nei[edgeI]].y());

                }
                // Loop over boundary patches
                forAll(w_.boundaryField(), patchI)
                {
                    const word& patchNameW(w_.boundaryField()[patchI].type());
                    const word& patchNameThX(thetaX_.boundaryField()[patchI].type());
                    const word& patchNameThY(thetaY_.boundaryField()[patchI].type());

                    const List<scalar>& pDelta(aMesh_.boundary()[patchI].deltaCoeffs());
                    const scalar Gamma = shearStrainStiffness_.value();
                    const vectorField pEdgeCentres(edgeCentres.boundaryField()[patchI]);

                    // This is wrong as pCellCentres and pEdgeCentres at boundaries are same
                    // const vectorField pCellCentres(cellCentres.boundaryField()[patchI]);

                    const scalarField nxb
                    (
                        edgeBiNormal.boundaryField()[patchI].component(vector::X)
                    );

                    const scalarField nyb
                    (
                        edgeBiNormal.boundaryField()[patchI].component(vector::Y)
                    );

                    // Procedure to grab the IDs of cells attached to boundary edges
                    // The information on owners of boundary edges is stored in
                    // aMesh_.edgeOwners(). The goal here is to extract the labels
                    // of individual boundary patches.

                    // const label nBouEdges = pEdgeCentres.size();
                    // const label nIntEdges = aMesh_.internalEdges().size();
                    // labelList cellCentreBou(nBouEdges, 0.0);
                    // const label startIndex(nIntEdges + patchI*nBouEdges);

                    // for (label i = 0; i <= nBouEdges; ++i)
                    // {
                    //     cellCentreBou[i] = edgeOwn[startIndex + i];
                    // }

                    // Loop over all faces of boundary patch
                    forAll(w_.boundaryField()[patchI], faceI)
                    {
                        // Boundary cell index
                        const label bI = faBouMesh[patchI].edgeFaces()[faceI];
                        const scalar leB = le.boundaryField()[patchI][faceI];
                        const scalar delB = pDelta[faceI];
                        // const label cellID = cellCentreBou[faceI];

                        const vector pDr = (pEdgeCentres[faceI] - cellCentres[bI]);
                        const scalar pThXCoeff = Gamma*nxb[faceI]*leB;
                        const scalar pThYCoeff = Gamma*nyb[faceI]*leB;

                        // Contribution of boundary edges to the diagonal of the matrix
                        matrix(3*bI, 3*bI) += lapWIntCoeffs[patchI][faceI];
                        matrix(3*bI + 1, 3*bI + 1) += lapThetaXIntCoeffs[patchI][faceI];
                        matrix(3*bI + 2, 3*bI + 2) += lapThetaYIntCoeffs[patchI][faceI];

                        // Explicit contribution of boundary edges to the source
                        source[3*bI] += lapWBouCoeffs[patchI][faceI];
                        source[3*bI + 1] += lapThetaXBouCoeffs[patchI][faceI];
                        source[3*bI + 2] += lapThetaYBouCoeffs[patchI][faceI];

                        if(patchNameW == "fixedValue")
                        {
                            // Boundary Coefficients of fam::grad(w)
                            matrix(3*bI + 1, 3*bI) -= Gamma*leB*delB*pDr.x();
                            matrix(3*bI + 2, 3*bI) -= Gamma*leB*delB*pDr.y();
                        }
                        else
                        {
                            FatalErrorIn
                            (
                                "bool mindlinDemirdzicPlateSolid::evolve()"
                            )   << "Block-Coupled Approach not implemented for BC: "
                                << patchNameW
                                << abort(FatalError);
                        }

                        if (patchNameThX == "fixedValue" && patchNameThY == "fixedValue")
                        {
                        // For clamped ends all sides, thetaX and thetaY are
                        // zero at BC. Hence zero contribution
                        }
                        else if (patchNameThX == "zeroGradient" && patchNameThY == "zeroGradient")
                        {
                            matrix(3*bI + 1, 3*bI + 1) -= pThXCoeff*pDr.x();
                            matrix(3*bI + 1, 3*bI + 2) -= pThYCoeff*pDr.x();

                            matrix(3*bI + 2, 3*bI + 1) -= pThXCoeff*pDr.y();
                            matrix(3*bI + 2, 3*bI + 2) -= pThYCoeff*pDr.y();
                        }
                        else
                        {
                            FatalErrorIn
                            (
                                "bool mindlinDemirdzicPlateSolid::evolve()"
                            )   << "Block-Coupled Approach not implemented for thetaX BC: "
                                << patchNameThX << " and thetaY BC: "
                                << patchNameThY
                                << abort(FatalError);
                        }

                    }
                }

                if(debug > 1)
                {
                    Info<< "\nBlock Matrix Coefficients: " << matrix.data() << endl;
                    Info<< "\nSource vector: " << source << endl;
                }

                // Solve the linear system of equations
                // Using Eigen SparseLU direct solver
                sparseMatrixTools::solveLinearSystemEigen
                (
                    matrix, source, solveWTheta, false, debug
                );

                // Retrieve solution
                for(label i = 0; i < nCells; ++i)
                {
                    w_[i] = solveWTheta[3*i];
                    thetaX_[i] = solveWTheta[3*i + 1];
                    thetaY_[i] = solveWTheta[3*i + 2];
                }

                // Correct the boundary conditions for w, thetaX and thetaY
                w_.correctBoundaryConditions();
                thetaX_.correctBoundaryConditions();
                thetaY_.correctBoundaryConditions();
            }
            while
            (
                !blockConverged
                (
                    iCorr,
                    w_,
                    thetaX_,
                    thetaY_
                )
                &&
                ++iCorr < nCorr()
            );
        }
        else
        {

            // const scalar scaleW(readScalar(solidModelDict().lookup("scaleW")));

            // Approach 2: Using segregated method of solving w, thetaX and thetaY equations separately
            // and then iteratively update them until the values fall below a solution tolerance

            Info<< "\nUsing segregated approach to solve for w, thetaX and thetaY eqns "
                << "separately and iteratively update them!" << endl;

            do
            {
                // w_ = scaleW*wAnalytical_;

                theta = thetaX_*vector(1, 0, 0) + thetaY_*vector(0, 1, 0);

                // Info<< "theta " << theta << endl;
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

                // w_.storeOldTime();
                // thetaX_.storeOldTime();
                // thetaY_.storeOldTime();
                // gradW_.storeOldTime();
                // gradThetaX_.storeOldTime();
                // gradThetaY_.storeOldTime();

                // Solve w equation
                // Also, "==" complains so we will move all terms to left
                // QUESTION - Is fac::div(shearStrainStiffness_*theta)
                // equivalent to physically looping over edges and putting
                // theta contributions interpolated at the edges in the wEqn?
                // I think it is correct because that is how fac::div
                // is calculated, but need to CHECK WITH IVAN.
                faScalarMatrix wEqn
                (
                    fam::laplacian(shearStrainStiffness_, w_)
                  - fac::div(shearStrainStiffness_*theta)
                  - p_
                );

                // d2dt2 can only take Euler as keyword, but if the user wants it to
                // be steadyState, it cannot happen. Hence check for ddtScheme
                // and add inertial terms for not steady state!!
                if(ddtSchemeName != "steadyState")
                {
                    wEqn -= rho_*h_*(fac::d2dt2(w_));
                }

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

                // Relax the linear system
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
                  - fam::Sp(shearStrainStiffness_, thetaX_)
                  + shearStrainStiffness_*thetaX_
                    // The fam:Sp terms above needs to be added implictly and
                    // removed explictly to increase diagonal dominance
                    // Otherwise the solver wont converge
                    // This trick in not mentioned in Demirdzic 1997 paper.
                );

                // Initialise thetaY equation with implicit laplacian terms
                faScalarMatrix thetaYEqn
                (
                    fam::laplacian(bendingStiffness_, thetaY_)
                  - fam::Sp(shearStrainStiffness_, thetaY_)
                  + shearStrainStiffness_*thetaY_
                    // Same comment as thetaX Eqn
                );

                // d2dt2 can only take Euler as keyword, but if the user wants it to
                // be steadyState, it cannot happen. Hence check for ddtScheme
                // and add inertial terms for not steady state!!
                if(ddtSchemeName != "steadyState")
                {
                    thetaXEqn -= rho_*pow(h_,3)*(fac::d2dt2(thetaX_))/12;
                    thetaYEqn -= rho_*pow(h_,3)*(fac::d2dt2(thetaY_))/12;
                }

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

                /*---------------------------------------------------------------*/
                /*---------------------------------------------------------------*/
                edgeVectorField gradWEdge(fac::interpolate(gradW_));

                // Avoid oscillations in gradient calculations
                if (compactEdgeNormalGrad)
                {
                    const edgeScalarField lnGradWEdge(fac::lnGrad(w_));
                    gradWEdge +=
                        lnGradWEdge*edgeBiNormal - (sqr(edgeBiNormal) & gradWEdge);
                }

                // const vectorField gradWEdgeI(gradWEdge.internalField());
                // const scalarField leI(le.internalField());

                // APPROACH 2 (LOOP OVER CELLS and for each cell LOOP over EDGES)
                const edgeVectorField shearForceEdge
                (
                    shearStrainStiffness_*(gradWEdge - thetaEdge)
                );

                // Function call to evaluate the shear force contribution to source
                const vectorField shearForceContrib =
                    calculateShearForceContribution(shearForceEdge);

                thetaXEqn.source() += shearForceContrib.component(vector::X);
                thetaYEqn.source() += shearForceContrib.component(vector::Y);

                /*---------------------------------------------------------------*/
                /*---------------------------------------------------------------*/
                // GRADIENT THETA TERMS OF THE SOURCE
                // The first three terms of s_{\phi \l} term of thetaX and thetaY variables
                // in Demirdzic 1997 plate paper.

                // Interpolate grad(thetaX), grad(thetaY) to edges
                edgeVectorField gradthetaXEdge(fac::interpolate(gradThetaX_));
                edgeVectorField gradthetaYEdge(fac::interpolate(gradThetaY_));

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

                // Function call to evaluate the gradient theta contributions to source
                const vectorField gradThetaContrib =
                    calculateGradientThetaContribution
                    (
                        gradthetaXEdge,
                        gradthetaYEdge
                    );

                thetaXEqn.source() -= gradThetaContrib.component(vector::X);
                thetaYEqn.source() -= gradThetaContrib.component(vector::Y);

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
        }


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
