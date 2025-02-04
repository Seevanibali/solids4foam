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

#include "freeEdgeDisplacementFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "kirchhoffPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(dict.lookupOrDefault<scalar>("relaxationFactor", 1.0))

{
   if (dict.found("value"))
   {
       faPatchField<scalar>::operator==(Field<scalar>("value", dict, p.size()));
   }
   else
   {
       faPatchField<scalar>::operator==(0.0);
   }

   gradient() = 0.0;
}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf
)
:
    fixedGradientFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementFaPatchScalarField::
freeEdgeDisplacementFaPatchScalarField
(
    const freeEdgeDisplacementFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeEdgeDisplacementFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Info<< nl << "------------------------------------------" << nl << endl;
    // Lookup angle of rotation field
    const faPatchField<vector>& theta =
        patch().lookupPatchField<areaVectorField, vector>("theta");

    // Lookup the gradient of rotation field
    const faPatchField<tensor>& gradTheta =
        patch().lookupPatchField<areaTensorField, tensor>("grad(theta)");

    // Calculate correction vectors
    const vectorField n(patch().edgeNormals());
    const vectorField delta(patch().delta());
    const vectorField k((I - sqr(n)) & delta);

    // Calculate the patch internal field and correction for non-orthogonality
    const scalarField nDotThetaPif
    (
        n & (theta.patchInternalField() + (k & gradTheta.patchInternalField()))
    );

    // Lookup angle of rotation field
    // const faPatchField<vector>& theta =
    //     patch().lookupPatchField<areaVectorField, vector>("theta");

    // if (!db().foundObject<areaVectorField>("thetaPrevIter"))
    // {
    //     Info<< "theta pre iter not found " << endl;
    //     return;
    // }

    
    // Lookup angle of rotation field
    // const faPatchField<vector>& thetaPrevIter =
    //     patch().lookupPatchField<areaVectorField, vector>("thetaPrevIter");

    // Lookup the gradient of rotation field
    // const faPatchField<tensor>& gradTheta =
    //     patch().lookupPatchField<areaTensorField, tensor>("grad(theta)");

    // // Calculate correction vectors
    // const vectorField n(patch().edgeNormals());
    // const vectorField delta(patch().delta());
    // const vectorField k = (I - sqr(n)) & delta;

    // // Calculate the patch internal field and correction for non-orthogonality
    // const scalarField nDotThetaPif =
    //  n & (theta.patchInternalField() + (k & gradTheta.patchInternalField()));

    // Info<< "theta pif: "
    //     << (theta.patchInternalField()) << nl << endl;

    // gradTheta patch internal field in normal direction
    // const vectorField nDotGradThetaPif = n & gradTheta.patchInternalField();

    // Lookup moment sum field
    const faPatchField<scalar>& M =
        patch().lookupPatchField<areaScalarField, scalar>("M");

    // Ivan's code
    const scalarField MI(M.patchInternalField());

    // const faPatchField<scalar>& MPrevIter =
    //     patch().lookupPatchField<areaScalarField, scalar>("MPrevIter");

    // // Lookup w field
    // const faPatchField<scalar>& w =
    //     patch().lookupPatchField<areaScalarField, scalar>("w");

    // const faPatchField<scalar>& wPrevIter =
    //     patch().lookupPatchField<areaScalarField, scalar>("wPrevIter");

    // Info<< "theta: " 
    //     << vectorField(theta) 
    //     << nl << endl;


    // Lookup fvMesh
    // Fix for FSI: this is only correct for
    const fvMesh* vmeshPtr = NULL;
    if (db().parent().foundObject<fvMesh>("solid"))
    {
        vmeshPtr = &db().parent().lookupObject<fvMesh>("solid");
    }
    else
    {
        vmeshPtr = &db().parent().lookupObject<fvMesh>("region0");
    }
    const fvMesh& vmesh = *vmeshPtr;

    // Lookup the solidModel object
    const solidModel& solMod = lookupSolidModel(vmesh);

    // Cast the solid model to a Kirchhoff plate solid model
    const solidModels::kirchhoffPlateSolid& plateSolid =
        refCast<const solidModels::kirchhoffPlateSolid>(solMod);

    // Lookup flexural stiffness from the solver
    const scalar D = plateSolid.bendingStiffness().value();
    const scalar nu = plateSolid.nu().value();

    // Update the gradient
    // gradient() +=
    //     relaxFac_*(1.0 + nu)*M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs());

    // Info<< "n: " << n << nl << endl;
    // Info<< "delta: " << delta << nl << endl;

    // gradient() = (n & theta)
    //              // + (1.0 - relaxFac_)*(n & thetaPrevIter) 
    //             + M*nu/(2*D*(1.0 - pow(nu,2))*(n & delta));
    //             + (1.0 - relaxFac_)*MPrevIter*nu/(2*D*(1.0 - nu)*(n & delta));

    // Info<< "g aft: " << gradient() << nl << endl;         
    //  gradient() =  relaxFac_*(n & -theta.patchInternalField())
    //             + (1.0 - relaxFac_)*(n & -thetaPrevIter.patchInternalField())
    //             + relaxFac_*M*nu/(2*D*(1 - nu)*mag(delta))
    //             + (1.0 - relaxFac_)*MPrevIter*nu/(2*D*(1 - nu)*mag(delta));

// Ivan's code
//     const scalarField A(n & theta);
//     const scalarField B((0.5*nu/(1.0-nu))*(MI/D)*(1/patch().deltaCoeffs()));

// //Ivan's code
// gradient() =
//         relaxFac_*(-A+B) + (1.0 - relaxFac_)*gradient();

// Seevani's version

    // Seevani's original code
    // gradient() = relaxFac_
    //             *(
    //                 (n & -theta.patchInternalField()) + M*nu/(2*D*(1 - nu)*mag(delta))
    //              )
    //             + (1.0 - relaxFac_)*gradient();

    // Seevani comparison with Ivan - 29th Oct

    const scalarField AA(n & theta.patchInternalField());

    // Correct version
    const scalarField BB((nu/(1.0-nu))*(M.patchInternalField()/D)*(1/patch().deltaCoeffs()));
    
    // Incorrect version
    // const scalarField BB((nu*M*mag(delta)/((1.0-nu)*D)));



     gradient() = relaxFac_
                *(-AA + BB)
                + (1.0 - relaxFac_)*gradient();



    // Info<< "gradient: " << gradient() << nl << endl;

    // gradient() =
    //     relaxFac_
    //    *(
    //        nDotThetaPif
    //     //  + (delta & nDotGradThetaPif)
    //      + M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs())
    //     )
    //   + (1.0 - relaxFac_)*gradient();

    // Info<< "new gradient(): " << gradient() << nl << endl;

    // Info<< "this (w): " << max(scalarField(*this)) << nl << endl;

    // Info<< "nDotThetaPif: " << nDotThetaPif << endl;
    // Info<< "(delta & nDotGradThetaPif) "
    //     << (delta & nDotGradThetaPif) << endl;
    // Info<< "3: "
    //     << M/(D*(1.0 - pow(nu, 2))*patch().deltaCoeffs()) << endl;

    fixedGradientFaPatchField<scalar>::updateCoeffs();

    //Info<< nl << "-------------------------------------------" << nl << endl;
}


void Foam::freeEdgeDisplacementFaPatchScalarField::write
(
    Ostream& os
) const
{
    faPatchField<scalar>::write(os);

    os.writeKeyword("relaxationFactor")
        << relaxFac_ << token::END_STATEMENT << endl;

    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFaPatchTypeField
    (
        faPatchScalarField,
        freeEdgeDisplacementFaPatchScalarField
    );
}

#endif // OPENFOAM_ORG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
