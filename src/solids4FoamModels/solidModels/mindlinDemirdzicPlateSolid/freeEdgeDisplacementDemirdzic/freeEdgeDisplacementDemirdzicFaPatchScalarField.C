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

#include "freeEdgeDisplacementDemirdzicFaPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "mindlinDemirdzicPlateSolid.H"
#include "lookupSolidModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::
freeEdgeDisplacementDemirdzicFaPatchScalarField
(
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(p, iF),
    relaxFac_(1.0)
{}


Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::
freeEdgeDisplacementDemirdzicFaPatchScalarField
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


Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::
freeEdgeDisplacementDemirdzicFaPatchScalarField
(
    const freeEdgeDisplacementDemirdzicFaPatchScalarField& ptf,
    const faPatch& p,
    const DimensionedField<scalar, areaMesh>& iF,
    const faPatchFieldMapper& mapper
)
:
    fixedGradientFaPatchField<scalar>(ptf, p, iF, mapper),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::
freeEdgeDisplacementDemirdzicFaPatchScalarField
(
    const freeEdgeDisplacementDemirdzicFaPatchScalarField& ptf
)
:
    fixedGradientFaPatchField<scalar>(ptf),
    relaxFac_(ptf.relaxFac_)
{}


Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::
freeEdgeDisplacementDemirdzicFaPatchScalarField
(
    const freeEdgeDisplacementDemirdzicFaPatchScalarField& ptf,
    const DimensionedField<scalar, areaMesh>& iF
)
:
    fixedGradientFaPatchField<scalar>(ptf, iF),
    relaxFac_(ptf.relaxFac_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    // Info<< nl << "------------------------------------------" << nl << endl;

    // Lookup angle of rotation field
    const faPatchField<scalar>& thetaX =
        patch().lookupPatchField<areaScalarField, scalar>("thetaX");

    // const faPatchField<vector>& gradW =
    // patch().lookupPatchField<areaVectorField, vector>("grad(w)");

    // const vectorField gradWI(gradW.patchInternalField());
    // Info<< "gradW " << gradW << nl << "gradW I " << gradWI << endl;
    // Lookup angle of rotation field
    // const faPatchField<scalar>& thetaY =
    //     patch().lookupPatchField<areaScalarField, scalar>("thetaY");

    // Patch edge binormal
    const vectorField n(patch().edgeLengths()/patch().magEdgeLengths());

    // Since the solver is 1-D with empty patches in y-direction,
    // gradient() = n & theta, becomes simply thetaX

    gradient() = thetaX;
    // gradient() = (gradWI & n);
    Info<< "gradient of w " << gradient() << endl;

    fixedGradientFaPatchField<scalar>::updateCoeffs();
}


void Foam::freeEdgeDisplacementDemirdzicFaPatchScalarField::write
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
        freeEdgeDisplacementDemirdzicFaPatchScalarField
    );
}

#endif // OPENFOAM_ORG

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
