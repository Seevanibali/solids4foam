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

Class
    Foam::simplySupportedMomentFaPatchField

Description
    Free plate edge boundary condition for transverse displacement (w)
    - Demirdzic Mindlin plate

Author
    Seevani Bali, UCD.

\*---------------------------------------------------------------------------*/

#ifndef freeEdgeDisplacementDemirdzicFaPatchField_H
#define freeEdgeDisplacementDemirdzicFaPatchField_H

#include "fixedGradientFaPatchField.H"
#include "areaFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class freeEdgeDisplacementDemirdzicFaPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class freeEdgeDisplacementDemirdzicFaPatchScalarField
:
    public fixedGradientFaPatchField<scalar>
{
    // Private data

        //- Shear force
        //const scalarField shearForce_;

        //- Relaxation factor
        const scalar relaxFac_;


public:

    //- Runtime type information
    TypeName("freeEdgeDisplacementDemirdzic");


    // Constructors

        //- Construct from patch and internal field
        freeEdgeDisplacementDemirdzicFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct from patch, internal field and dictionary
        freeEdgeDisplacementDemirdzicFaPatchScalarField
        (
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const dictionary&
        );

        //- Construct by mapping given patch field onto a new patch
        freeEdgeDisplacementDemirdzicFaPatchScalarField
        (
            const freeEdgeDisplacementDemirdzicFaPatchScalarField&,
            const faPatch&,
            const DimensionedField<scalar, areaMesh>&,
            const faPatchFieldMapper&
        );

        //- Construct as copy
        freeEdgeDisplacementDemirdzicFaPatchScalarField
        (
            const freeEdgeDisplacementDemirdzicFaPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<faPatchField<scalar> > clone() const
        {
            return tmp<faPatchField<scalar> >
            (
                new freeEdgeDisplacementDemirdzicFaPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        freeEdgeDisplacementDemirdzicFaPatchScalarField
        (
            const freeEdgeDisplacementDemirdzicFaPatchScalarField&,
            const DimensionedField<scalar, areaMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<faPatchField<scalar> > clone
        (
            const DimensionedField<scalar, areaMesh>& iF
        ) const
        {
            return tmp<faPatchField<scalar> >
            (
                new freeEdgeDisplacementDemirdzicFaPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
