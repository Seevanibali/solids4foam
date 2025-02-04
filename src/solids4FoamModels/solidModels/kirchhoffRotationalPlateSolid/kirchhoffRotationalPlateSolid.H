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
    kirchhoffRotationalPlateSolid

Description
    Rotational Kirchhoff plate solver.

    Based on "rotational formulation" of Torlak (2006).

Author
    Philip Cardiff, UCD.  All rights reserved.
    Seevani Bali, UCD.  All rights reserved.

SourceFiles
    kirchhoffRotationalPlateSolid.C

\*---------------------------------------------------------------------------*/

#ifndef kirchhoffRotationalPlateSolid_H
#define kirchhoffRotationalPlateSolid_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "areaFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
                        Class kirchhoffRotationalPlateSolid Declaration
\*---------------------------------------------------------------------------*/

class kirchhoffRotationalPlateSolid
:
    public solidModel
{
    // Private data

        //- Finite Area mesh
        faMesh aMesh_;

        //- Transverse displacement
        //  displacement prependicular to the plane of the plate
        areaScalarField w_;

        //- Transverse displacement volume field for visualisation
        volScalarField wVf_;

        //- Moment sum
        areaScalarField M_;

        //- Moment sum volume field for visualisation
        volScalarField MVf_;

        //- Gradient of Moment sum
        areaVectorField gradM_;

        //- Net pressure field on the plate (perpendicular to the plate normal)
        areaScalarField p_;

        //- Net pressure field on the plate volume field
        volScalarField pVf_;
        
        //- Angle of rotation thetaX (about y - axis)
        areaScalarField thetaX_;

        //- Angle of rotation volume field for visualisation
        volScalarField thetaXVf_;

        //- Angle of rotation thetaY (about x - axis)
        areaScalarField thetaY_;

        //- Angle of rotation volume field for visualisation
        volScalarField thetaYVf_;

        //- Zero scalar Field required to set the z-comp of theta to zero
        areaScalarField zeroField_;
        
        // //- Angle of rotation 
        // areaVectorField theta_;

        // //- Angle of rotation volume field for visualisation
        // volVectorField thetaVf_;

        // //- Gradient of angle of rotation
        // areaTensorField gradTheta_;

        //- Gradient of angle of rotation thetaX
        areaVectorField gradThetaX_;

        //- Gradient of angle of rotation thetaY
        areaVectorField gradThetaY_;

        //- Gradient of thetaX - x component
        areaScalarField gradThXX_;

        //- Gradient of thetaX - y component
        areaScalarField gradThXY_;

        //- Gradient of thetaY - x component
        areaScalarField gradThYX_;

        //- Gradient of thetaY - y component
        areaScalarField gradThYY_;

        //- Source vector Q (x-component) in thetaX Eqn 
        areaScalarField QThetaXX_;

        //- Source vector Q (y-component) in thetaX Eqn 
        areaScalarField QThetaXY_;

        //- Source vector Q (x-component) in thetaY Eqn 
        areaScalarField QThetaYX_;

        //- Source vector Q (y-component) in thetaY Eqn 
        areaScalarField QThetaYY_;

        //- Source vector Q in thetaX Eqn 
        areaVectorField QThetaX_;

        //- Source vector Q in thetaY Eqn 
        areaVectorField QThetaY_;

        //- Plate density
        dimensionedScalar rho_;

        //- Plate Young's modulus
        dimensionedScalar E_;

        //- Plate Poisson's ratio
        dimensionedScalar nu_;

        //- Plate thickness
        const dimensionedScalar h_;

        //- Plate bending stiffness
        dimensionedScalar bendingStiffness_;

        //- Plate torsional stiffness
        dimensionedScalar torsionalStiffness_;

        //- areaPatch index
        //  This is the patch on the fvMesh mesh that corresponds to the
        //  areaMesh
        //  Note: it is assumed that this is one patch
        mutable label areaPatchID_;

        //- areaShadowPatch index
        //  This is the patch on the fvMesh mesh that is directly opposite from
        //  the areaPatch
        //  Note: it is assumed that the fvMesh is one cell thick from the
        //  areaPatch to the areaShadowPatch
        mutable label areaShadowPatchID_;

        //- Boolean for block-coupled -- SB
        const bool coupled_;

        //- Debug flag
        const bool debug_;

    // Private Member Functions

        //- Check if the momentum equation is converged
        bool converged
        (
        const int iCorr,
    #ifdef OPENFOAM_NOT_EXTEND
        const SolverPerformance<scalar>& solverPerfM,
        const SolverPerformance<scalar>& solverPerfThetaX,
        const SolverPerformance<scalar>&  solverPerfThetaY,
    #else
        const lduSolverPerformance& solverPerfM,
        const lduSolverPerformance& solverPerfThetaX,
        const lduSolverPerformance& solverPerfThetaY,
    #endif
        const areaScalarField& M,
        const areaScalarField& thetaX,
        const areaScalarField& thetaY,
        const areaScalarField& w
        );

        //- Return a reference to the areaPatch
        const fvPatch& areaPatch() const;

        //- Return a reference to the areaShadowPatch
        const fvPatch& areaShadowPatch() const;

        //- Calculate and set areaPatchID and shadowPatchID
        void calcAreaPatches() const;

        //- Map an area geometric field to a volume geometric field, where we
        //  assume the volume mesh to be one layer thick i.e. there is a
        //  one-to-one relationship between the number of area faces to the
        //  number of volume cells
        template<class Type>
        void mapAreaFieldToSingleLayerVolumeField
        (
            const GeometricField<Type, faPatchField, areaMesh>& af,
            GeometricField<Type, fvPatchField, volMesh>& vf
        ) const;

        //- Disallow default bitwise copy construct
        kirchhoffRotationalPlateSolid(const kirchhoffRotationalPlateSolid&);

        //- Disallow default bitwise assignment
        void operator=(const kirchhoffRotationalPlateSolid&);


protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("kirchhoffRotationalPlate");

    // Constructors

        //- Construct from components
        kirchhoffRotationalPlateSolid
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~kirchhoffRotationalPlateSolid()
        {}


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented(type());

                // Keep compiler happy
                return D();
            }

            //- Plate bending stiffness
            const dimensionedScalar& bendingStiffness() const
            {
                return bendingStiffness_;
            }

            //- Plate Poisson's ration
            const dimensionedScalar& nu() const
            {
                return nu_;
            }


        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const;

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "kirchhoffRotationalPlateSolidTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
