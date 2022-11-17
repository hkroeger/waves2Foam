/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "relaxationShapePatchDist.H"
#include "addToRunTimeSelectionTable.H"

#include "patchDistMethod.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace relaxationShapes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(relaxationShapePatchDist, 0);
addToRunTimeSelectionTable
(
    relaxationShape,
    relaxationShapePatchDist,
    dictionary
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


relaxationShapePatchDist::relaxationShapePatchDist
(
    const word& subDictName,
    const fvMesh& mesh_
)
:
    relaxationShape(subDictName, mesh_),
    relaxationPatches_(coeffDict_.lookup("relaxationPatches")),
      y_(
          IOobject
          (
              "yRelax",
              mesh_.time().timeName(),
              mesh_
          ),
          mesh_,
          dimensionedScalar("yRelax", dimLength, SMALL),
          patchDistMethod::patchTypes<scalar>(
              mesh_,
              mesh_.boundaryMesh().patchSet(
                  relaxationPatches_
                  ))
      ),
      width_(readScalar(coeffDict_.lookup("width")))
{
    
    ym_.reset(patchDistMethod::New(
                 coeffDict_.subDict("patchDist"),
                 mesh_,
                 mesh_.boundaryMesh().patchSet(
                     relaxationPatches_
                     ) ) );
    ym_->correct(y_);

    // Find computational cells inside the relaxation-shape
    findComputationalCells();

    // Computate the sigma coordinate
    computeSigmaCoordinate();
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void relaxationShapePatchDist::findComputationalCells()
{
    ym_->correct(y_);
    
    const vectorField& cc = mesh_.C();

    cells_.setSize(5000);
    label count(0);

    forAll (cc, celli)
    {
        if (y_[celli]>=0. && (y_[celli] < width_))
        {
            cells_[count++] = celli;

            if (count == cells_.size())
            {
                cells_.setSize( static_cast<label>( count*1.1 ) );
            }
        }
    }

    cells_.setSize(count);
}


void relaxationShapePatchDist::computeSigmaCoordinate()
{
    ym_->correct(y_);
    
    sigma_.setSize(cells_.size(), 0);

    forAll (cells_, celli)
    {
        label ci=cells_[celli];
        sigma_[celli] = 1. - max(0., min(1., y_[ci]/width_));
        

//         if (relaxType_ == "INLET")
//         {
//             sigma_[celli] = Foam::mag( sigma_[celli] - 1.0 );
//         }
    }
}



const pointField& relaxationShapePatchDist::pointSet()
{
    notImplemented("pointSet is not implemented for this shape");
}


scalar relaxationShapePatchDist::interpolation
(
    const scalarField& source,
    const point& p0
) const
{
    notImplemented("interpolation is not implemented for this shape");
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace relaxationShapes
} // End namespace Foam

// ************************************************************************* //
