/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "nonConformalPolyFacesFvsPatchLabelField.H"
#include "fieldMapper.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::nonConformalPolyFacesFvsPatchLabelField::
nonConformalPolyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(p, iF)
{
    if (p.size())
    {
        FatalErrorInFunction
            << "Default construction of a " << typeName << " field should "
            << "only occur on a patch with no faces"
            << exit(FatalError);
    }
}


Foam::nonConformalPolyFacesFvsPatchLabelField::
nonConformalPolyFacesFvsPatchLabelField
(
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const dictionary& dict
)
:
    fvsPatchLabelField(p, iF, dict, false)
{
    labelField::operator=
    (
        labelField("value", dict, dict.lookup<label>("size"))
    );
}


Foam::nonConformalPolyFacesFvsPatchLabelField::
nonConformalPolyFacesFvsPatchLabelField
(
    const nonConformalPolyFacesFvsPatchLabelField& ptf,
    const fvPatch& p,
    const DimensionedField<label, surfaceMesh>& iF,
    const fieldMapper& mapper
)
:
    fvsPatchLabelField(ptf, p, iF, mapper)
{
    if (p.size())
    {
        FatalErrorInFunction
            << "Mapped construction of a " << typeName << " field should "
            << "only occur on a patch with no faces"
            << exit(FatalError);
    }
}


Foam::nonConformalPolyFacesFvsPatchLabelField::
nonConformalPolyFacesFvsPatchLabelField
(
    const nonConformalPolyFacesFvsPatchLabelField& ptf,
    const DimensionedField<label, surfaceMesh>& iF
)
:
    fvsPatchLabelField(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::nonConformalPolyFacesFvsPatchLabelField::write(Ostream& os) const
{
    fvsPatchLabelField::write(os);
    writeEntry(os, "size", this->size());
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makeFvsPatchTypeField
    (
        fvsPatchLabelField,
        nonConformalPolyFacesFvsPatchLabelField
    );
}


// ************************************************************************* //
