/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

#include "nonConformalProcessorCyclicFvsPatchField.H"

// * * * * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
nonConformalProcessorCyclicFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    processorCyclicFvsPatchField<Type>(p, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
nonConformalProcessorCyclicFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const Field<Type>& f
)
:
    processorCyclicFvsPatchField<Type>(p, iF, f),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
nonConformalProcessorCyclicFvsPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const dictionary& dict
)
:
    processorCyclicFvsPatchField<Type>(p, iF, dict),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
nonConformalProcessorCyclicFvsPatchField
(
    const nonConformalProcessorCyclicFvsPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, surfaceMesh>& iF,
    const fieldMapper& mapper
)
:
    processorCyclicFvsPatchField<Type>(ptf, p, iF, mapper),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(p))
{}


template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
nonConformalProcessorCyclicFvsPatchField
(
    const nonConformalProcessorCyclicFvsPatchField<Type>& ptf,
    const DimensionedField<Type, surfaceMesh>& iF
)
:
    processorCyclicFvsPatchField<Type>(ptf, iF),
    procPatch_(refCast<const nonConformalProcessorCyclicFvPatch>(ptf.patch()))
{}


// * * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * //

template<class Type>
Foam::nonConformalProcessorCyclicFvsPatchField<Type>::
~nonConformalProcessorCyclicFvsPatchField()
{}


// ************************************************************************* //
